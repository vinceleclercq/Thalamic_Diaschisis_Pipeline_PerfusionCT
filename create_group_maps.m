function summary_table = register_tmax_to_template(cfg)
%REGISTER_TMAX_TO_TEMPLATE Coregister and normalize Tmax maps.
%
% The workflow reproduces the original processing strategy:
%   1. find anatomical CT and Tmax NIfTI files;
%   2. run CTseg if no deformation field exists;
%   3. rewrite Tmax as float32;
%   4. create a native validity mask;
%   5. reinitialize the Tmax affine near the anatomical CT;
%   6. estimate rigid coregistration to CT;
%   7. reslice Tmax and mask to the CT grid;
%   8. warp both to CTseg template space;
%   9. copy final maps to Group_results;
%  10. create group maps.

fprintf('\n============================================================\n');
fprintf('STEP 2 - TMAX REGISTRATION AND NORMALIZATION\n');
fprintf('============================================================\n');

subjects = list_subject_directories(cfg);
n = numel(subjects);

Subject = cell(n,1);
Status = cell(n,1);
AnatFile = cell(n,1);
TmaxFile = cell(n,1);
DeformationFile = cell(n,1);
FloatFile = cell(n,1);
ReslicedFile = cell(n,1);
WarpedFile = cell(n,1);
MaskFile = cell(n,1);
GroupWarpedFile = cell(n,1);
GroupMaskFile = cell(n,1);
RawMin = NaN(n,1);
RawMax = NaN(n,1);
WarpedMin = NaN(n,1);
WarpedMax = NaN(n,1);
WarpedNonzero = NaN(n,1);
ValidMaskVoxels = NaN(n,1);

for i = 1:n
    subject = subjects(i).name;
    subject_dir = fullfile(cfg.root_dir, subject);

    Subject{i} = subject;
    Status{i} = 'OK';

    group_warped = fullfile(cfg.group_dir, [subject '_wfTmax.nii']);
    group_mask = fullfile(cfg.group_dir, [subject '_mask.nii']);
    GroupWarpedFile{i} = group_warped;
    GroupMaskFile{i} = group_mask;

    fprintf('\n[%d/%d] %s\n', i, n, subject);

    try
        if exist(group_warped, 'file') && exist(group_mask, 'file') ...
                && ~cfg.overwrite_processing
            fprintf('  Existing normalized outputs retained.\n');
            Status{i} = 'SKIPPED_EXISTING';

            [WarpedMin(i), WarpedMax(i), WarpedNonzero(i)] = ...
                image_qc(group_warped);
            [~, ~, ValidMaskVoxels(i)] = image_qc(group_mask);
            WarpedFile{i} = group_warped;
            MaskFile{i} = group_mask;
            continue;
        end

        anat_file = find_subject_file(subject_dir, subject, 'anat');
        tmax_file = find_subject_file(subject_dir, subject, 'tmax');

        AnatFile{i} = anat_file;
        TmaxFile{i} = tmax_file;

        fprintf('  Anatomical CT: %s\n', anat_file);
        fprintf('  Tmax map     : %s\n', tmax_file);

        deformation_file = get_or_create_deformation( ...
            subject_dir, anat_file);
        DeformationFile{i} = deformation_file;

        Va = spm_vol(anat_file);
        Vt = spm_vol(tmax_file);
        assert_single_volume(Va, anat_file);
        assert_single_volume(Vt, tmax_file);

        Yt = spm_read_vols(Vt);
        RawMin(i) = min(Yt(:), [], 'omitnan');
        RawMax(i) = max(Yt(:), [], 'omitnan');

        [tmax_path, tmax_name, tmax_ext] = fileparts(tmax_file);
        float_file = fullfile(tmax_path, ['f_' tmax_name tmax_ext]);
        native_mask_file = fullfile(tmax_path, ...
            ['mask_f_' tmax_name tmax_ext]);

        FloatFile{i} = float_file;

        write_float_tmax(Vt, Yt, float_file);
        write_native_mask(float_file, native_mask_file, ...
            cfg.native_mask_lower_bound);

        reinitialize_affine(anat_file, float_file, native_mask_file);
        estimate_coregistration(anat_file, float_file, native_mask_file);

        resliced_file = reslice_to_reference( ...
            anat_file, float_file, 1, 'r');
        resliced_mask = reslice_to_reference( ...
            anat_file, native_mask_file, 0, 'r');

        ReslicedFile{i} = resliced_file;

        warped_file = warp_to_ctseg_template( ...
            deformation_file, cfg.mu_file, resliced_file, 1, subject_dir);
        warped_mask_raw = warp_to_ctseg_template( ...
            deformation_file, cfg.mu_file, resliced_mask, 0, subject_dir);

        WarpedFile{i} = warped_file;

        final_mask_file = fullfile(subject_dir, [subject '_mask.nii']);
        create_final_binary_mask(warped_mask_raw, final_mask_file);
        MaskFile{i} = final_mask_file;

        [WarpedMin(i), WarpedMax(i), WarpedNonzero(i)] = ...
            image_qc(warped_file);
        [~, ~, ValidMaskVoxels(i)] = image_qc(final_mask_file);

        check_matching_geometry(warped_file, final_mask_file, ...
            cfg.affine_tolerance);

        copyfile(warped_file, group_warped, 'f');
        copyfile(final_mask_file, group_mask, 'f');

        fprintf('  Final Tmax: %s\n', group_warped);
        fprintf('  Final mask: %s\n', group_mask);

        if cfg.cleanup_intermediate_files
            delete_if_exists(float_file);
            delete_if_exists(native_mask_file);
            delete_if_exists(resliced_file);
            delete_if_exists(resliced_mask);
            delete_if_exists(warped_mask_raw);
        end

    catch ME
        Status{i} = ['ERROR: ' ME.message];
        fprintf(2, '  ERROR: %s\n', ME.message);
    end
end

summary_table = table(Subject, Status, AnatFile, TmaxFile, ...
    DeformationFile, FloatFile, ReslicedFile, WarpedFile, MaskFile, ...
    GroupWarpedFile, GroupMaskFile, RawMin, RawMax, WarpedMin, ...
    WarpedMax, WarpedNonzero, ValidMaskVoxels);

summary_file = fullfile(cfg.group_dir, 'processing_summary.csv');
writetable(summary_table, summary_file, 'Delimiter', cfg.csv_delimiter);
fprintf('\nProcessing summary written to:\n%s\n', summary_file);

if cfg.create_group_maps
    create_group_maps(cfg);
end
end

function deformation_file = get_or_create_deformation(subject_dir, anat_file)
candidates = spm_select('FPList', subject_dir, '^y_.*\.nii$');

if isempty(candidates)
    fprintf('  Running CTseg...\n');
    [~, ~] = spm_CTseg(anat_file, subject_dir, ...
        true, true, true, true, 1.0);

    candidates = spm_select('FPList', subject_dir, '^y_.*\.nii$');
    if isempty(candidates)
        error('CTseg completed but no y_*.nii deformation field was found.');
    end
else
    fprintf('  Existing CTseg deformation field retained.\n');
end

deformation_file = deblank(candidates(1,:));
end

function write_float_tmax(Vt, Yt, output_file)
Vout = Vt;
Vout.fname = output_file;
Vout.dt = [16 0];     % float32
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, single(Yt));
end

function write_native_mask(float_file, mask_file, lower_bound)
V = spm_vol(float_file);
Y = spm_read_vols(V);
valid = isfinite(Y) & Y > lower_bound;

Vmask = V;
Vmask.fname = mask_file;
Vmask.dt = [2 0];     % uint8
Vmask.pinfo = [1; 0; 0];
spm_write_vol(Vmask, uint8(valid));
end

function reinitialize_affine(anat_file, float_file, mask_file)
Va = spm_vol(anat_file);
Vf = spm_vol(float_file);

voxel_size = sqrt(sum(Vf.mat(1:3,1:3).^2, 1));
rotation = Va.mat(1:3,1:3);
column_norms = sqrt(sum(rotation.^2, 1));

if any(column_norms == 0) || any(voxel_size == 0)
    error('Invalid affine matrix encountered during reinitialization.');
end

rotation_unit = bsxfun(@rdivide, rotation, column_norms);
rotation_new = bsxfun(@times, rotation_unit, voxel_size);

centre_anat_vox = [(Va.dim(1)+1)/2; (Va.dim(2)+1)/2; ...
                   (Va.dim(3)+1)/2; 1];
centre_tmax_vox = [(Vf.dim(1)+1)/2; (Vf.dim(2)+1)/2; ...
                   (Vf.dim(3)+1)/2; 1];
centre_anat_mm = Va.mat * centre_anat_vox;

new_affine = eye(4);
new_affine(1:3,1:3) = rotation_new;
new_affine(1:3,4) = centre_anat_mm(1:3) ...
    - rotation_new * centre_tmax_vox(1:3);

spm_get_space(float_file, new_affine);
spm_get_space(mask_file, new_affine);
end

function estimate_coregistration(anat_file, source_file, other_file)
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estimate.ref = ...
    {[anat_file ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source = ...
    {[source_file ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = ...
    {[other_file ',1']};

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
    [0.0200 0.0200 0.0200 ...
     0.0010 0.0010 0.0010 ...
     0.0100 0.0100 0.0100 ...
     0.0010 0.0010 0.0010];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

spm_jobman('run', matlabbatch);
end

function output_file = reslice_to_reference(ref_file, source_file, interp, prefix)
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.write.ref = {[ref_file ',1']};
matlabbatch{1}.spm.spatial.coreg.write.source = {[source_file ',1']};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = prefix;
spm_jobman('run', matlabbatch);

[path_name, file_name, extension] = fileparts(source_file);
output_file = fullfile(path_name, [prefix file_name extension]);

if ~exist(output_file, 'file')
    error('Expected resliced file was not created: %s', output_file);
end
end

function output_file = warp_to_ctseg_template( ...
    deformation_file, template_file, input_file, interp, output_dir)

matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = ...
    {deformation_file};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = ...
    {template_file};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = ...
    {input_file};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = ...
    {output_dir};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = interp;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
spm_jobman('run', matlabbatch);

[~, file_name, extension] = fileparts(input_file);
output_file = fullfile(output_dir, ['w' file_name extension]);

if ~exist(output_file, 'file')
    error('Expected warped file was not created: %s', output_file);
end
end

function create_final_binary_mask(input_mask, output_mask)
V = spm_vol(input_mask);
Y = spm_read_vols(V);
binary_mask = isfinite(Y) & Y > 0.5;

Vout = V;
Vout.fname = output_mask;
Vout.dt = [2 0];
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, uint8(binary_mask));
end

function [minimum, maximum, nonzero] = image_qc(file_name)
V = spm_vol(file_name);
assert_single_volume(V, file_name);
Y = spm_read_vols(V);
minimum = min(Y(:), [], 'omitnan');
maximum = max(Y(:), [], 'omitnan');
nonzero = nnz(isfinite(Y) & Y ~= 0);
end

function check_matching_geometry(file_a, file_b, tolerance)
Va = spm_vol(file_a);
Vb = spm_vol(file_b);
assert_single_volume(Va, file_a);
assert_single_volume(Vb, file_b);

if any(Va.dim ~= Vb.dim)
    error('Image dimensions do not match:\n%s\n%s', file_a, file_b);
end
if max(abs(Va.mat(:) - Vb.mat(:))) > tolerance
    error('Image affine matrices do not match:\n%s\n%s', file_a, file_b);
end
end

function assert_single_volume(V, file_name)
if numel(V) ~= 1
    error('Expected a single 3D volume but found %d in %s.', ...
          numel(V), file_name);
end
end

function file_name = find_subject_file(subject_dir, subject, kind)
switch lower(kind)
    case 'anat'
        patterns = { ...
            ['^' regexptranslate('escape', subject) '.*_Anat\.nii$'], ...
            '.*_Anat\.nii$', ...
            '.*CT.*\.nii$', ...
            '.*anat.*\.nii$'};
    case 'tmax'
        patterns = { ...
            ['^' regexptranslate('escape', subject) '.*_Perf_T\.nii$'], ...
            '.*_Perf_T\.nii$', ...
            '.*Tmax.*\.nii$', ...
            '.*Perf_T.*\.nii$', ...
            '.*tmax.*\.nii$'};
    otherwise
        error('Unknown file kind: %s', kind);
end

file_name = '';
for p = 1:numel(patterns)
    candidate = spm_select('FPList', subject_dir, patterns{p});
    if ~isempty(candidate)
        file_name = deblank(candidate(1,:));
        return;
    end
end

error('No %s NIfTI file found in %s.', kind, subject_dir);
end

function delete_if_exists(file_name)
if exist(file_name, 'file')
    delete(file_name);
end
end

function subjects = list_subject_directories(cfg)
subjects = dir(fullfile(cfg.root_dir, cfg.subject_pattern));
subjects = subjects([subjects.isdir]);
subjects = subjects(~ismember({subjects.name}, {'.', '..'}));

if isempty(subjects)
    error('No subject directories matching "%s" were found in %s.', ...
          cfg.subject_pattern, cfg.root_dir);
end

[~, order] = sort({subjects.name});
subjects = subjects(order);
end
