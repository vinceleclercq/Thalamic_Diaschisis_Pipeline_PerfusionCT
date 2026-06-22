function PerfCT_Reg_MNI_001()
%PERFCT_REG_MNI_001 Convert, coregister, normalize, and summarize Tmax maps.
%
% This function:
%   1. Converts ANAT and Perf_T DICOM series to NIfTI.
%   2. Runs CTseg on the anatomical CT when required.
%   3. Rewrites Tmax as float32.
%   4. Creates a native validity mask.
%   5. Reinitializes the Tmax affine near the anatomical CT.
%   6. Coregisters and reslices Tmax and the mask to the CT grid.
%   7. Warps Tmax and the mask to CTseg template space.
%   8. Copies final subject maps to Group_results.
%   9. Creates group sum, count, mean, and coverage maps.
%
% Requirements:
%   MATLAB, SPM12, CTseg, and mu_CTseg.nii.
%
% The function reads local paths from config/local_config.m.

repo_root = fileparts(mfilename('fullpath'));
cfg = load_local_config(repo_root);
initialize_spm(cfg);

group_dir = fullfile(cfg.root_dir, 'Group_results');
if ~exist(group_dir, 'dir')
    mkdir(group_dir);
end

subjects = list_subjects(cfg);
fprintf('\nFound %d subject folders.\n', numel(subjects));

%% Step 1: DICOM conversion
conversion_rows = cell(numel(subjects), 5);

for i = 1:numel(subjects)
    subject = subjects(i).name;
    subject_dir = fullfile(cfg.root_dir, subject);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('DICOM conversion: %s (%d/%d)\n', subject, i, numel(subjects));
    fprintf('------------------------------------------------------------\n');

    anat_output = fullfile(subject_dir, [subject '_Anat.nii']);
    tmax_output = fullfile(subject_dir, [subject '_Perf_T.nii']);

    anat_status = convert_dicom_series( ...
        fullfile(subject_dir, cfg.anat_dicom_folder), ...
        anat_output, cfg.overwrite_nifti);

    tmax_status = convert_dicom_series( ...
        fullfile(subject_dir, cfg.tmax_dicom_folder), ...
        tmax_output, cfg.overwrite_nifti);

    conversion_rows(i,:) = {subject, anat_status, tmax_status, ...
                            anat_output, tmax_output};
end

conversion_table = cell2table(conversion_rows, ...
    'VariableNames', {'Subject','AnatStatus','TmaxStatus', ...
                      'AnatOutput','TmaxOutput'});

conversion_csv = fullfile(group_dir, 'dicom_conversion_summary.csv');
writetable(conversion_table, conversion_csv, ...
    'Delimiter', cfg.csv_delimiter);

fprintf('\nDICOM conversion summary:\n%s\n', conversion_csv);

%% Step 2: Registration and normalization
Subject = cell(numel(subjects),1);
Status = cell(numel(subjects),1);
AnatFile = cell(numel(subjects),1);
TmaxFile = cell(numel(subjects),1);
DeformationFile = cell(numel(subjects),1);
FinalTmaxFile = cell(numel(subjects),1);
FinalMaskFile = cell(numel(subjects),1);
RawMin = NaN(numel(subjects),1);
RawMax = NaN(numel(subjects),1);
WarpedMin = NaN(numel(subjects),1);
WarpedMax = NaN(numel(subjects),1);
WarpedNonzero = NaN(numel(subjects),1);
ValidMaskVoxels = NaN(numel(subjects),1);

for i = 1:numel(subjects)
    subject = subjects(i).name;
    subject_dir = fullfile(cfg.root_dir, subject);
    Subject{i} = subject;
    Status{i} = 'OK';

    fprintf('\n============================================================\n');
    fprintf('Image processing: %s (%d/%d)\n', subject, i, numel(subjects));
    fprintf('============================================================\n');

    final_group_tmax = fullfile(group_dir, [subject '_wfTmax.nii']);
    final_group_mask = fullfile(group_dir, [subject '_mask.nii']);

    try
        if exist(final_group_tmax, 'file') && ...
                exist(final_group_mask, 'file') && ...
                ~cfg.overwrite_processing
            fprintf('Existing final outputs retained.\n');
            Status{i} = 'SKIPPED_EXISTING';
            FinalTmaxFile{i} = final_group_tmax;
            FinalMaskFile{i} = final_group_mask;

            [WarpedMin(i), WarpedMax(i), WarpedNonzero(i)] = ...
                image_statistics(final_group_tmax);
            [~, ~, ValidMaskVoxels(i)] = ...
                image_statistics(final_group_mask);
            continue;
        end

        anat_file = fullfile(subject_dir, [subject '_Anat.nii']);
        tmax_file = fullfile(subject_dir, [subject '_Perf_T.nii']);

        if ~exist(anat_file, 'file')
            error('Anatomical NIfTI not found: %s', anat_file);
        end
        if ~exist(tmax_file, 'file')
            error('Tmax NIfTI not found: %s', tmax_file);
        end

        AnatFile{i} = anat_file;
        TmaxFile{i} = tmax_file;

        deformation_file = get_or_create_ctseg_deformation( ...
            anat_file, subject_dir);
        DeformationFile{i} = deformation_file;

        Va = spm_vol(anat_file);
        Vt = spm_vol(tmax_file);
        assert_single_volume(Va, anat_file);
        assert_single_volume(Vt, tmax_file);

        Yt = spm_read_vols(Vt);
        RawMin(i) = min(Yt(:));
        RawMax(i) = max(Yt(:));

        fprintf('Raw Tmax min/max: %g / %g\n', RawMin(i), RawMax(i));

        [tmax_path, tmax_name, tmax_ext] = fileparts(tmax_file);
        float_tmax_file = fullfile(tmax_path, ...
            ['f_' tmax_name tmax_ext]);
        native_mask_file = fullfile(tmax_path, ...
            ['mask_f_' tmax_name tmax_ext]);

        write_float32_image(Vt, Yt, float_tmax_file);
        write_validity_mask(float_tmax_file, native_mask_file, ...
            cfg.native_mask_threshold);

        reinitialize_affine_to_anatomical( ...
            anat_file, float_tmax_file, native_mask_file);

        estimate_coregistration( ...
            anat_file, float_tmax_file, native_mask_file);

        resliced_tmax = reslice_to_reference( ...
            anat_file, float_tmax_file, 1);
        resliced_mask = reslice_to_reference( ...
            anat_file, native_mask_file, 0);

        warped_tmax = warp_to_ctseg_space( ...
            deformation_file, cfg.mu_file, ...
            resliced_tmax, subject_dir, 1);

        warped_mask = warp_to_ctseg_space( ...
            deformation_file, cfg.mu_file, ...
            resliced_mask, subject_dir, 0);

        subject_mask = fullfile(subject_dir, [subject '_mask.nii']);
        binarize_mask(warped_mask, subject_mask);

        check_same_geometry(warped_tmax, subject_mask);

        [WarpedMin(i), WarpedMax(i), WarpedNonzero(i)] = ...
            image_statistics(warped_tmax);
        [~, ~, ValidMaskVoxels(i)] = image_statistics(subject_mask);

        copyfile(warped_tmax, final_group_tmax, 'f');
        copyfile(subject_mask, final_group_mask, 'f');

        FinalTmaxFile{i} = final_group_tmax;
        FinalMaskFile{i} = final_group_mask;

        fprintf('Final normalized Tmax: %s\n', final_group_tmax);
        fprintf('Final validity mask  : %s\n', final_group_mask);

    catch ME
        Status{i} = ['ERROR: ' ME.message];
        fprintf(2, 'ERROR for %s: %s\n', subject, ME.message);
    end
end

processing_table = table(Subject, Status, AnatFile, TmaxFile, ...
    DeformationFile, FinalTmaxFile, FinalMaskFile, RawMin, RawMax, ...
    WarpedMin, WarpedMax, WarpedNonzero, ValidMaskVoxels);

processing_csv = fullfile(group_dir, 'processing_summary.csv');
writetable(processing_table, processing_csv, ...
    'Delimiter', cfg.csv_delimiter);

fprintf('\nProcessing summary:\n%s\n', processing_csv);

create_group_maps(group_dir, Status, Subject, ...
    FinalTmaxFile, FinalMaskFile);

fprintf('\nPerfCT_Reg_MNI_001 completed.\n');
end

%% Local helper functions

function cfg = load_local_config(repo_root)
config_file = fullfile(repo_root, 'config', 'local_config.m');
example_file = fullfile(repo_root, 'config', 'example_config.m');

if ~exist(config_file, 'file')
    error(['Missing configuration file:\n%s\n\n' ...
           'Copy %s to local_config.m and edit the paths.'], ...
           config_file, example_file);
end

cfg = struct(); %#ok<NASGU>
run(config_file);

required_fields = {'root_dir','spm_dir','ctseg_dir','fieldtrip_dir'};
for i = 1:numel(required_fields)
    field_name = required_fields{i};
    if ~isfield(cfg, field_name) || isempty(cfg.(field_name))
        error('Missing configuration field: cfg.%s', field_name);
    end
end

if ~isfield(cfg, 'subject_pattern')
    cfg.subject_pattern = 'Sub*';
end
if ~isfield(cfg, 'anat_dicom_folder')
    cfg.anat_dicom_folder = 'ANAT';
end
if ~isfield(cfg, 'tmax_dicom_folder')
    cfg.tmax_dicom_folder = 'Perf_T';
end
if ~isfield(cfg, 'overwrite_nifti')
    cfg.overwrite_nifti = false;
end
if ~isfield(cfg, 'overwrite_processing')
    cfg.overwrite_processing = false;
end
if ~isfield(cfg, 'native_mask_threshold')
    cfg.native_mask_threshold = -100;
end
if ~isfield(cfg, 'csv_delimiter')
    cfg.csv_delimiter = ';';
end
if ~isfield(cfg, 'mu_file') || isempty(cfg.mu_file)
    cfg.mu_file = fullfile(cfg.ctseg_dir, 'mu_CTseg.nii');
end

if ~exist(cfg.root_dir, 'dir')
    error('Root data directory not found: %s', cfg.root_dir);
end
if ~exist(cfg.spm_dir, 'dir')
    error('SPM directory not found: %s', cfg.spm_dir);
end
if ~exist(cfg.ctseg_dir, 'dir')
    error('CTseg directory not found: %s', cfg.ctseg_dir);
end
if ~exist(cfg.mu_file, 'file')
    error('mu_CTseg.nii not found: %s', cfg.mu_file);
end
end

function initialize_spm(cfg)
addpath(cfg.spm_dir);
addpath(cfg.ctseg_dir);

if exist('spm', 'file') ~= 2
    error('SPM was not found after adding %s', cfg.spm_dir);
end
if exist('spm_CTseg', 'file') ~= 2
    error('spm_CTseg.m was not found after adding %s', cfg.ctseg_dir);
end

spm('defaults', 'fmri');
spm_jobman('initcfg');
end

function subjects = list_subjects(cfg)
subjects = dir(fullfile(cfg.root_dir, cfg.subject_pattern));
subjects = subjects([subjects.isdir]);
subjects = subjects(~ismember({subjects.name}, {'.','..'}));

if isempty(subjects)
    error('No subject folders found in %s', cfg.root_dir);
end

[~, order] = sort({subjects.name});
subjects = subjects(order);
end

function status = convert_dicom_series(dicom_dir, output_file, overwrite)
if ~exist(dicom_dir, 'dir')
    status = 'MISSING_DICOM_FOLDER';
    fprintf(2, 'Missing DICOM folder: %s\n', dicom_dir);
    return;
end

if exist(output_file, 'file') && ~overwrite
    status = 'SKIPPED_EXISTING';
    fprintf('Existing NIfTI retained: %s\n', output_file);
    return;
end

file_list = spm_select('FPList', dicom_dir, '.*');
if isempty(file_list)
    status = 'NO_DICOM_FILES';
    fprintf(2, 'No files found in: %s\n', dicom_dir);
    return;
end

temporary_dir = [tempname(fileparts(output_file)) '_dicom'];
mkdir(temporary_dir);
cleanup_object = onCleanup(@() remove_directory(temporary_dir)); %#ok<NASGU>

current_dir = pwd;
cleanup_directory = onCleanup(@() cd(current_dir)); %#ok<NASGU>
cd(temporary_dir);

try
    headers = spm_dicom_headers(file_list);
    spm_dicom_convert(headers, 'all', 'flat', 'nii');

    nifti_files = dir(fullfile(temporary_dir, '*.nii'));
    if isempty(nifti_files)
        status = 'NO_NIFTI_CREATED';
        fprintf(2, 'No NIfTI created from %s\n', dicom_dir);
        return;
    end

    [~, largest_index] = max([nifti_files.bytes]);
    source_file = fullfile(temporary_dir, ...
        nifti_files(largest_index).name);

    if exist(output_file, 'file')
        delete(output_file);
    end
    copyfile(source_file, output_file, 'f');

    status = 'OK';
    fprintf('Created: %s\n', output_file);

catch ME
    status = ['ERROR: ' ME.message];
    fprintf(2, 'DICOM conversion error: %s\n', ME.message);
end
end

function remove_directory(directory_name)
if exist(directory_name, 'dir')
    try
        rmdir(directory_name, 's');
    catch
        warning('Could not remove temporary folder: %s', directory_name);
    end
end
end

function deformation_file = get_or_create_ctseg_deformation( ...
    anat_file, subject_dir)

deformation_candidates = spm_select( ...
    'FPList', subject_dir, '^y_.*\.nii$');

if isempty(deformation_candidates)
    fprintf('No deformation field found. Running CTseg...\n');
    spm_CTseg(anat_file, subject_dir, ...
        true, true, true, true, 1.0);

    deformation_candidates = spm_select( ...
        'FPList', subject_dir, '^y_.*\.nii$');

    if isempty(deformation_candidates)
        error('CTseg did not create a y_*.nii deformation field.');
    end
else
    fprintf('Existing CTseg deformation field retained.\n');
end

deformation_file = deblank(deformation_candidates(1,:));
fprintf('Deformation field: %s\n', deformation_file);
end

function assert_single_volume(V, file_name)
if numel(V) ~= 1
    error('Expected one 3D volume in %s, found %d.', ...
          file_name, numel(V));
end
end

function write_float32_image(V, data, output_file)
Vout = V;
Vout.fname = output_file;
Vout.dt = [16 0];
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, single(data));
end

function write_validity_mask(float_file, mask_file, threshold)
V = spm_vol(float_file);
Y = spm_read_vols(V);
valid = isfinite(Y) & Y > threshold;

Vout = V;
Vout.fname = mask_file;
Vout.dt = [2 0];
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, uint8(valid));

fprintf('Native validity mask: %s (%d voxels)\n', ...
        mask_file, nnz(valid));
end

function reinitialize_affine_to_anatomical( ...
    anat_file, tmax_file, mask_file)

Va = spm_vol(anat_file);
Vt = spm_vol(tmax_file);

voxel_size_tmax = sqrt(sum(Vt.mat(1:3,1:3).^2, 1));
rotation_anat = Va.mat(1:3,1:3);
rotation_norm = sqrt(sum(rotation_anat.^2, 1));

if any(rotation_norm == 0) || any(voxel_size_tmax == 0)
    error('Invalid affine matrix encountered.');
end

rotation_unit = bsxfun(@rdivide, ...
    rotation_anat, rotation_norm);
rotation_new = bsxfun(@times, ...
    rotation_unit, voxel_size_tmax);

centre_anat_voxel = [(Va.dim(1)+1)/2; ...
                     (Va.dim(2)+1)/2; ...
                     (Va.dim(3)+1)/2; 1];

centre_tmax_voxel = [(Vt.dim(1)+1)/2; ...
                     (Vt.dim(2)+1)/2; ...
                     (Vt.dim(3)+1)/2; 1];

centre_anat_mm = Va.mat * centre_anat_voxel;

new_affine = eye(4);
new_affine(1:3,1:3) = rotation_new;
new_affine(1:3,4) = centre_anat_mm(1:3) ...
    - rotation_new * centre_tmax_voxel(1:3);

spm_get_space(tmax_file, new_affine);
spm_get_space(mask_file, new_affine);

fprintf('Manual affine reinitialization applied.\n');
end

function estimate_coregistration(anat_file, tmax_file, mask_file)
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estimate.ref = ...
    {[anat_file ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source = ...
    {[tmax_file ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = ...
    {[mask_file ',1']};

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

function output_file = reslice_to_reference( ...
    reference_file, source_file, interpolation)

matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.write.ref = ...
    {[reference_file ',1']};
matlabbatch{1}.spm.spatial.coreg.write.source = ...
    {[source_file ',1']};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = ...
    interpolation;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);

[source_path, source_name, source_extension] = ...
    fileparts(source_file);
output_file = fullfile(source_path, ...
    ['r' source_name source_extension]);

if ~exist(output_file, 'file')
    error('Resliced file was not created: %s', output_file);
end
end

function output_file = warp_to_ctseg_space( ...
    deformation_file, mu_file, input_file, output_dir, interpolation)

matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = ...
    {deformation_file};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = ...
    {mu_file};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = ...
    {input_file};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = ...
    {output_dir};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = ...
    interpolation;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';

spm_jobman('run', matlabbatch);

[~, input_name, input_extension] = fileparts(input_file);
output_file = fullfile(output_dir, ...
    ['w' input_name input_extension]);

if ~exist(output_file, 'file')
    error('Warped file was not created: %s', output_file);
end
end

function binarize_mask(input_mask, output_mask)
V = spm_vol(input_mask);
Y = spm_read_vols(V);
binary_mask = isfinite(Y) & Y > 0.5;

Vout = V;
Vout.fname = output_mask;
Vout.dt = [2 0];
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, uint8(binary_mask));

fprintf('Final mask: %s (%d voxels)\n', ...
        output_mask, nnz(binary_mask));
end

function check_same_geometry(file_a, file_b)
Va = spm_vol(file_a);
Vb = spm_vol(file_b);
assert_single_volume(Va, file_a);
assert_single_volume(Vb, file_b);

if any(Va.dim ~= Vb.dim)
    error('Image dimensions do not match.');
end
if max(abs(Va.mat(:) - Vb.mat(:))) > 1e-4
    error('Image affine matrices do not match.');
end
end

function [minimum_value, maximum_value, nonzero_count] = ...
    image_statistics(file_name)

V = spm_vol(file_name);
assert_single_volume(V, file_name);
Y = spm_read_vols(V);

finite_values = Y(isfinite(Y));
if isempty(finite_values)
    minimum_value = NaN;
    maximum_value = NaN;
else
    minimum_value = min(finite_values);
    maximum_value = max(finite_values);
end

nonzero_count = nnz(isfinite(Y) & Y ~= 0);
end

function create_group_maps(group_dir, Status, Subject, ...
    FinalTmaxFile, FinalMaskFile)

fprintf('\n============================================================\n');
fprintf('Creating group maps\n');
fprintf('============================================================\n');

ok = strcmp(Status, 'OK') | strcmp(Status, 'SKIPPED_EXISTING');
ok_tmax = FinalTmaxFile(ok);
ok_mask = FinalMaskFile(ok);
ok_subject = Subject(ok);

valid_entries = ~cellfun(@isempty, ok_tmax) & ...
                ~cellfun(@isempty, ok_mask);
ok_tmax = ok_tmax(valid_entries);
ok_mask = ok_mask(valid_entries);
ok_subject = ok_subject(valid_entries);

if isempty(ok_tmax)
    warning('No successfully processed subjects available.');
    return;
end

Vreference = spm_vol(ok_tmax{1});
Yreference = spm_read_vols(Vreference);

sum_map = zeros(size(Yreference), 'double');
count_map = zeros(size(Yreference), 'double');
included_subjects = {};

for i = 1:numel(ok_tmax)
    if ~exist(ok_tmax{i}, 'file') || ~exist(ok_mask{i}, 'file')
        warning('Missing final output for %s; excluded.', ...
                ok_subject{i});
        continue;
    end

    Vt = spm_vol(ok_tmax{i});
    Vm = spm_vol(ok_mask{i});

    if any(Vt.dim ~= Vreference.dim) || ...
            max(abs(Vt.mat(:) - Vreference.mat(:))) > 1e-4
        warning('Geometry mismatch for %s; excluded.', ok_subject{i});
        continue;
    end
    if any(Vm.dim ~= Vt.dim) || ...
            max(abs(Vm.mat(:) - Vt.mat(:))) > 1e-4
        warning('Mask mismatch for %s; excluded.', ok_subject{i});
        continue;
    end

    Yt = spm_read_vols(Vt);
    Ym = spm_read_vols(Vm);

    valid = isfinite(Yt) & Ym > 0.5;
    sum_map(valid) = sum_map(valid) + Yt(valid);
    count_map(valid) = count_map(valid) + 1;
    included_subjects{end+1,1} = ok_subject{i}; %#ok<AGROW>
end

if isempty(included_subjects)
    warning('No subjects passed the group-map checks.');
    return;
end

mean_map = zeros(size(sum_map), 'double');
contributing = count_map > 0;
mean_map(contributing) = ...
    sum_map(contributing) ./ count_map(contributing);

coverage_map = count_map / numel(included_subjects);

write_group_image(Vreference, sum_map, ...
    fullfile(group_dir, 'group_sum_wfTmax.nii'));
write_group_image(Vreference, count_map, ...
    fullfile(group_dir, 'group_count_valid.nii'));
write_group_image(Vreference, mean_map, ...
    fullfile(group_dir, 'group_mean_wfTmax.nii'));
write_group_image(Vreference, coverage_map, ...
    fullfile(group_dir, 'group_coverage_fraction.nii'));

subjects_file = fullfile(group_dir, 'group_subjects.txt');
file_id = fopen(subjects_file, 'w');
if file_id < 0
    error('Could not create %s', subjects_file);
end
cleanup_file = onCleanup(@() fclose(file_id)); %#ok<NASGU>
fprintf(file_id, '%s\n', included_subjects{:});

fprintf('Group maps created from %d subjects.\n', ...
        numel(included_subjects));
end

function write_group_image(reference_volume, data, output_file)
Vout = reference_volume;
Vout.fname = output_file;
Vout.dt = [16 0];
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, single(data));
end
