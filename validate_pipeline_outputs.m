function output_table = extract_thalamic_tmax(cfg)
%EXTRACT_THALAMIC_TMAX Extract bilateral thalamic Tmax measurements.
%
% Inputs in Group_results:
%   <Subject>_wfTmax.nii
%   <Subject>_mask.nii
%
% Outputs:
%   thalamus_CT_values.csv
%   thalamus_CT_results.mat

fprintf('\n============================================================\n');
fprintf('STEP 3 - BILATERAL THALAMIC TMAX EXTRACTION\n');
fprintf('============================================================\n');

atlas_file = resolve_atlas_file(cfg);
atlas = ft_read_atlas(atlas_file);

left_id = find_exact_atlas_label(atlas, cfg.aal_left_label);
right_id = find_exact_atlas_label(atlas, cfg.aal_right_label);

fprintf('AAL atlas: %s\n', atlas_file);
fprintf('Left thalamus label : %s (ID %d)\n', ...
    cfg.aal_left_label, left_id);
fprintf('Right thalamus label: %s (ID %d)\n', ...
    cfg.aal_right_label, right_id);

maps = dir(fullfile(cfg.group_dir, '*_wfTmax.nii'));
if isempty(maps)
    error('No normalized Tmax maps found in %s.', cfg.group_dir);
end
[~, order] = sort({maps.name});
maps = maps(order);
n = numel(maps);

Subject = cell(n,1);
TmaxFile = cell(n,1);
MaskFile = cell(n,1);
VoxelSizeX_mm = NaN(n,1);
VoxelSizeY_mm = NaN(n,1);
VoxelSizeZ_mm = NaN(n,1);
Thal_L_Mean = NaN(n,1);
Thal_L_Median = NaN(n,1);
Thal_L_Std = NaN(n,1);
Thal_L_nAtlasVoxels = NaN(n,1);
Thal_L_nValidVoxels = NaN(n,1);
Thal_L_Coverage = NaN(n,1);
Thal_R_Mean = NaN(n,1);
Thal_R_Median = NaN(n,1);
Thal_R_Std = NaN(n,1);
Thal_R_nAtlasVoxels = NaN(n,1);
Thal_R_nValidVoxels = NaN(n,1);
Thal_R_Coverage = NaN(n,1);
Status = cell(n,1);

left_values = cell(n,1);
right_values = cell(n,1);

cached_dim = [];
cached_mat = [];
cached_left_mask = [];
cached_right_mask = [];

for i = 1:n
    subject = erase(maps(i).name, '_wfTmax.nii');
    tmax_file = fullfile(cfg.group_dir, maps(i).name);
    mask_file = fullfile(cfg.group_dir, [subject '_mask.nii']);

    Subject{i} = subject;
    TmaxFile{i} = tmax_file;
    MaskFile{i} = mask_file;
    Status{i} = 'OK';

    fprintf('\n[%d/%d] %s\n', i, n, subject);

    try
        if ~exist(mask_file, 'file')
            error('Validity mask not found: %s', mask_file);
        end

        Vt = spm_vol(tmax_file);
        Vm = spm_vol(mask_file);

        if numel(Vt) ~= 1 || numel(Vm) ~= 1
            error('Tmax and mask files must each contain one 3D volume.');
        end
        if any(Vt.dim ~= Vm.dim) || ...
                max(abs(Vt.mat(:) - Vm.mat(:))) > cfg.affine_tolerance
            error('Tmax and mask geometries do not match.');
        end

        voxel_size = sqrt(sum(Vt.mat(1:3,1:3).^2, 1));
        VoxelSizeX_mm(i) = voxel_size(1);
        VoxelSizeY_mm(i) = voxel_size(2);
        VoxelSizeZ_mm(i) = voxel_size(3);

        if isempty(cached_dim) || any(Vt.dim ~= cached_dim) || ...
                max(abs(Vt.mat(:) - cached_mat(:))) ...
                > cfg.affine_tolerance
            fprintf('  Resampling AAL thalamic labels to target grid...\n');
            [cached_left_mask, cached_right_mask] = ...
                resample_thalamic_masks_nearest( ...
                    atlas, left_id, right_id, Vt);
            cached_dim = Vt.dim;
            cached_mat = Vt.mat;
        end

        Yt = double(spm_read_vols(Vt));
        Ym = double(spm_read_vols(Vm));

        valid_data = isfinite(Yt) & Ym > 0.5;
        if cfg.extract_positive_tmax_only
            valid_data = valid_data & Yt > 0;
        end

        left_atlas_n = nnz(cached_left_mask);
        right_atlas_n = nnz(cached_right_mask);

        left_valid_mask = cached_left_mask & valid_data;
        right_valid_mask = cached_right_mask & valid_data;

        lv = Yt(left_valid_mask);
        rv = Yt(right_valid_mask);

        if isempty(lv)
            warning('No valid left-thalamic voxels for %s.', subject);
        else
            Thal_L_Mean(i) = mean(lv);
            Thal_L_Median(i) = median(lv);
            Thal_L_Std(i) = std(lv);
            left_values{i} = lv;
        end

        if isempty(rv)
            warning('No valid right-thalamic voxels for %s.', subject);
        else
            Thal_R_Mean(i) = mean(rv);
            Thal_R_Median(i) = median(rv);
            Thal_R_Std(i) = std(rv);
            right_values{i} = rv;
        end

        Thal_L_nAtlasVoxels(i) = left_atlas_n;
        Thal_L_nValidVoxels(i) = numel(lv);
        Thal_L_Coverage(i) = safe_fraction(numel(lv), left_atlas_n);

        Thal_R_nAtlasVoxels(i) = right_atlas_n;
        Thal_R_nValidVoxels(i) = numel(rv);
        Thal_R_Coverage(i) = safe_fraction(numel(rv), right_atlas_n);

        fprintf('  Left : mean %.4f, valid voxels %d/%d\n', ...
            Thal_L_Mean(i), numel(lv), left_atlas_n);
        fprintf('  Right: mean %.4f, valid voxels %d/%d\n', ...
            Thal_R_Mean(i), numel(rv), right_atlas_n);

    catch ME
        Status{i} = ['ERROR: ' ME.message];
        fprintf(2, '  ERROR: %s\n', ME.message);
    end
end

output_table = table(Subject, Status, TmaxFile, MaskFile, ...
    VoxelSizeX_mm, VoxelSizeY_mm, VoxelSizeZ_mm, ...
    Thal_L_Mean, Thal_L_Median, Thal_L_Std, ...
    Thal_L_nAtlasVoxels, Thal_L_nValidVoxels, Thal_L_Coverage, ...
    Thal_R_Mean, Thal_R_Median, Thal_R_Std, ...
    Thal_R_nAtlasVoxels, Thal_R_nValidVoxels, Thal_R_Coverage);

csv_file = fullfile(cfg.group_dir, 'thalamus_CT_values.csv');
writetable(output_table, csv_file, 'Delimiter', cfg.csv_delimiter);
fprintf('\nThalamic summary written to:\n%s\n', csv_file);

if cfg.save_voxel_values
    mat_file = fullfile(cfg.group_dir, 'thalamus_CT_results.mat');
    metadata = struct();
    metadata.atlas_file = atlas_file;
    metadata.left_label = cfg.aal_left_label;
    metadata.right_label = cfg.aal_right_label;
    metadata.extract_positive_tmax_only = ...
        cfg.extract_positive_tmax_only;
    metadata.created = datestr(now, 30); %#ok<DATST>
    save(mat_file, 'output_table', 'left_values', 'right_values', ...
        'metadata', '-v7.3');
    fprintf('Voxel-level local MAT output written to:\n%s\n', mat_file);
end
end

function atlas_file = resolve_atlas_file(cfg)
if ~isempty(cfg.aal_atlas_file)
    atlas_file = cfg.aal_atlas_file;
else
    atlas_file = fullfile(fileparts(which('ft_defaults')), ...
        'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');
end

if ~exist(atlas_file, 'file')
    error('AAL atlas not found: %s', atlas_file);
end
end

function label_id = find_exact_atlas_label(atlas, label_name)
matches = find(strcmpi(atlas.tissuelabel, label_name));
if numel(matches) ~= 1
    available = strjoin(atlas.tissuelabel( ...
        contains(atlas.tissuelabel, 'Thalamus', ...
        'IgnoreCase', true)), ', ');
    error(['Expected one exact atlas label "%s" but found %d. ' ...
           'Available thalamic labels: %s'], ...
          label_name, numel(matches), available);
end
label_id = matches(1);
end

function [left_mask, right_mask] = resample_thalamic_masks_nearest( ...
    atlas, left_id, right_id, Vtarget)

target_dim = Vtarget.dim;
left_mask = false(target_dim);
right_mask = false(target_dim);

atlas_inverse = inv(atlas.transform);
[x_grid, y_grid] = ndgrid(1:target_dim(1), 1:target_dim(2));
n_xy = numel(x_grid);

for z = 1:target_dim(3)
    target_voxels = [ ...
        x_grid(:)'; ...
        y_grid(:)'; ...
        repmat(z, 1, n_xy); ...
        ones(1, n_xy)];

    world_mm = Vtarget.mat * target_voxels;
    atlas_voxels = atlas_inverse * world_mm;

    ax = round(atlas_voxels(1,:));
    ay = round(atlas_voxels(2,:));
    az = round(atlas_voxels(3,:));

    in_bounds = ax >= 1 & ax <= atlas.dim(1) & ...
                ay >= 1 & ay <= atlas.dim(2) & ...
                az >= 1 & az <= atlas.dim(3);

    labels = zeros(1, n_xy);
    linear_indices = sub2ind(atlas.dim, ...
        ax(in_bounds), ay(in_bounds), az(in_bounds));
    labels(in_bounds) = double(atlas.tissue(linear_indices));

    left_mask(:,:,z) = reshape(labels == left_id, ...
        target_dim(1), target_dim(2));
    right_mask(:,:,z) = reshape(labels == right_id, ...
        target_dim(1), target_dim(2));
end

if ~any(left_mask(:)) || ~any(right_mask(:))
    error(['AAL resampling produced an empty thalamic ROI. ' ...
           'Check the target affine and atlas coordinate system.']);
end
end

function value = safe_fraction(numerator, denominator)
if denominator > 0
    value = numerator / denominator;
else
    value = NaN;
end
end
