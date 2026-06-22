function Extract_ROI_Thalamus_001()
%EXTRACT_ROI_THALAMUS_001 Extract bilateral thalamic Tmax measurements.
%
% This function:
%   1. Loads the AAL atlas distributed with FieldTrip, unless another
%      atlas path is specified in config/local_config.m.
%   2. Identifies the Thalamus_L and Thalamus_R labels.
%   3. Resamples these labels to each normalized Tmax grid using
%      nearest-neighbour interpolation.
%   4. Applies the subject-specific validity mask.
%   5. Extracts mean, median, standard deviation, voxel count, and
%      coverage for both thalami.
%   6. Writes thalamus_CT_values.csv and, optionally, a MAT file.
%
% Requirements:
%   MATLAB, SPM12, and FieldTrip.
%
% The function reads local paths from config/local_config.m.

repo_root = fileparts(mfilename('fullpath'));
cfg = load_local_config(repo_root);
initialize_dependencies(cfg);

group_dir = fullfile(cfg.root_dir, 'Group_results');
if ~exist(group_dir, 'dir')
    error('Group_results directory not found: %s', group_dir);
end

atlas_file = resolve_atlas_file(cfg);
atlas = ft_read_atlas(atlas_file);

left_label_id = find_atlas_label(atlas, cfg.aal_left_label);
right_label_id = find_atlas_label(atlas, cfg.aal_right_label);

fprintf('\nAAL atlas: %s\n', atlas_file);
fprintf('Left thalamus : %s (label %d)\n', ...
    cfg.aal_left_label, left_label_id);
fprintf('Right thalamus: %s (label %d)\n', ...
    cfg.aal_right_label, right_label_id);

tmax_files = dir(fullfile(group_dir, '*_wfTmax.nii'));
if isempty(tmax_files)
    error('No *_wfTmax.nii files found in %s', group_dir);
end

[~, order] = sort({tmax_files.name});
tmax_files = tmax_files(order);
n_subjects = numel(tmax_files);

Subject = cell(n_subjects,1);
Status = cell(n_subjects,1);
TmaxFile = cell(n_subjects,1);
MaskFile = cell(n_subjects,1);
VoxelSizeX_mm = NaN(n_subjects,1);
VoxelSizeY_mm = NaN(n_subjects,1);
VoxelSizeZ_mm = NaN(n_subjects,1);

Thal_L_Mean = NaN(n_subjects,1);
Thal_L_Median = NaN(n_subjects,1);
Thal_L_Std = NaN(n_subjects,1);
Thal_L_nAtlasVoxels = NaN(n_subjects,1);
Thal_L_nValidVoxels = NaN(n_subjects,1);
Thal_L_Coverage = NaN(n_subjects,1);

Thal_R_Mean = NaN(n_subjects,1);
Thal_R_Median = NaN(n_subjects,1);
Thal_R_Std = NaN(n_subjects,1);
Thal_R_nAtlasVoxels = NaN(n_subjects,1);
Thal_R_nValidVoxels = NaN(n_subjects,1);
Thal_R_Coverage = NaN(n_subjects,1);

left_voxel_values = cell(n_subjects,1);
right_voxel_values = cell(n_subjects,1);

cached_dimensions = [];
cached_affine = [];
cached_left_mask = [];
cached_right_mask = [];

for i = 1:n_subjects
    subject = erase(tmax_files(i).name, '_wfTmax.nii');
    tmax_file = fullfile(group_dir, tmax_files(i).name);
    mask_file = fullfile(group_dir, [subject '_mask.nii']);

    Subject{i} = subject;
    Status{i} = 'OK';
    TmaxFile{i} = tmax_file;
    MaskFile{i} = mask_file;

    fprintf('\n------------------------------------------------------------\n');
    fprintf('ROI extraction: %s (%d/%d)\n', ...
        subject, i, n_subjects);
    fprintf('------------------------------------------------------------\n');

    try
        if ~exist(mask_file, 'file')
            error('Validity mask not found: %s', mask_file);
        end

        Vt = spm_vol(tmax_file);
        Vm = spm_vol(mask_file);
        assert_single_volume(Vt, tmax_file);
        assert_single_volume(Vm, mask_file);
        check_same_geometry(Vt, Vm);

        voxel_size = sqrt(sum(Vt.mat(1:3,1:3).^2, 1));
        VoxelSizeX_mm(i) = voxel_size(1);
        VoxelSizeY_mm(i) = voxel_size(2);
        VoxelSizeZ_mm(i) = voxel_size(3);

        if isempty(cached_dimensions) || ...
                any(Vt.dim ~= cached_dimensions) || ...
                max(abs(Vt.mat(:) - cached_affine(:))) > 1e-4

            fprintf('Resampling AAL thalamic labels to target grid...\n');

            [cached_left_mask, cached_right_mask] = ...
                resample_atlas_labels_nearest( ...
                    atlas, left_label_id, ...
                    right_label_id, Vt);

            cached_dimensions = Vt.dim;
            cached_affine = Vt.mat;
        end

        Yt = double(spm_read_vols(Vt));
        Ym = double(spm_read_vols(Vm));

        valid_data = isfinite(Yt) & Ym > 0.5;
        if cfg.extract_positive_values_only
            valid_data = valid_data & Yt > 0;
        end

        left_atlas_voxels = nnz(cached_left_mask);
        right_atlas_voxels = nnz(cached_right_mask);

        left_valid_mask = cached_left_mask & valid_data;
        right_valid_mask = cached_right_mask & valid_data;

        left_values = Yt(left_valid_mask);
        right_values = Yt(right_valid_mask);

        if isempty(left_values)
            warning('No valid left-thalamic voxels for %s.', subject);
        else
            Thal_L_Mean(i) = mean(left_values);
            Thal_L_Median(i) = median(left_values);
            Thal_L_Std(i) = std(left_values);
            left_voxel_values{i} = left_values;
        end

        if isempty(right_values)
            warning('No valid right-thalamic voxels for %s.', subject);
        else
            Thal_R_Mean(i) = mean(right_values);
            Thal_R_Median(i) = median(right_values);
            Thal_R_Std(i) = std(right_values);
            right_voxel_values{i} = right_values;
        end

        Thal_L_nAtlasVoxels(i) = left_atlas_voxels;
        Thal_L_nValidVoxels(i) = numel(left_values);
        Thal_L_Coverage(i) = safe_divide( ...
            numel(left_values), left_atlas_voxels);

        Thal_R_nAtlasVoxels(i) = right_atlas_voxels;
        Thal_R_nValidVoxels(i) = numel(right_values);
        Thal_R_Coverage(i) = safe_divide( ...
            numel(right_values), right_atlas_voxels);

        fprintf('Left thalamus : mean %.4f; valid %d/%d voxels\n', ...
            Thal_L_Mean(i), numel(left_values), left_atlas_voxels);
        fprintf('Right thalamus: mean %.4f; valid %d/%d voxels\n', ...
            Thal_R_Mean(i), numel(right_values), right_atlas_voxels);

    catch ME
        Status{i} = ['ERROR: ' ME.message];
        fprintf(2, 'ERROR for %s: %s\n', subject, ME.message);
    end
end

results_table = table(Subject, Status, TmaxFile, MaskFile, ...
    VoxelSizeX_mm, VoxelSizeY_mm, VoxelSizeZ_mm, ...
    Thal_L_Mean, Thal_L_Median, Thal_L_Std, ...
    Thal_L_nAtlasVoxels, Thal_L_nValidVoxels, Thal_L_Coverage, ...
    Thal_R_Mean, Thal_R_Median, Thal_R_Std, ...
    Thal_R_nAtlasVoxels, Thal_R_nValidVoxels, Thal_R_Coverage);

output_csv = fullfile(group_dir, 'thalamus_CT_values.csv');
writetable(results_table, output_csv, ...
    'Delimiter', cfg.csv_delimiter);

fprintf('\nThalamic summary:\n%s\n', output_csv);

if cfg.save_voxel_values
    output_mat = fullfile(group_dir, 'thalamus_CT_results.mat');

    metadata = struct();
    metadata.atlas_file = atlas_file;
    metadata.left_label = cfg.aal_left_label;
    metadata.right_label = cfg.aal_right_label;
    metadata.extract_positive_values_only = ...
        cfg.extract_positive_values_only;
    metadata.created = datestr(now, 30); %#ok<DATST>

    save(output_mat, 'results_table', ...
        'left_voxel_values', 'right_voxel_values', ...
        'metadata', '-v7.3');

    fprintf('Voxel-level local results:\n%s\n', output_mat);
end

if isfield(cfg, 'expected_subject_count') && ...
        ~isempty(cfg.expected_subject_count) && ...
        n_subjects ~= cfg.expected_subject_count
    warning('Found %d subjects; expected %d.', ...
        n_subjects, cfg.expected_subject_count);
end

fprintf('\nExtract_ROI_Thalamus_001 completed.\n');
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

required_fields = {'root_dir','spm_dir', ...
                   'fieldtrip_dir'};
for i = 1:numel(required_fields)
    field_name = required_fields{i};
    if ~isfield(cfg, field_name) || isempty(cfg.(field_name))
        error('Missing configuration field: cfg.%s', field_name);
    end
end

if ~isfield(cfg, 'aal_atlas_file')
    cfg.aal_atlas_file = '';
end
if ~isfield(cfg, 'aal_left_label')
    cfg.aal_left_label = 'Thalamus_L';
end
if ~isfield(cfg, 'aal_right_label')
    cfg.aal_right_label = 'Thalamus_R';
end
if ~isfield(cfg, 'extract_positive_values_only')
    cfg.extract_positive_values_only = true;
end
if ~isfield(cfg, 'save_voxel_values')
    cfg.save_voxel_values = true;
end
if ~isfield(cfg, 'csv_delimiter')
    cfg.csv_delimiter = ';';
end

if ~exist(cfg.root_dir, 'dir')
    error('Root data directory not found: %s', cfg.root_dir);
end
if ~exist(cfg.spm_dir, 'dir')
    error('SPM directory not found: %s', cfg.spm_dir);
end
if ~exist(cfg.fieldtrip_dir, 'dir')
    error('FieldTrip directory not found: %s', ...
          cfg.fieldtrip_dir);
end
end

function initialize_dependencies(cfg)
addpath(cfg.spm_dir);
if exist('spm', 'file') ~= 2
    error('SPM was not found after adding %s', cfg.spm_dir);
end

spm('defaults', 'fmri');
spm_jobman('initcfg');

addpath(cfg.fieldtrip_dir);
if exist('ft_defaults', 'file') ~= 2
    error('FieldTrip was not found after adding %s', ...
          cfg.fieldtrip_dir);
end
ft_defaults;
end

function atlas_file = resolve_atlas_file(cfg)
if isempty(cfg.aal_atlas_file)
    atlas_file = fullfile(fileparts(which('ft_defaults')), ...
        'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');
else
    atlas_file = cfg.aal_atlas_file;
end

if ~exist(atlas_file, 'file')
    error('AAL atlas not found: %s', atlas_file);
end
end

function label_id = find_atlas_label(atlas, label_name)
matches = find(strcmpi(atlas.tissuelabel, label_name));

if numel(matches) ~= 1
    thalamic_labels = atlas.tissuelabel( ...
        contains(atlas.tissuelabel, ...
        'Thalamus', 'IgnoreCase', true));

    error(['Expected exactly one atlas label "%s", found %d. ' ...
           'Available thalamic labels: %s'], ...
          label_name, numel(matches), ...
          strjoin(thalamic_labels, ', '));
end

label_id = matches(1);
end

function assert_single_volume(V, file_name)
if numel(V) ~= 1
    error('Expected one 3D volume in %s, found %d.', ...
          file_name, numel(V));
end
end

function check_same_geometry(Va, Vb)
if any(Va.dim ~= Vb.dim)
    error('Tmax and mask dimensions do not match.');
end
if max(abs(Va.mat(:) - Vb.mat(:))) > 1e-4
    error('Tmax and mask affine matrices do not match.');
end
end

function [left_mask, right_mask] = ...
    resample_atlas_labels_nearest( ...
        atlas, left_label_id, right_label_id, Vtarget)

target_dimensions = Vtarget.dim;
left_mask = false(target_dimensions);
right_mask = false(target_dimensions);

atlas_inverse = inv(atlas.transform);

[x_grid, y_grid] = ndgrid( ...
    1:target_dimensions(1), ...
    1:target_dimensions(2));

points_per_slice = numel(x_grid);

for z_index = 1:target_dimensions(3)
    target_voxels = [ ...
        x_grid(:)'; ...
        y_grid(:)'; ...
        repmat(z_index, 1, points_per_slice); ...
        ones(1, points_per_slice)];

    world_coordinates = Vtarget.mat * target_voxels;
    atlas_coordinates = atlas_inverse * world_coordinates;

    atlas_x = round(atlas_coordinates(1,:));
    atlas_y = round(atlas_coordinates(2,:));
    atlas_z = round(atlas_coordinates(3,:));

    in_bounds = atlas_x >= 1 & atlas_x <= atlas.dim(1) & ...
                atlas_y >= 1 & atlas_y <= atlas.dim(2) & ...
                atlas_z >= 1 & atlas_z <= atlas.dim(3);

    labels = zeros(1, points_per_slice);

    valid_indices = sub2ind(atlas.dim, ...
        atlas_x(in_bounds), ...
        atlas_y(in_bounds), ...
        atlas_z(in_bounds));

    labels(in_bounds) = ...
        double(atlas.tissue(valid_indices));

    left_mask(:,:,z_index) = reshape( ...
        labels == left_label_id, ...
        target_dimensions(1), target_dimensions(2));

    right_mask(:,:,z_index) = reshape( ...
        labels == right_label_id, ...
        target_dimensions(1), target_dimensions(2));
end

if ~any(left_mask(:))
    error('The resampled left-thalamic ROI is empty.');
end
if ~any(right_mask(:))
    error('The resampled right-thalamic ROI is empty.');
end
end

function result = safe_divide(numerator, denominator)
if denominator > 0
    result = numerator / denominator;
else
    result = NaN;
end
end
