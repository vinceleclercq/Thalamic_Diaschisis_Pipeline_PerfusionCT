function cfg = load_pipeline_config(repo_root)
%LOAD_PIPELINE_CONFIG Load and validate config/local_config.m.

if nargin < 1 || isempty(repo_root)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
end

config_file = fullfile(repo_root, 'config', 'local_config.m');
example_file = fullfile(repo_root, 'config', 'example_config.m');

if ~exist(config_file, 'file')
    error(['Missing local configuration file:\n  %s\n\n' ...
           'Copy:\n  %s\n\nto:\n  %s\n\nand edit the local paths.'], ...
           config_file, example_file, config_file);
end

cfg = struct(); %#ok<NASGU>
run(config_file);

if ~exist('cfg', 'var') || ~isstruct(cfg)
    error('local_config.m must create a structure named cfg.');
end

cfg.repo_root = repo_root;

required = {'root_dir', 'spm_dir', 'ctseg_dir', 'fieldtrip_dir'};
for k = 1:numel(required)
    field_name = required{k};
    if ~isfield(cfg, field_name) || isempty(cfg.(field_name))
        error('Missing required configuration field: cfg.%s', field_name);
    end
end

cfg = set_default(cfg, 'subject_pattern', 'Sub*');
cfg = set_default(cfg, 'anat_dicom_subdir', 'ANAT');
cfg = set_default(cfg, 'tmax_dicom_subdir', 'Perf_T');
cfg = set_default(cfg, 'group_dir_name', 'Group_results');
cfg = set_default(cfg, 'mu_file', '');
cfg = set_default(cfg, 'aal_atlas_file', '');
cfg = set_default(cfg, 'aal_left_label', 'Thalamus_L');
cfg = set_default(cfg, 'aal_right_label', 'Thalamus_R');
cfg = set_default(cfg, 'overwrite_conversion', false);
cfg = set_default(cfg, 'overwrite_processing', false);
cfg = set_default(cfg, 'create_group_maps', true);
cfg = set_default(cfg, 'cleanup_intermediate_files', false);
cfg = set_default(cfg, 'native_mask_lower_bound', -100);
cfg = set_default(cfg, 'extract_positive_tmax_only', true);
cfg = set_default(cfg, 'save_voxel_values', true);
cfg = set_default(cfg, 'csv_delimiter', ';');
cfg = set_default(cfg, 'expected_subject_count', []);
cfg = set_default(cfg, 'affine_tolerance', 1e-4);
cfg = set_default(cfg, 'verbose', true);

cfg.root_dir = char(cfg.root_dir);
cfg.spm_dir = char(cfg.spm_dir);
cfg.ctseg_dir = char(cfg.ctseg_dir);
cfg.fieldtrip_dir = char(cfg.fieldtrip_dir);
cfg.group_dir = fullfile(cfg.root_dir, cfg.group_dir_name);

if isempty(cfg.mu_file)
    cfg.mu_file = fullfile(cfg.ctseg_dir, 'mu_CTseg.nii');
end

if ~exist(cfg.root_dir, 'dir')
    error('Data root directory not found: %s', cfg.root_dir);
end
if ~exist(cfg.spm_dir, 'dir')
    error('SPM directory not found: %s', cfg.spm_dir);
end
if ~exist(cfg.ctseg_dir, 'dir')
    error('CTseg directory not found: %s', cfg.ctseg_dir);
end
if ~exist(cfg.fieldtrip_dir, 'dir')
    error('FieldTrip directory not found: %s', cfg.fieldtrip_dir);
end
if ~exist(cfg.mu_file, 'file')
    error('CTseg template not found: %s', cfg.mu_file);
end

if ~exist(cfg.group_dir, 'dir')
    mkdir(cfg.group_dir);
end
end

function cfg = set_default(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end
