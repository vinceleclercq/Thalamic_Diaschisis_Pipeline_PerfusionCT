%% EXAMPLE LOCAL CONFIGURATION
% Copy this file to:
%   config/local_config.m
%
% Then edit the paths below.
% Do not upload local_config.m to GitHub.

cfg = struct();

%% Required local paths
cfg.root_dir      = 'D:\PATH\TO\CT_perf';
cfg.spm_dir       = 'C:\PATH\TO\spm12';
cfg.ctseg_dir     = 'C:\PATH\TO\CTseg';
cfg.fieldtrip_dir = 'C:\PATH\TO\fieldtrip';

%% Input folder conventions
cfg.subject_pattern    = 'Sub*';
cfg.anat_dicom_subdir  = 'ANAT';
cfg.tmax_dicom_subdir  = 'Perf_T';

%% Output location
cfg.group_dir_name = 'Group_results';

%% CTseg template
% Leave empty to use fullfile(cfg.ctseg_dir, 'mu_CTseg.nii').
cfg.mu_file = '';

%% AAL atlas
% Leave empty to use the atlas distributed with FieldTrip.
cfg.aal_atlas_file = '';
cfg.aal_left_label  = 'Thalamus_L';
cfg.aal_right_label = 'Thalamus_R';

%% Processing behavior
cfg.overwrite_conversion = false;
cfg.overwrite_processing = false;
cfg.create_group_maps     = true;
cfg.cleanup_intermediate_files = false;

%% Validity mask and ROI extraction
% This threshold reproduces the original script:
%   native_mask = Tmax > -100
cfg.native_mask_lower_bound = -100;

% The original ROI script excluded values <= 0.
% Keep true to reproduce the original analysis.
cfg.extract_positive_tmax_only = true;

% Store subject-level voxel values in a local MAT file.
% This output remains local and must not be uploaded to GitHub.
cfg.save_voxel_values = true;

%% CSV formatting
cfg.csv_delimiter = ';';

%% Expected study cohort size
% Used only as a warning in the validation report.
cfg.expected_subject_count = 62;

%% Numerical tolerances
cfg.affine_tolerance = 1e-4;

%% Verbosity
cfg.verbose = true;
