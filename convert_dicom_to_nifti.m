function setup_environment(cfg)
%SETUP_ENVIRONMENT Add dependencies and initialize SPM and FieldTrip.

addpath(cfg.spm_dir);
addpath(cfg.ctseg_dir);

if exist('spm', 'file') ~= 2
    error('SPM was not found after adding: %s', cfg.spm_dir);
end
if exist('spm_CTseg', 'file') ~= 2
    error('spm_CTseg.m was not found after adding: %s', cfg.ctseg_dir);
end

spm('defaults', 'fmri');
spm_jobman('initcfg');

% FieldTrip recommends adding only its root directory before ft_defaults.
addpath(cfg.fieldtrip_dir);
if exist('ft_defaults', 'file') ~= 2
    error('FieldTrip ft_defaults.m was not found after adding: %s', ...
          cfg.fieldtrip_dir);
end
ft_defaults;

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
