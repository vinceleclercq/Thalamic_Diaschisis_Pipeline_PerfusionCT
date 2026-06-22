function PerfCT_Reg_MNI_001()
%PERFCT_REG_MNI_001 Historical wrapper for conversion and registration.
%
% This wrapper preserves the original script name while using the
% refactored, configurable pipeline.

repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'src'));

cfg = load_pipeline_config(repo_root);
setup_environment(cfg);

convert_dicom_to_nifti(cfg);
register_tmax_to_template(cfg);
validate_pipeline_outputs(cfg);
end
