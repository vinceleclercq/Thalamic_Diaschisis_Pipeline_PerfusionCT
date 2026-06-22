function run_all()
%RUN_ALL Run the complete CT-perfusion thalamic extraction pipeline.
%
% Usage:
%   cd('PATH_TO_REPOSITORY')
%   run_all

repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'src'));

cfg = load_pipeline_config(repo_root);
setup_environment(cfg);

fprintf('\n============================================================\n');
fprintf('THALAMIC DIASCHISIS CT-PERFUSION PIPELINE\n');
fprintf('Repository: %s\n', repo_root);
fprintf('Data root : %s\n', cfg.root_dir);
fprintf('============================================================\n');

convert_dicom_to_nifti(cfg);
register_tmax_to_template(cfg);
extract_thalamic_tmax(cfg);
validate_pipeline_outputs(cfg);

fprintf('\nPipeline completed.\n');
fprintf('Review all processing summaries and perform visual QC before analysis.\n');
end
