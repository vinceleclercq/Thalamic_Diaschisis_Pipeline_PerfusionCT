function run_all()
%RUN_ALL Run the complete CT-perfusion and thalamic ROI pipeline.
%
% Before running:
%   1. Copy config/example_config.m to config/local_config.m
%   2. Edit all local paths in config/local_config.m
%
% Usage:
%   cd('PATH_TO_REPOSITORY')
%   run_all

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);

fprintf('\n============================================================\n');
fprintf('THALAMIC DIASCHISIS PIPELINE\n');
fprintf('============================================================\n');

PerfCT_Reg_MNI_001;
Extract_ROI_Thalamus_001;

fprintf('\nPipeline completed.\n');
fprintf('Review all CSV summaries and perform visual image quality control.\n');
end
