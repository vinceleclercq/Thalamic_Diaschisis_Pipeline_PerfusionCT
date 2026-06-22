function summary_table = convert_dicom_to_nifti(cfg)
%CONVERT_DICOM_TO_NIFTI Convert anatomical CT and Tmax DICOM series.
%
% Output in each subject directory:
%   <Subject>_Anat.nii
%   <Subject>_Perf_T.nii
%
% Summary:
%   Group_results/dicom_conversion_summary.csv

fprintf('\n============================================================\n');
fprintf('STEP 1 - DICOM TO NIFTI CONVERSION\n');
fprintf('============================================================\n');

subjects = list_subject_directories(cfg);
n = numel(subjects);

Subject = cell(n,1);
AnatStatus = cell(n,1);
TmaxStatus = cell(n,1);
AnatOutput = cell(n,1);
TmaxOutput = cell(n,1);

for i = 1:n
    subject = subjects(i).name;
    subject_dir = fullfile(cfg.root_dir, subject);

    Subject{i} = subject;
    anat_output = fullfile(subject_dir, [subject '_Anat.nii']);
    tmax_output = fullfile(subject_dir, [subject '_Perf_T.nii']);

    fprintf('\n[%d/%d] %s\n', i, n, subject);

    try
        AnatStatus{i} = convert_one_series( ...
            fullfile(subject_dir, cfg.anat_dicom_subdir), ...
            anat_output, cfg.overwrite_conversion);
    catch ME
        AnatStatus{i} = ['ERROR: ' ME.message];
        fprintf(2, '  ANAT error: %s\n', ME.message);
    end

    try
        TmaxStatus{i} = convert_one_series( ...
            fullfile(subject_dir, cfg.tmax_dicom_subdir), ...
            tmax_output, cfg.overwrite_conversion);
    catch ME
        TmaxStatus{i} = ['ERROR: ' ME.message];
        fprintf(2, '  Tmax error: %s\n', ME.message);
    end

    AnatOutput{i} = anat_output;
    TmaxOutput{i} = tmax_output;
end

summary_table = table(Subject, AnatStatus, TmaxStatus, AnatOutput, TmaxOutput);
summary_file = fullfile(cfg.group_dir, 'dicom_conversion_summary.csv');
writetable(summary_table, summary_file, 'Delimiter', cfg.csv_delimiter);

fprintf('\nDICOM conversion summary written to:\n%s\n', summary_file);
end

function status = convert_one_series(dicom_dir, output_file, overwrite)
if ~exist(dicom_dir, 'dir')
    error('DICOM directory not found: %s', dicom_dir);
end

if exist(output_file, 'file') && ~overwrite
    fprintf('  Existing output retained: %s\n', output_file);
    status = 'SKIPPED_EXISTING';
    return;
end

dicom_files = spm_select('FPList', dicom_dir, '.*');
if isempty(dicom_files)
    error('No files found in DICOM directory: %s', dicom_dir);
end

subject_dir = fileparts(output_file);
tmp_dir = fullfile(subject_dir, '_dicom_conversion_tmp');

if exist(tmp_dir, 'dir')
    rmdir(tmp_dir, 's');
end
mkdir(tmp_dir);
cleanup_tmp = onCleanup(@() cleanup_directory(tmp_dir)); %#ok<NASGU>

old_dir = pwd;
cleanup_pwd = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(tmp_dir);

headers = spm_dicom_headers(dicom_files);
spm_dicom_convert(headers, 'all', 'flat', 'nii');

created = dir(fullfile(tmp_dir, '*.nii'));
if isempty(created)
    error('SPM created no NIfTI file from: %s', dicom_dir);
end

[~, largest_idx] = max([created.bytes]);
source_file = fullfile(tmp_dir, created(largest_idx).name);

if exist(output_file, 'file')
    delete(output_file);
end
copyfile(source_file, output_file, 'f');

fprintf('  Created: %s\n', output_file);
status = 'OK';
end

function cleanup_directory(path_name)
if exist(path_name, 'dir')
    try
        rmdir(path_name, 's');
    catch
        warning('Could not remove temporary directory: %s', path_name);
    end
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
