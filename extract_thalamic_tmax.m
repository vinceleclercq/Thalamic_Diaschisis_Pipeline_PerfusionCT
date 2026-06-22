function create_group_maps(cfg)
%CREATE_GROUP_MAPS Create sum, count, mean, and coverage Tmax maps.

fprintf('\n============================================================\n');
fprintf('CREATING GROUP MAPS\n');
fprintf('============================================================\n');

maps = dir(fullfile(cfg.group_dir, '*_wfTmax.nii'));
if isempty(maps)
    warning('No normalized subject Tmax maps found in %s.', cfg.group_dir);
    return;
end

[~, order] = sort({maps.name});
maps = maps(order);

included_subjects = {};
reference_volume = [];
sum_map = [];
count_map = [];

for i = 1:numel(maps)
    tmax_file = fullfile(cfg.group_dir, maps(i).name);
    subject = erase(maps(i).name, '_wfTmax.nii');
    mask_file = fullfile(cfg.group_dir, [subject '_mask.nii']);

    if ~exist(mask_file, 'file')
        warning('Mask missing for %s; subject excluded from group maps.', ...
                subject);
        continue;
    end

    Vt = spm_vol(tmax_file);
    Vm = spm_vol(mask_file);

    if numel(Vt) ~= 1 || numel(Vm) ~= 1
        warning('Non-3D file found for %s; subject excluded.', subject);
        continue;
    end

    if isempty(reference_volume)
        reference_volume = Vt;
        sum_map = zeros(Vt.dim, 'double');
        count_map = zeros(Vt.dim, 'double');
    else
        if any(Vt.dim ~= reference_volume.dim) || ...
                max(abs(Vt.mat(:) - reference_volume.mat(:))) ...
                > cfg.affine_tolerance
            warning('Geometry mismatch for %s; subject excluded.', subject);
            continue;
        end
    end

    if any(Vm.dim ~= Vt.dim) || ...
            max(abs(Vm.mat(:) - Vt.mat(:))) > cfg.affine_tolerance
        warning('Mask geometry mismatch for %s; subject excluded.', subject);
        continue;
    end

    Yt = spm_read_vols(Vt);
    Ym = spm_read_vols(Vm);

    valid = isfinite(Yt) & Ym > 0.5;
    sum_map(valid) = sum_map(valid) + Yt(valid);
    count_map(valid) = count_map(valid) + 1;
    included_subjects{end+1,1} = subject; %#ok<AGROW>
end

if isempty(included_subjects)
    warning('No subjects were eligible for group-map creation.');
    return;
end

mean_map = zeros(size(sum_map), 'double');
valid_group = count_map > 0;
mean_map(valid_group) = sum_map(valid_group) ...
    ./ count_map(valid_group);

coverage_map = count_map / numel(included_subjects);

write_map(reference_volume, sum_map, ...
    fullfile(cfg.group_dir, 'group_sum_wfTmax.nii'));
write_map(reference_volume, count_map, ...
    fullfile(cfg.group_dir, 'group_count_valid.nii'));
write_map(reference_volume, mean_map, ...
    fullfile(cfg.group_dir, 'group_mean_wfTmax.nii'));
write_map(reference_volume, coverage_map, ...
    fullfile(cfg.group_dir, 'group_coverage_fraction.nii'));

subjects_file = fullfile(cfg.group_dir, 'group_subjects.txt');
fid = fopen(subjects_file, 'w');
if fid < 0
    error('Could not create: %s', subjects_file);
end
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', included_subjects{:});

fprintf('Group maps created from %d subjects.\n', numel(included_subjects));
end

function write_map(reference_volume, data, output_file)
Vout = reference_volume;
Vout.fname = output_file;
Vout.dt = [16 0];     % float32
Vout.pinfo = [1; 0; 0];
spm_write_vol(Vout, single(data));
end
