% Adaptation analyses - new - integration, schematization, grouping, remapping

% data_parent_dir = '/Users/mpeer/Desktop/Segmentation_subjects_data';
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';


names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'framewise_displacement', 'csf', 'white_matter'};      % Confounds to use
num_confounds = length(names);
% Reading the subjects file, with the path of all logfiles for each subject
[~, ~, subjects_file_data]= xlsread(fullfile(data_parent_dir, 'fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx'));
subjects_file_data = subjects_file_data(4:end, :);
for i = 1:size(subjects_file_data,1), if isnan(subjects_file_data{i,1}), subjects_file_data = subjects_file_data(1:i-1,:); break; end,end       % Getting rid of rows without subjects data
subj_data_dirs = subjects_file_data(:,29);
num_betas = 16;


%% Creating the design files for adaptation analyses
for s = 1:size(subjects_file_data, 1)        % Looping over subjects
    for i = 17:20       % Object viewing runs
        segmentation_create_SPM_design_objview_adaptation_integration(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_integration.mat']));
        segmentation_create_SPM_design_objview_adaptation_overlay(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_overlay.mat']));
        segmentation_create_SPM_design_objview_adaptation_group_remap(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_group_remap.mat']));
        segmentation_create_SPM_design_objview_adaptation_grmp_overlay(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_grmp_overlay.mat']));
    end
    for i = 21:23       % JRD runs
        if ~isnan(subjects_file_data{s, i})
            segmentation_create_SPM_design_jrd_adaptation_integration(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_integration.mat']));
            segmentation_create_SPM_design_jrd_adaptation_overlay(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_overlay.mat']));
            segmentation_create_SPM_design_jrd_adaptation_group_remap(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_group_remap.mat']));
            segmentation_create_SPM_design_jrd_adaptation_grmp_overlay(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation_grmp_overlay.mat']));
        end
    end
end


% Creating analysis directories for the adaptation analyses, and copying design matrix files
for subj = 1:size(subjects_file_data, 1)
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs{subj});
    
    % Defining subject analysis directories
    curr_subj_analysis_dir = fullfile(curr_subj_parent_dir, 'Analysis');    
    curr_subj_objview_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_objview_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    curr_subj_session_1_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day1/func']);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);
    % Subdirectories for different adaptation types
    curr_subj_objview_adaptation_integration_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Integration');
    curr_subj_objview_adaptation_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Overlay');
    curr_subj_objview_adaptation_group_remap_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap');
    curr_subj_objview_adaptation_grmp_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap_overlay');
    curr_subj_jrd_adaptation_integration_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Integration');
    curr_subj_jrd_adaptation_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Overlay');
    curr_subj_jrd_adaptation_group_remap_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap');
    curr_subj_jrd_adaptation_grmp_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap_overlay');
    
    % Creating the directories, if they don't already exist
    if ~exist(curr_subj_objview_adaptation_integration_dir, 'dir'), mkdir(curr_subj_objview_adaptation_integration_dir); end
    if ~exist(curr_subj_objview_adaptation_overlay_dir, 'dir'), mkdir(curr_subj_objview_adaptation_overlay_dir); end
    if ~exist(curr_subj_objview_adaptation_group_remap_dir, 'dir'), mkdir(curr_subj_objview_adaptation_group_remap_dir); end
    if ~exist(curr_subj_objview_adaptation_grmp_overlay_dir, 'dir'), mkdir(curr_subj_objview_adaptation_grmp_overlay_dir); end
    if ~exist(curr_subj_jrd_adaptation_integration_dir, 'dir'), mkdir(curr_subj_jrd_adaptation_integration_dir); end
    if ~exist(curr_subj_jrd_adaptation_overlay_dir, 'dir'), mkdir(curr_subj_jrd_adaptation_overlay_dir); end
    if ~exist(curr_subj_jrd_adaptation_group_remap_dir, 'dir'), mkdir(curr_subj_jrd_adaptation_group_remap_dir); end
    if ~exist(curr_subj_jrd_adaptation_grmp_overlay_dir, 'dir'), mkdir(curr_subj_jrd_adaptation_grmp_overlay_dir); end
    
    % Copying the design matrix files to each relevant directory
    for i=17:20     % 4 object viewing runs
        curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_integration.mat']); copyfile(curr_file, curr_subj_objview_adaptation_integration_dir);  % Object viewing adaptation files
        curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_overlay.mat']); copyfile(curr_file, curr_subj_objview_adaptation_overlay_dir);
        curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_group_remap.mat']); copyfile(curr_file, curr_subj_objview_adaptation_group_remap_dir);
        curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_grmp_overlay.mat']); copyfile(curr_file, curr_subj_objview_adaptation_grmp_overlay_dir);
    end
    for i=21:23     % 3 JRD runs
        if ~isnan(subjects_file_data{subj, i}) % JRD adaptation files, if the subject has JRD runs
            curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_integration.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_integration_dir);
            curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_overlay.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_overlay_dir);
            curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_group_remap.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_group_remap_dir);
            curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, i}, '_adaptation_grmp_overlay.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_grmp_overlay_dir);
        end
    end
    
    
    % Reading the confounds files and saving them in the right format
    % Object viewing
    object_viewing_confound_files = dir(fullfile(curr_subj_session_1_dir, '*object_viewing*confounds*.tsv'));
    object_viewing_confound_files = [object_viewing_confound_files; dir(fullfile(curr_subj_session_2_dir, '*object_viewing*confounds*.tsv'))];
    for f = 1:length(object_viewing_confound_files)
        conf = tdfread(fullfile(object_viewing_confound_files(f).folder, object_viewing_confound_files(f).name));
        fd = zeros(length(conf.csf), 1); for i = 2:length(conf.csf), fd(i) = str2double(conf.framewise_displacement(i, :)); end     % Converting framewise displacement from string to numbers
        R = [conf.trans_x, conf.trans_y, conf.trans_z, conf.rot_x, conf.rot_y, conf.rot_z, fd, conf.csf, conf.white_matter];
        save(fullfile(curr_subj_objview_adaptation_integration_dir, ['object_viewing', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_objview_adaptation_overlay_dir, ['object_viewing', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_objview_adaptation_group_remap_dir, ['object_viewing', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_objview_adaptation_grmp_overlay_dir, ['object_viewing', num2str(f), '_confounds.mat']), 'R', 'names');
    end
    % JRD
    jrd_confound_files = dir(fullfile(curr_subj_session_2_dir, '*jrd*confounds*.tsv'));
    for f = 1:length(jrd_confound_files)
        conf = tdfread(fullfile(curr_subj_session_2_dir, jrd_confound_files(f).name));
        fd = zeros(length(conf.csf), 1); for i = 2:length(conf.csf), fd(i) = str2double(conf.framewise_displacement(i, :)); end     % Converting framewise displacement from string to numbers
        R = [conf.trans_x, conf.trans_y, conf.trans_z, conf.rot_x, conf.rot_y, conf.rot_z, fd, conf.csf, conf.white_matter];
        save(fullfile(curr_subj_jrd_adaptation_integration_dir, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_jrd_adaptation_overlay_dir, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_jrd_adaptation_group_remap_dir, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
    end
end



%% Running the GLM on each run (estimating betas)
spm('defaults', 'FMRI');
for subj = 1:length(subj_data_dirs)
    disp(subj)
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs{subj});
    curr_subj_analysis_dir = fullfile(curr_subj_parent_dir, 'Analysis');
    curr_subj_session_1_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day1/func']);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);
    curr_subj_objview_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    
    curr_subj_objview_adaptation_integration_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Integration');
    curr_subj_objview_adaptation_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Overlay');
    curr_subj_objview_adaptation_group_remap_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap');
    curr_subj_objview_adaptation_grmp_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap_overlay');
    curr_subj_jrd_adaptation_integration_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Integration');
    curr_subj_jrd_adaptation_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Overlay');
    curr_subj_jrd_adaptation_group_remap_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap');
    curr_subj_jrd_adaptation_grmp_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap_overlay');

    object_viewing_adaptation_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*object_viewing*preproc_bold.nii'));
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    
    objview_adapt_dirs = {curr_subj_objview_adaptation_integration_dir, curr_subj_objview_adaptation_overlay_dir, curr_subj_objview_adaptation_group_remap_dir, curr_subj_objview_adaptation_grmp_overlay_dir};
    jrd_adapt_dirs = {curr_subj_jrd_adaptation_integration_dir, curr_subj_jrd_adaptation_overlay_dir, curr_subj_jrd_adaptation_group_remap_dir, curr_subj_jrd_adaptation_grmp_overlay_dir};
    
    for d=1:length(objview_adapt_dirs)
        % Defining design and running GLM for object viewing runs - adaptation analysis - day 2 only
        clear matlabbatch;
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess = struct();
        design_file = dir(fullfile(objview_adapt_dirs{d}, '*Segmentation*.mat'));
        confounds_file = dir(fullfile(objview_adapt_dirs{d}, '*confounds.mat'));
        for r = 1:2     % Defining specific parameters for each run
            n = nifti(fullfile(object_viewing_adaptation_data_files(r).folder, object_viewing_adaptation_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.dir = {objview_adapt_dirs{d}};    % The output directory
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(objview_adapt_dirs{d}, design_file(r+2).name)};      % The design file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(objview_adapt_dirs{d}, confounds_file(r+2).name)};      % The confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(object_viewing_adaptation_data_files(r).folder, object_viewing_adaptation_data_files(r).name), ',', num2str(v)];
            end
        end
        if ~exist(fullfile(objview_adapt_dirs{d}, 'SPM.mat'), 'file')        % Checking if SPM.mat file already exists, running if not
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(objview_adapt_dirs{d}, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
            % Running the design creation and estimation script
            spm_jobman('run', matlabbatch);
        end

        % Defining design and running GLM for JRD runs - adaptation analysis
        clear matlabbatch;
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        design_file = dir(fullfile(jrd_adapt_dirs{d}, '*Segmentation*.mat'));
        confounds_file = dir(fullfile(jrd_adapt_dirs{d}, '*confounds.mat'));
        exist_jrd = 0;
        for r = 1:length(jrd_data_files)     % Defining specific parameters for each run
            if exist(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name))
                exist_jrd = 1;
                n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
                matlabbatch{1}.spm.stats.fmri_spec.dir = {jrd_adapt_dirs{d}};    % The output directory
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(jrd_adapt_dirs{d}, design_file(r).name)};      % The design file
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(jrd_adapt_dirs{d}, confounds_file(r).name)};      % The confounds file
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
                for v = 1:num_volumes
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
                end
            end
        end
        if ~exist(fullfile(jrd_adapt_dirs{d}, 'SPM.mat'), 'file') & exist_jrd == 1        % Checking if SPM.mat file already exists, running if not
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(jrd_adapt_dirs{d}, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
            % Running the design creation and estimation script
            spm_jobman('run', matlabbatch);
        end

    end
end




%% Creating contrasts and computing T-values
for subj = 1:length(subj_data_dirs)
    disp(subj)
    
    % Defining directories
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'Analysis');
    curr_subj_objview_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_objview_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    
    curr_subj_objview_adaptation_integration_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Integration');
    curr_subj_objview_adaptation_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Overlay');
    curr_subj_objview_adaptation_group_remap_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap');
    curr_subj_objview_adaptation_grmp_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap_overlay');
    curr_subj_jrd_adaptation_integration_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Integration');
    curr_subj_jrd_adaptation_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Overlay');
    curr_subj_jrd_adaptation_group_remap_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap');
    curr_subj_jrd_adaptation_grmp_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap_overlay');
    
    
    % Integration model contrasts
    % Object viewing 
    SPMmat_filename = fullfile(curr_subj_objview_adaptation_integration_dir, 'SPM.mat');
    % Deleting existing contrasts
    load(SPMmat_filename);
    if ~isempty(SPM.xCon)
        SPM.xCon = []; save(SPMmat_filename, 'SPM');
        temp_con_files = dir(fullfile(curr_subj_objview_adaptation_integration_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objview_adaptation_integration_dir, temp_con_files(i).name)); end
        temp_t_files = dir(fullfile(curr_subj_objview_adaptation_integration_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objview_adaptation_integration_dir, temp_t_files(i).name)); end
    end
    % Creating the new contrasts and T images
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Distance effect parametric modulation';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 0, 0, 1, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'All stimuli regressors vs rest';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([1, 1, 1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    spm_jobman('run', matlabbatch);
    % JRD
    SPMmat_filename = fullfile(curr_subj_jrd_adaptation_integration_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_adaptation_integration_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_adaptation_integration_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_adaptation_integration_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_adaptation_integration_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_adaptation_integration_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        clear matlabbatch;
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Distance effect parametric modulation';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 0, 1, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'All stimuli regressors vs rest';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([1, 1, 0, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end
    
    
    
    % Overlay model contrasts
    % Object viewing 
    SPMmat_filename = fullfile(curr_subj_objview_adaptation_overlay_dir, 'SPM.mat');
    % Deleting existing contrasts
    load(SPMmat_filename);
    if ~isempty(SPM.xCon)
        SPM.xCon = []; save(SPMmat_filename, 'SPM');
        temp_con_files = dir(fullfile(curr_subj_objview_adaptation_overlay_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objview_adaptation_overlay_dir, temp_con_files(i).name)); end
        temp_t_files = dir(fullfile(curr_subj_objview_adaptation_overlay_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objview_adaptation_overlay_dir, temp_t_files(i).name)); end
    end
    % Creating the new contrasts and T images
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Distance effect parametric modulation';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 0, 0, 1, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'All stimuli regressors vs rest';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([1, 1, 1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    % JRD
    SPMmat_filename = fullfile(curr_subj_jrd_adaptation_overlay_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_adaptation_overlay_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_adaptation_overlay_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_adaptation_overlay_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_adaptation_overlay_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_adaptation_overlay_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        clear matlabbatch;
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Distance effect parametric modulation';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 0, 1, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'All stimuli regressors vs rest';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([1, 1, 0, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        spm_jobman('run', matlabbatch);
    end

    
    
    % Grouping and remapping model contrasts
    % Object viewing 
    SPMmat_filename = fullfile(curr_subj_objview_adaptation_group_remap_dir, 'SPM.mat');
    % Deleting existing contrasts
    load(SPMmat_filename);
    if ~isempty(SPM.xCon)
        SPM.xCon = []; save(SPMmat_filename, 'SPM');
        temp_con_files = dir(fullfile(curr_subj_objview_adaptation_group_remap_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objview_adaptation_group_remap_dir, temp_con_files(i).name)); end
        temp_t_files = dir(fullfile(curr_subj_objview_adaptation_group_remap_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objview_adaptation_group_remap_dir, temp_t_files(i).name)); end
    end
    % Creating the new contrasts and T images
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Grouping effect';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 1, 0, -1, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Distance coding between segments';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([0, 1, 0, 0, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Distance coding within segment';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [repmat([0, 0, 0, 1, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Euclidean distance effect';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [repmat([0, 1, 0, 1, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'All stimuli regressors vs rest';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [repmat([1, 0, 1, 0, 1, 1, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    spm_jobman('run', matlabbatch);
    % JRD
    SPMmat_filename = fullfile(curr_subj_jrd_adaptation_group_remap_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_adaptation_group_remap_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_adaptation_group_remap_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_adaptation_group_remap_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_adaptation_group_remap_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_adaptation_group_remap_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        clear matlabbatch;
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Grouping effect';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 1, 0, -1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Distance coding between segments';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([0, 1, 0, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Distance coding within segment';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [repmat([0, 0, 0, 1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Euclidean distance effect';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [repmat([0, 1, 0, 1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'All stimuli regressors vs rest';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [repmat([1, 0, 1, 0, 1, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        spm_jobman('run', matlabbatch);
    end
    
    
    
    % Grouping and remapping model contrasts - overlay
    % Object viewing 
    SPMmat_filename = fullfile(curr_subj_objview_adaptation_grmp_overlay_dir, 'SPM.mat');
    % Deleting existing contrasts
    load(SPMmat_filename);
    if ~isempty(SPM.xCon)
        SPM.xCon = []; save(SPMmat_filename, 'SPM');
        temp_con_files = dir(fullfile(curr_subj_objview_adaptation_grmp_overlay_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objview_adaptation_grmp_overlay_dir, temp_con_files(i).name)); end
        temp_t_files = dir(fullfile(curr_subj_objview_adaptation_grmp_overlay_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objview_adaptation_grmp_overlay_dir, temp_t_files(i).name)); end
    end
    % Creating the new contrasts and T images
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Grouping effect';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 1, 0, -1, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Distance coding between segments';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([0, 1, 0, 0, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Distance coding within segment';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [repmat([0, 0, 0, 1, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Euclidean distance effect';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [repmat([0, 1, 0, 1, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'All stimuli regressors vs rest';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [repmat([1, 0, 1, 0, 1, 1, zeros(1, num_confounds)], 1, 2), 0, 0];
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    spm_jobman('run', matlabbatch);
    % JRD
    SPMmat_filename = fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_adaptation_grmp_overlay_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        clear matlabbatch;
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Grouping effect';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 1, 0, -1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Distance coding between segments';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([0, 1, 0, 0, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Distance coding within segment';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [repmat([0, 0, 0, 1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Euclidean distance effect';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [repmat([0, 1, 0, 1, 0, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'All stimuli regressors vs rest';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [repmat([1, 0, 1, 0, 1, zeros(1, num_confounds)], 1, 2), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        spm_jobman('run', matlabbatch);
    end
    
end





%% Analyzing ROI results

% Segmentation adaptation analysis

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC', 'lHCpost','rHCpost','lHCant','rHCant'};
num_ROIs = length(all_ROI_names);

% Identify subjects and create empty variables
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
means_tvals_all = cell(6,1);
nums_contrasts = [2,2,5,2,2,5];
for d=1:6
    means_tvals_all{d} = nan(length(subj_data_dirs), num_ROIs, nums_contrasts(d));
end
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    
    curr_subj_objview_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    curr_subj_objview_adaptation_integration_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Integration');
    curr_subj_objview_adaptation_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Overlay');
    curr_subj_objview_adaptation_group_remap_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap');
    curr_subj_objview_adaptation_grmp_overlay_dir = fullfile(curr_subj_objview_adaptation_analysis_dir, 'Group_remap_overlay');
    
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    curr_subj_jrd_adaptation_integration_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Integration');
    curr_subj_jrd_adaptation_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Overlay');
    curr_subj_jrd_adaptation_group_remap_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap');
    curr_subj_jrd_adaptation_grmp_overlay_dir = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'Group_remap_overlay');
    
    % Loading the subject specific ROIs
    all_ROIs_individual = {};
    for i = 1:length(all_ROI_names)
        all_ROIs_individual{i} = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
    end    
    num_ROIs_individual = length(all_ROIs_individual);
    
    current_subj_adaptation_dirs = {curr_subj_objview_adaptation_integration_dir, curr_subj_objview_adaptation_overlay_dir, curr_subj_objview_adaptation_group_remap_dir, curr_subj_objview_adaptation_grmp_overlay_dir, ...
        curr_subj_jrd_adaptation_integration_dir, curr_subj_jrd_adaptation_overlay_dir, curr_subj_jrd_adaptation_group_remap_dir, curr_subj_jrd_adaptation_grmp_overlay_dir};
    for d = 1:length(current_subj_adaptation_dirs)
        % Finding the contrast t-value files
        curr_tval_files = dir(fullfile(current_subj_adaptation_dirs{d}, 'spmT*.nii'));
        if ~isempty(curr_tval_files)
            % Adaptation analysis results averaging in each ROI
            for roi = 1:num_ROIs_individual        % Going over ROIs
                curr_ROI = find(all_ROIs_individual{roi}(:)==1);
                for t = 1:length(curr_tval_files)
                    curr_tvals = spm_read_vols(spm_vol(fullfile(curr_tval_files(t).folder, curr_tval_files(t).name)));
                    means_tvals_all{d}(subj,roi,t) = nanmean(curr_tvals(curr_ROI));
                end
            end
        end
    end
end

for d = 1:length(current_subj_adaptation_dirs)
    for t = 1:size(means_tvals_all{d}, 3)
        for roi = 1:num_ROIs_individual        % Going over ROIs
            pvals{d}(t,roi) = signrank(means_tvals_all{d}(:,roi,t),[],'tail','right');
        end
    end
end
% Cells - 1.objview_integration; 2.objview_overlay; 3.objview_group_remap; 4.jrd_integration; 5.jrd_overlay; 6.jrd_group_remap
% Contrasts - integration - 1.distance effect assuming integration; 2.all stimuli vs. rest
% Contrasts - overlay - 1.distance effect assuming overlay; 2.all stimuli vs. rest
% Contrasts - group_remap - 1.remapping effect (between-segment beta larger than within); 2.between-seg distance coding; 3.within-seg distance coding; 4. overall distance effect; 5. all vs. rest
