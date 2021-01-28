%% Initial data definitions
% data_parent_dir = '/Users/mpeer/Desktop/Segmentation_subjects_data';
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data';
names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'framewise_displacement', 'csf', 'white_matter'};      % Confounds to use
num_confounds = length(names);
% Reading the subjects file, with the path of all logfiles for each subject
[~, ~, subjects_file_data]= xlsread(fullfile(data_parent_dir, 'fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx')); 
subjects_file_data = subjects_file_data(4:end, :);
for i = 1:size(subjects_file_data,1), if isnan(subjects_file_data{i,1}), subjects_file_data = subjects_file_data(1:i-1,:); break; end,end       % Getting rid of rows without subjects data
subj_data_dirs = subjects_file_data(:,29);
num_betas = 16;


%% Creating design matrix files for SPM analysis
for s = 1:size(subjects_file_data, 1)        % Looping over subjects
    for i = 15:16       % Functional localizer runs
        % segmentation_create_SPM_design_functional_localizer(fullfile(a{s, 8}, a{s, i}, '.csv'), fullfile(a{s, 8}, a{s, i}, '.mat'));    % PC directory
        segmentation_create_SPM_design_functional_localizer(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.mat']));   % Mac directory
    end
    for i = 17:20       % Object viewing runs
        segmentation_create_SPM_design_object_viewing(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.mat']));
        segmentation_create_SPM_design_object_viewing_adaptation(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation.mat']));
    end
    for i = 21:23       % JRD runs
        if ~isnan(subjects_file_data{s, i})
            segmentation_create_SPM_design_JRD(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.mat']));
            segmentation_create_SPM_design_jrd_adaptation(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_adaptation.mat']));
            segmentation_create_SPM_design_JRD_univariate(fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '.csv']), fullfile(subjects_file_data{s, 27}, [subjects_file_data{s, i}, '_univariate.mat']));
        end
    end
end

% Creating analysis directories and copying design matrix files
for subj = 1:size(subjects_file_data, 1)
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs{subj});
    
    % Creating subject analysis directories
    curr_subj_analysis_dir = fullfile(curr_subj_parent_dir, 'Analysis');
    if ~exist(curr_subj_analysis_dir, 'dir'), mkdir(curr_subj_analysis_dir); end
    
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_objectviewing_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    curr_subj_jrd_univariate_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_univariate');
    if ~exist(curr_subj_objectviewing_analysis_dir, 'dir'), mkdir(curr_subj_objectviewing_analysis_dir); end
    if ~exist(curr_subj_localizer_analysis_dir, 'dir'), mkdir(curr_subj_localizer_analysis_dir); end
    if ~exist(curr_subj_jrd_analysis_dir, 'dir'), mkdir(curr_subj_jrd_analysis_dir); end
    if ~exist(curr_subj_objectviewing_adaptation_analysis_dir, 'dir'), mkdir(curr_subj_objectviewing_adaptation_analysis_dir); end
    if ~exist(curr_subj_jrd_adaptation_analysis_dir, 'dir'), mkdir(curr_subj_jrd_adaptation_analysis_dir); end
    if ~exist(curr_subj_jrd_univariate_analysis_dir, 'dir'), mkdir(curr_subj_jrd_univariate_analysis_dir); end
    
    curr_subj_objectviewing_analysis_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run4')};
    for i = 1:length(curr_subj_objectviewing_analysis_dirs)
        if ~exist(curr_subj_objectviewing_analysis_dirs{i}, 'dir')
            mkdir(curr_subj_objectviewing_analysis_dirs{i});
        end
    end
    
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    for i = 1:length(curr_subj_jrd_analysis_dirs)
        if ~exist(curr_subj_jrd_analysis_dirs{i}, 'dir')
            mkdir(curr_subj_jrd_analysis_dirs{i});
        end
    end
    
    curr_subj_objectviewing_combined_day1_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day1');
    curr_subj_objectviewing_combined_day2_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day2');
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');
    if ~exist(curr_subj_objectviewing_combined_day1_dir, 'dir'), mkdir(curr_subj_objectviewing_combined_day1_dir); end
    if ~exist(curr_subj_objectviewing_combined_day2_dir, 'dir'), mkdir(curr_subj_objectviewing_combined_day2_dir); end
    if ~exist(curr_subj_jrd_combined_dir, 'dir'), mkdir(curr_subj_jrd_combined_dir); end

    
    % Copying the design matrix files to each relevant directory
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 15}, '.mat']); copyfile(curr_file, curr_subj_localizer_analysis_dir);  % Functional localizer files
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 16}, '.mat']); copyfile(curr_file, curr_subj_localizer_analysis_dir);
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 17}, '.mat']); copyfile(curr_file, curr_subj_objectviewing_analysis_dirs{1});  % Object viewing files
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 18}, '.mat']); copyfile(curr_file, curr_subj_objectviewing_analysis_dirs{2});
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 19}, '.mat']); copyfile(curr_file, curr_subj_objectviewing_analysis_dirs{3});
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 20}, '.mat']); copyfile(curr_file, curr_subj_objectviewing_analysis_dirs{4});
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 17}, '_adaptation.mat']); copyfile(curr_file, curr_subj_objectviewing_adaptation_analysis_dir);  % Object viewing adaptation files
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 18}, '_adaptation.mat']); copyfile(curr_file, curr_subj_objectviewing_adaptation_analysis_dir);
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 19}, '_adaptation.mat']); copyfile(curr_file, curr_subj_objectviewing_adaptation_analysis_dir);
    curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 20}, '_adaptation.mat']); copyfile(curr_file, curr_subj_objectviewing_adaptation_analysis_dir);
    if ~isnan(subjects_file_data{subj, 21}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 21}, '.mat']); copyfile(curr_file, curr_subj_jrd_analysis_dirs{1}); end % JRD runs
    if ~isnan(subjects_file_data{subj, 22}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 22}, '.mat']); copyfile(curr_file, curr_subj_jrd_analysis_dirs{2}); end
    if ~isnan(subjects_file_data{subj, 23}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 23}, '.mat']); copyfile(curr_file, curr_subj_jrd_analysis_dirs{3}); end
    if ~isnan(subjects_file_data{subj, 21}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 21}, '_adaptation.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_analysis_dir); end  % JRD adaptation files
    if ~isnan(subjects_file_data{subj, 22}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 22}, '_adaptation.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_analysis_dir); end
    if ~isnan(subjects_file_data{subj, 23}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 23}, '_adaptation.mat']); copyfile(curr_file, curr_subj_jrd_adaptation_analysis_dir); end
    if ~isnan(subjects_file_data{subj, 21}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 21}, '_univariate.mat']); copyfile(curr_file, curr_subj_jrd_univariate_analysis_dir); end  % JRD univariate analysis files
    if ~isnan(subjects_file_data{subj, 22}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 22}, '_univariate.mat']); copyfile(curr_file, curr_subj_jrd_univariate_analysis_dir); end
    if ~isnan(subjects_file_data{subj, 23}), curr_file = fullfile(subjects_file_data{subj, 27}, [subjects_file_data{subj, 23}, '_univariate.mat']); copyfile(curr_file, curr_subj_jrd_univariate_analysis_dir); end
    jrd_info_files = dir(fullfile(subjects_file_data{subj, 27},'all_info_*.mat')); for i = 1:length(jrd_info_files), copyfile(fullfile(jrd_info_files(i).folder, jrd_info_files(i).name), curr_subj_jrd_analysis_dirs{i}); end
    
    % Extracting files from gz format to nii
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs{subj});
    curr_subj_session_1_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day1/func']);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);
    session1_gz_files = dir(fullfile(curr_subj_session_1_dir, '*.gz'));
    for i = 1:length(session1_gz_files)
        curr_gz_file = fullfile(curr_subj_session_1_dir, session1_gz_files(i).name);
        if ~exist(curr_gz_file(1:end-3), 'file')
            gunzip(fullfile(curr_subj_session_1_dir, session1_gz_files(i).name), curr_subj_session_1_dir);
            delete(fullfile(curr_subj_session_1_dir, session1_gz_files(i).name), curr_subj_session_1_dir);
        end
    end
    session2_gz_files = dir(fullfile(curr_subj_session_2_dir, '*.gz'));
    for i = 1:length(session2_gz_files)
        curr_gz_file = fullfile(curr_subj_session_2_dir, session2_gz_files(i).name);
        if ~exist(curr_gz_file(1:end-3), 'file')
            gunzip(fullfile(curr_subj_session_2_dir, session2_gz_files(i).name), curr_subj_session_2_dir);
            delete(fullfile(curr_subj_session_2_dir, session2_gz_files(i).name), curr_subj_session_2_dir);
        end
    end
    
    % Reading the confounds files and saving them in the right format
    % Functional localizer
    localizer_confound_files = dir(fullfile(curr_subj_session_1_dir, '*func_localizer*confounds*.tsv'));
    for f = 1:length(localizer_confound_files)
        conf = tdfread(fullfile(curr_subj_session_1_dir, localizer_confound_files(f).name));
        fd = zeros(length(conf.csf), 1); for i = 2:length(conf.csf), fd(i) = str2double(conf.framewise_displacement(i, :)); end     % Converting framewise displacement from string to numbers
        R = [conf.trans_x, conf.trans_y, conf.trans_z, conf.rot_x, conf.rot_y, conf.rot_z, fd, conf.csf, conf.white_matter];
        save(fullfile(curr_subj_localizer_analysis_dir, ['localizer', num2str(f), '_confounds.mat']), 'R', 'names');
    end
    % Object viewing
    object_viewing_confound_files = dir(fullfile(curr_subj_session_1_dir, '*object_viewing*confounds*.tsv'));
    object_viewing_confound_files = [object_viewing_confound_files; dir(fullfile(curr_subj_session_2_dir, '*object_viewing*confounds*.tsv'))];
    for f = 1:length(object_viewing_confound_files)
        conf = tdfread(fullfile(object_viewing_confound_files(f).folder, object_viewing_confound_files(f).name));
        fd = zeros(length(conf.csf), 1); for i = 2:length(conf.csf), fd(i) = str2double(conf.framewise_displacement(i, :)); end     % Converting framewise displacement from string to numbers
        R = [conf.trans_x, conf.trans_y, conf.trans_z, conf.rot_x, conf.rot_y, conf.rot_z, fd, conf.csf, conf.white_matter];
        save(fullfile(curr_subj_objectviewing_analysis_dirs{f}, ['object_viewing', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, ['object_viewing', num2str(f), '_confounds.mat']), 'R', 'names');
    end
    % JRD
    jrd_confound_files = dir(fullfile(curr_subj_session_2_dir, '*jrd*confounds*.tsv'));
    for f = 1:length(jrd_confound_files)
        conf = tdfread(fullfile(curr_subj_session_2_dir, jrd_confound_files(f).name));
        fd = zeros(length(conf.csf), 1); for i = 2:length(conf.csf), fd(i) = str2double(conf.framewise_displacement(i, :)); end     % Converting framewise displacement from string to numbers
        R = [conf.trans_x, conf.trans_y, conf.trans_z, conf.rot_x, conf.rot_y, conf.rot_z, fd, conf.csf, conf.white_matter];
        save(fullfile(curr_subj_jrd_analysis_dirs{f}, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_jrd_adaptation_analysis_dir, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
        save(fullfile(curr_subj_jrd_univariate_analysis_dir, ['jrd', num2str(f), '_confounds.mat']), 'R', 'names');
    end
    
end


%% Smoothing the functional data
spm('defaults', 'FMRI');
for subj = 1:length(subj_data_dirs)
    disp(subj)
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs{subj});
    curr_subj_analysis_dir = fullfile(curr_subj_parent_dir, 'Analysis');
    curr_subj_session_1_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day1/func']);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);
    
    % Defining the spm batch for running a GLM
    matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'sm_';
    
    % Identifying all of the functional runs
    all_func_data_files = dir(fullfile(curr_subj_session_1_dir, '*preproc_bold.nii'));
    all_func_data_files = [all_func_data_files; dir(fullfile(curr_subj_session_2_dir, '*preproc_bold.nii'))];
    all_func_data_files = all_func_data_files(cellfun(@isempty, strfind({all_func_data_files.name},'sm_')));    % Removing already-smoothed files from the list

    for r = 1:length(all_func_data_files)
        n = nifti(fullfile(all_func_data_files(r).folder, all_func_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;     % Finding out the number of runs
        matlabbatch{1}.spm.spatial.smooth.data = cell(num_volumes, 1);     % The actual data
        for v = 1:num_volumes
            matlabbatch{1}.spm.spatial.smooth.data{v} = [fullfile(all_func_data_files(r).folder, all_func_data_files(r).name), ',', num2str(v)];
        end
        % Checking if the smoothed data already exists, running smoothing if not
        if ~exist(fullfile(all_func_data_files(r).folder, ['sm_', all_func_data_files(r).name]), 'file')
            spm_jobman('run', matlabbatch);
        end
    end

    clear matlabbatch;
end




%% Running the GLM on each run (estimating betas)
spm('defaults', 'FMRI');
for subj = 1:length(subj_data_dirs)
    disp(subj)
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs{subj});
    curr_subj_analysis_dir = fullfile(curr_subj_parent_dir, 'Analysis');
    curr_subj_session_1_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day1/func']);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);
    
    % Defining the spm batch for running a GLM
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % Defining design and running GLM for object viewing runs (separately for each run)
    object_viewing_data_files = dir(fullfile(curr_subj_session_1_dir, 'sm_*object_viewing*preproc_bold.nii'));
    object_viewing_data_files = [object_viewing_data_files; dir(fullfile(curr_subj_session_2_dir, 'sm_*object_viewing*preproc_bold.nii'))];
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objectviewing_analysis_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run4')};
    for r = 1:length(curr_subj_objectviewing_analysis_dirs)     % Defining specific parameters for each run
        design_file = dir(fullfile(curr_subj_objectviewing_analysis_dirs{r}, '*Segmentation*.mat'));
        confounds_file = dir(fullfile(curr_subj_objectviewing_analysis_dirs{r}, '*confounds.mat'));
        n = nifti(fullfile(object_viewing_data_files(r).folder, object_viewing_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
        matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_objectviewing_analysis_dirs{r}};    % The output directory
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(curr_subj_objectviewing_analysis_dirs{r}, design_file(1).name)};      % The design file
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(curr_subj_objectviewing_analysis_dirs{r}, confounds_file(1).name)};      % The confounds file
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cell(num_volumes, 1);     % The actual data
        for v = 1:num_volumes
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans{v} = [fullfile(object_viewing_data_files(r).folder, object_viewing_data_files(r).name), ',', num2str(v)];
        end
        if ~exist(fullfile(curr_subj_objectviewing_analysis_dirs{r}, 'SPM.mat'), 'file')        % Checking if SPM.mat file already exists, running if not
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_objectviewing_analysis_dirs{r}, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
            % Running the design creation and estimation script
            spm_jobman('run', matlabbatch);
        end
    end
    
    % Defining design and running GLM for JRD runs (separately for each run)
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    for r = 1:length(curr_subj_jrd_analysis_dirs)     % Defining specific parameters for each run
        design_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, '*Segmentation*.mat'));
        if ~isempty(design_file)
            confounds_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, '*confounds.mat'));
            n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_jrd_analysis_dirs{r}};    % The output directory
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(curr_subj_jrd_analysis_dirs{r}, design_file(1).name)};      % The design file
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(curr_subj_jrd_analysis_dirs{r}, confounds_file(1).name)};      % The confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cell(num_volumes, 1);     % The actual data
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
            end
            if ~exist(fullfile(curr_subj_jrd_analysis_dirs{r}, 'SPM.mat'), 'file')        % Checking if SPM.mat file already exists, running if not
                matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_jrd_analysis_dirs{r}, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
                % Running the design creation and estimation script
                spm_jobman('run', matlabbatch);
            end
        end
    end
    
    % Defining design and running GLM for Functional localizer runs
    localizer_data_files = dir(fullfile(curr_subj_session_1_dir, 'sm_*localizer*preproc_bold.nii'));
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    matlabbatch{1}.spm.stats.fmri_spec.sess = struct();
    design_file = dir(fullfile(curr_subj_localizer_analysis_dir, '*psychopy*.mat'));
    confounds_file = dir(fullfile(curr_subj_localizer_analysis_dir, '*confounds.mat'));
    for r = 1:2     % Defining specific parameters for each of the two localizer runs
        n = nifti(fullfile(localizer_data_files(r).folder, localizer_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
        matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_localizer_analysis_dir};    % The output directory
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_localizer_analysis_dir, design_file(r).name)};      % The design file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_localizer_analysis_dir, confounds_file(r).name)};      % The confounds file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
        for v = 1:num_volumes
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(localizer_data_files(r).folder, localizer_data_files(r).name), ',', num2str(v)];
        end
    end
    if ~exist(fullfile(curr_subj_localizer_analysis_dir, 'SPM.mat'), 'file')
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_localizer_analysis_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
        % Running the design creation and estimation script
        spm_jobman('run', matlabbatch);
    end
    
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
    object_viewing_adaptation_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*object_viewing*preproc_bold.nii'));
    curr_subj_objectviewing_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    matlabbatch{1}.spm.stats.fmri_spec.sess = struct();
    design_file = dir(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, '*Segmentation*.mat'));
    confounds_file = dir(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, '*confounds.mat'));
    for r = 1:2     % Defining specific parameters for each run
        n = nifti(fullfile(object_viewing_adaptation_data_files(r).folder, object_viewing_adaptation_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
        matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_objectviewing_adaptation_analysis_dir};    % The output directory
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_objectviewing_adaptation_analysis_dir, design_file(r+2).name)};      % The design file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_objectviewing_adaptation_analysis_dir, confounds_file(r+2).name)};      % The confounds file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
        for v = 1:num_volumes
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(object_viewing_adaptation_data_files(r).folder, object_viewing_adaptation_data_files(r).name), ',', num2str(v)];
        end
    end
    if ~exist(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, 'SPM.mat'), 'file')        % Checking if SPM.mat file already exists, running if not
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_objectviewing_adaptation_analysis_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
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
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    design_file = dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*Segmentation*.mat'));
    confounds_file = dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat'));
    exist_jrd = 0;
    for r = 1:length(jrd_data_files)     % Defining specific parameters for each run
        if exist(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name))
            exist_jrd = 1;
            n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_jrd_adaptation_analysis_dir};    % The output directory
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_jrd_adaptation_analysis_dir, design_file(r).name)};      % The design file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_jrd_adaptation_analysis_dir, confounds_file(r).name)};      % The confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
            end
        end
    end
    if ~exist(fullfile(curr_subj_jrd_adaptation_analysis_dir, 'SPM.mat'), 'file') & exist_jrd == 1        % Checking if SPM.mat file already exists, running if not
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_jrd_adaptation_analysis_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
        % Running the design creation and estimation script
        spm_jobman('run', matlabbatch);
    end
    
    % Defining design and running GLM for JRD runs - univariate analysis
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
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    curr_subj_jrd_univariate_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_univariate');
    design_file = dir(fullfile(curr_subj_jrd_univariate_analysis_dir, '*Segmentation*.mat'));
    confounds_file = dir(fullfile(curr_subj_jrd_univariate_analysis_dir, '*confounds.mat'));
    exist_jrd = 0;
    for r = 1:length(jrd_data_files)     % Defining specific parameters for each run
        if exist(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name))
            exist_jrd = 1;
            n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_jrd_univariate_analysis_dir};    % The output directory
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_jrd_univariate_analysis_dir, design_file(r).name)};      % The design file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_jrd_univariate_analysis_dir, confounds_file(r).name)};      % The confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
            end
        end
    end
    if ~exist(fullfile(curr_subj_jrd_univariate_analysis_dir, 'SPM.mat'), 'file') & exist_jrd == 1        % Checking if SPM.mat file already exists, running if not
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_jrd_univariate_analysis_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
        % Running the design creation and estimation script
        spm_jobman('run', matlabbatch);
    end
    
    
    % Single trial beta estimates (Mumford 2012 NeuroImage - each trial as an independent GLM with it as one regressor and the rest of the trials as another)
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
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    for r = 1:length(curr_subj_jrd_analysis_dirs)
        % Defining and creating directories
        curr_subj_jrd_singletrials_analysis_dirs{r} = fullfile(curr_subj_jrd_analysis_dirs{r}, 'single_trials');
        curr_subj_jrd_singletrials_temp_dirs{r} = fullfile(curr_subj_jrd_singletrials_analysis_dirs{r}, 'temp');
        if ~exist(curr_subj_jrd_singletrials_analysis_dirs{r}, 'dir')
            mkdir(curr_subj_jrd_singletrials_analysis_dirs{r});
            mkdir(curr_subj_jrd_singletrials_temp_dirs{r});
        end
        
        % Reading the design file for this run
        design_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, '*Segmentation*.mat'));
        if ~isempty(design_file)
            design_file = load(fullfile(design_file(1).folder, design_file(1).name));
            [all_onsets, indices] = sort(cell2mat(design_file.onsets));                                     % Sorting by the order of appearance
            all_durations = cell2mat(design_file.durations); all_durations = all_durations(indices);        % Sorting by the order of the onsets sort
            % Defining general SPM batch parameters
            matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_jrd_singletrials_temp_dirs{r}};    % The output directory
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_jrd_singletrials_temp_dirs{r}, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
            % Entering the functional data filenames
            n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cell(num_volumes, 1);
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
            end
            % Running the estimation for each run separately
            for i = 1:length(all_onsets)
                if exist(fullfile(curr_subj_jrd_singletrials_temp_dirs{r}, 'SPM.mat')), delete(fullfile(curr_subj_jrd_singletrials_temp_dirs{r}, 'SPM.mat')); end
                % Defining all of the regressors - one for the specific trial, the other for all the rest
                onsets = {}; durations = {};
                onsets{1} = all_onsets(i);
                onsets{2} = all_onsets; onsets{2}(i) = [];    % A regressor for all the other trials
                durations{1} = all_durations(i);     % A regressor only with this specific trial
                durations{2} = all_durations; durations{2}(i) = [];
                names = {'current_trial','all_other_trials'};
                % Saving the design file
                temp_design_file_name = fullfile(curr_subj_jrd_singletrials_temp_dirs{r}, 'temp_design.mat');
                save(temp_design_file_name, 'names','onsets','durations');
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {temp_design_file_name};      % Defining the design file in the SPM matlab batch
                % Running the design creation and estimation script
                spm_jobman('run', matlabbatch);
                % Copying the resulting beta1 (the single-trial estimate) to the parent directory
                copyfile(fullfile(curr_subj_jrd_singletrials_temp_dirs{r}, 'beta_0001.nii'), fullfile(curr_subj_jrd_singletrials_analysis_dirs{r}, ['singletrial_beta_time_',num2str(onsets{1},'%03.0f'),'.nii']));
            end
        end
    end
    
    
    
    % Defining design and running GLM for object viewing runs - combined runs
    clear matlabbatch;
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
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
    % Day 1 GLM
    curr_subj_objectviewing_analysis_dirs_day1 = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2')};
    curr_subj_objectviewing_combined_day1_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day1');
    matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_objectviewing_combined_day1_dir};    % The output directory
    object_viewing_data_files = dir(fullfile(curr_subj_session_1_dir, 'sm_*object_viewing*preproc_bold.nii'));
    for r = 1:length(object_viewing_data_files)     % Defining specific parameters for each run in day1
        design_file = dir(fullfile(curr_subj_objectviewing_analysis_dirs_day1{r}, '*Segmentation*.mat'));
        confounds_file = dir(fullfile(curr_subj_objectviewing_analysis_dirs_day1{r}, '*confounds.mat'));
        n = nifti(fullfile(object_viewing_data_files(r).folder, object_viewing_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_objectviewing_analysis_dirs_day1{r}, design_file(1).name)};      % The design file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_objectviewing_analysis_dirs_day1{r}, confounds_file(1).name)};      % The confounds file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
        for v = 1:num_volumes
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(object_viewing_data_files(r).folder, object_viewing_data_files(r).name), ',', num2str(v)];
        end
    end
    if ~exist(fullfile(curr_subj_objectviewing_combined_day1_dir, 'SPM.mat'), 'file')        % Checking if SPM.mat file already exists, running if not
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_objectviewing_combined_day1_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
        % Running the design creation and estimation script
        spm_jobman('run', matlabbatch);
    end
    % Day 2 GLM
    curr_subj_objectviewing_analysis_dirs_day2 = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run4')};
    curr_subj_objectviewing_combined_day2_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day2');
    matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_objectviewing_combined_day2_dir};    % The output directory
    object_viewing_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*object_viewing*preproc_bold.nii'));
    for r = 1:length(object_viewing_data_files)     % Defining specific parameters for each run in day2
        design_file = dir(fullfile(curr_subj_objectviewing_analysis_dirs_day2{r}, '*Segmentation*.mat'));
        confounds_file = dir(fullfile(curr_subj_objectviewing_analysis_dirs_day2{r}, '*confounds.mat'));
        n = nifti(fullfile(object_viewing_data_files(r).folder, object_viewing_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_objectviewing_analysis_dirs_day2{r}, design_file(1).name)};      % The design file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_objectviewing_analysis_dirs_day2{r}, confounds_file(1).name)};      % The confounds file
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
        for v = 1:num_volumes
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(object_viewing_data_files(r).folder, object_viewing_data_files(r).name), ',', num2str(v)];
        end
    end
    if ~exist(fullfile(curr_subj_objectviewing_combined_day2_dir, 'SPM.mat'), 'file')        % Checking if SPM.mat file already exists, running if not
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_objectviewing_combined_day2_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
        % Running the design creation and estimation script
        spm_jobman('run', matlabbatch);
    end
    
    % Defining design and running GLM for JRD runs - combined runs
    clear matlabbatch;
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
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
    matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_jrd_combined_dir};    % The output directory
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    exist_jrd = 0;
    for r = 1:length(curr_subj_jrd_analysis_dirs)     % Defining specific parameters for each run
        design_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, '*Segmentation*.mat'));
        if ~isempty(design_file)
            exist_jrd = 1;
            confounds_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, '*confounds.mat'));
            n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {fullfile(curr_subj_jrd_analysis_dirs{r}, design_file(1).name)};      % The design file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(curr_subj_jrd_analysis_dirs{r}, confounds_file(1).name)};      % The confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cell(num_volumes, 1);     % The actual data
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
            end
        end
    end
    if ~exist(fullfile(curr_subj_jrd_combined_dir, 'SPM.mat'), 'file') & exist_jrd == 1       % Checking if SPM.mat file already exists, running if not
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_jrd_combined_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
        % Running the design creation and estimation script
        spm_jobman('run', matlabbatch);
    end

    clear matlabbatch;
end


%% Creating contrasts and computing T-values
for subj = 1:length(subj_data_dirs)
    disp(subj)
    
    % Defining directories
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'Analysis');
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objectviewing_analysis_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run4')};
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    curr_subj_objectviewing_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    curr_subj_jrd_univariate_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_univariate');
    curr_subj_objectviewing_combined_day1_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day1');
    curr_subj_objectviewing_combined_day2_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day2');
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');

    % Functional localizer contrasts
    SPMmat_filename = fullfile(curr_subj_localizer_analysis_dir, 'SPM.mat');
    % Deleting existing contrasts
    load(SPMmat_filename);
    if ~isempty(SPM.xCon)
        SPM.xCon = []; save(SPMmat_filename, 'SPM');
        temp_con_files = dir(fullfile(curr_subj_localizer_analysis_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_localizer_analysis_dir, temp_con_files(i).name)); end
        temp_t_files = dir(fullfile(curr_subj_localizer_analysis_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_localizer_analysis_dir, temp_t_files(i).name)); end
    end
    % Creating the new contrasts and T images
    matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'places vs. objects';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1, 0, -1, 0, zeros(1, num_confounds), 1, 0, -1, 0, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'places vs. everything';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [3, -1, -1, -1, zeros(1, num_confounds), 3, -1, -1, -1, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'objects vs. faces';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0, -1, 1, 0, zeros(1, num_confounds), 0, -1, 1, 0, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'objects vs. everything';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [-1, -1, 3, -1, zeros(1, num_confounds), -1, -1, 3, -1, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    
    % Object viewing - one contrast for each object (to get the t-value instead of beta)
    for r = 1:length(curr_subj_objectviewing_analysis_dirs)
        SPMmat_filename = fullfile(curr_subj_objectviewing_analysis_dirs{r}, 'SPM.mat');
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_objectviewing_analysis_dirs{r}, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objectviewing_analysis_dirs{r}, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_objectviewing_analysis_dirs{r}, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objectviewing_analysis_dirs{r}, temp_t_files(i).name)); end
        end
        % Creating the new contrasts and T images
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        for o = 1:num_betas
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.name = ['object', num2str(o)];
            current_contrast = zeros(1, num_betas); current_contrast(o) = 1;
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.weights = [current_contrast, zeros(1, num_confounds + 1)];
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.sessrep = 'none';
        end
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end
    
    % JRD - one contrast for each object (to get the t-value)
    for r = 1:length(curr_subj_jrd_analysis_dirs)
        SPMmat_filename = fullfile(curr_subj_jrd_analysis_dirs{r}, 'SPM.mat');
        if exist(SPMmat_filename)
            % Deleting existing contrasts
            load(SPMmat_filename);
            if ~isempty(SPM.xCon)
                SPM.xCon = []; save(SPMmat_filename, 'SPM');
                temp_con_files = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_analysis_dirs{r}, temp_con_files(i).name)); end
                temp_t_files = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_analysis_dirs{r}, temp_t_files(i).name)); end
            end
            % Creating the new contrasts and T images
            matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
            matlabbatch{1}.spm.stats.con.delete = 0;
            for o = 1:num_betas
                matlabbatch{1}.spm.stats.con.consess{o}.tcon.name = ['object', num2str(o)];
                current_contrast = zeros(1, num_betas); current_contrast(o) = 1;
                matlabbatch{1}.spm.stats.con.consess{o}.tcon.weights = [current_contrast, zeros(1, num_confounds + 1)];
                matlabbatch{1}.spm.stats.con.consess{o}.tcon.sessrep = 'none';
            end
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
        end
    end
    
    % Object viewing adaptation contrasts
    SPMmat_filename = fullfile(curr_subj_objectviewing_adaptation_analysis_dir, 'SPM.mat');
    % Deleting existing contrasts
    load(SPMmat_filename);
    if ~isempty(SPM.xCon)
        SPM.xCon = []; save(SPMmat_filename, 'SPM');
        temp_con_files = dir(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, temp_con_files(i).name)); end
        temp_t_files = dir(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objectviewing_adaptation_analysis_dir, temp_t_files(i).name)); end
    end
    % Creating the new contrasts and T images
    matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'All stimuli regressors vs rest';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0, 1, 1, 1, 0, zeros(1, num_confounds), 0, 1, 1, 1, 0, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Different segment > same segment';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1, 0, 0, 0, 0, zeros(1, num_confounds), 1, 0, 0, 0, 0, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Distance effect parametric modulation';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0, 0, 0, 0, 1, zeros(1, num_confounds), 0, 0, 0, 0, 1, zeros(1, num_confounds + 2)];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % JRD adaptation contrasts
    SPMmat_filename = fullfile(curr_subj_jrd_adaptation_analysis_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_adaptation_analysis_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_adaptation_analysis_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'All stimuli regressors vs rest';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([0, 1, 1, 0, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Different segment > same segment';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([1, 0, 0, 0, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Distance effect parametric modulation';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [repmat([0, 0, 0, 1, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end

    % JRD univariate contrasts
    SPMmat_filename = fullfile(curr_subj_jrd_univariate_analysis_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_univariate_analysis_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_univariate_analysis_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_univariate_analysis_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_univariate_analysis_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_univariate_analysis_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'All stimuli regressors vs rest';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [repmat([1, 1, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Between_segments > within segments';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [repmat([-1, 1, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Within segments vs rest';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [repmat([1, 0, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Between_segments vs rest';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [repmat([0, 1, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end
    
    
    % Object viewing combined runs contrasts - one per object, combined across runs
    % Day 1
    SPMmat_filename = fullfile(curr_subj_objectviewing_combined_day1_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_objectviewing_combined_day1_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objectviewing_combined_day1_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_objectviewing_combined_day1_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objectviewing_combined_day1_dir, temp_t_files(i).name)); end
        end
        % Creating the new contrasts and T images
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        for o = 1:num_betas
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.name = ['object', num2str(o)];
            current_contrast = zeros(1, num_betas); current_contrast(o) = 1;
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.weights = [current_contrast, zeros(1, num_confounds), current_contrast, zeros(1, num_confounds), 0, 0];
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.sessrep = 'none';
        end
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end
    % Day2
    SPMmat_filename = fullfile(curr_subj_objectviewing_combined_day2_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_objectviewing_combined_day2_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_objectviewing_combined_day2_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_objectviewing_combined_day2_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_objectviewing_combined_day2_dir, temp_t_files(i).name)); end
        end
        % Creating the new contrasts and T images
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        for o = 1:num_betas
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.name = ['object', num2str(o)];
            current_contrast = zeros(1, num_betas); current_contrast(o) = 1;
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.weights = [current_contrast, zeros(1, num_confounds), current_contrast, zeros(1, num_confounds), 0, 0];
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.sessrep = 'none';
        end
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end

    % JRD combined runs contrasts - one per object, combined across runs
    SPMmat_filename = fullfile(curr_subj_jrd_combined_dir, 'SPM.mat');
    if exist(SPMmat_filename)
        % Deleting existing contrasts
        load(SPMmat_filename);
        if ~isempty(SPM.xCon)
            SPM.xCon = []; save(SPMmat_filename, 'SPM');
            temp_con_files = dir(fullfile(curr_subj_jrd_combined_dir, 'con*.nii')); for i = 1:length(temp_con_files), delete(fullfile(curr_subj_jrd_combined_dir, temp_con_files(i).name)); end
            temp_t_files = dir(fullfile(curr_subj_jrd_combined_dir, 'spmT*.nii')); for i = 1:length(temp_t_files), delete(fullfile(curr_subj_jrd_combined_dir, temp_t_files(i).name)); end
        end
        curr_num_JRD_sessions = length(dir(fullfile(curr_subj_jrd_univariate_analysis_dir, '*confounds.mat')));
        % Creating the new contrasts and T images
        matlabbatch{1}.spm.stats.con.spmmat = {SPMmat_filename};
        matlabbatch{1}.spm.stats.con.delete = 0;
        for o = 1:num_betas
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.name = ['object', num2str(o)];
            current_contrast = zeros(1, num_betas); current_contrast(o) = 1;
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.weights = [repmat([current_contrast, zeros(1, num_confounds)], 1, curr_num_JRD_sessions), 0, 0];
            matlabbatch{1}.spm.stats.con.consess{o}.tcon.sessrep = 'none';
        end
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end

end





%% Defining subject specific ROIs
% Read group-level ROIs (Julian et al. 2012)
% betas_file1 = '/Users/mpeer/Desktop/Segmentation_subjects_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
lRSC = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/lRSC.nii', betas_file1, 0);
rRSC = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/rRSC.nii', betas_file1, 0);
lPPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/lPPA.nii', betas_file1, 0);
rPPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/rPPA.nii', betas_file1, 0);
lOPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/lTOS.nii', betas_file1, 0);
rOPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/rTOS.nii', betas_file1, 0);
bilat_RSC = lRSC + rRSC; bilat_PPA = lPPA + rPPA; bilat_OPA = lOPA + rOPA;
all_ROIs = {lRSC, rRSC, lPPA, rPPA, lOPA, rOPA};
localizer_contrasts_to_use = {'2','2','2','2','2','2'};         % Which localizer contrast to use for each ROI - 1 is places vs. objects, 2 places vs. everything, 3 objects vs. faces, 4 objects vs. everything

% Defining ROIs for each subject, based on anatomical parcellation and functional localizer contrasts
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');

    % Defining subject specific functional ROIs (from the localizer run, masked with the Julian et al. ROIs) and combining bilaterally
    all_ROIs_individual_temp = all_ROIs;
    for i = 1:length(all_ROIs)
        localizer_contrast = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, ['spmT_000',localizer_contrasts_to_use{i},'.nii'])));
        current_threshold = min(maxk(localizer_contrast(all_ROIs{i} == 1), 100));       % Choosing the 100 voxels with the maximal t value inside the ROI
        all_ROIs_individual_temp{i} = all_ROIs{i} & (localizer_contrast >= current_threshold);     % masking the first ROIs with the individual localizer images
    end
    % Combining ROIs bilaterally
    all_ROIs_individual = {};
    for i = 1:2:length(all_ROIs_individual_temp)
        all_ROIs_individual{(i+1)/2} = all_ROIs_individual_temp{i} + all_ROIs_individual_temp{i+1};
    end
    % Adding unilateral ROIs
    all_ROIs_individual = [all_ROIs_individual, all_ROIs_individual_temp];
    
    % Getting hippocampal and entorhinal anatomical ROIs of this subject from freesurfer parcellation
    curr_subject_functional_directory = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'fmriprep');
    curr_subject_functional_directory = fullfile(fullfile(curr_subject_functional_directory, subj_data_dirs{subj}), 'ses-day1/func/');
    curr_parcellation_file = dir(fullfile(curr_subject_functional_directory, '*aparcaseg*.nii'));
    curr_parcellation = reslice_data(fullfile(curr_subject_functional_directory, curr_parcellation_file(1).name), betas_file1, 0);
    left_HC = curr_parcellation == 17; right_HC = curr_parcellation == 53; bilat_HC = left_HC + right_HC;
    left_EC = curr_parcellation == 1006; right_EC = curr_parcellation == 2006; bilat_EC = left_EC + right_EC;
    all_ROIs_individual = [left_EC, right_EC, bilat_EC, left_HC, right_HC, bilat_HC, all_ROIs_individual];
    
    % Saving the individual ROIs to files
    all_ROI_names = {'lEC', 'rEC', 'bilat_EC', 'lHC', 'rHC', 'bilat_HC', 'bilat_RSC', 'bilat_PPA', 'bilat_OPA', 'lRSC', 'rRSC', 'lPPA', 'rPPA', 'lOPA', 'rOPA'};
    for i = 1:length(all_ROIs_individual)
        save_mat_to_nifti(betas_file1, all_ROIs_individual{i}, fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii']));
    end
end

% Defining ROIs for each subject, based on maximum activation in tasks
all_ROI_names = {'lRSC', 'rRSC', 'lPPA', 'rPPA', 'lOPA', 'rOPA', 'lHC', 'rHC'};
for subj = 1:length(subj_data_dirs)
    disp(subj)
    
    % Defining current subject directories
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objectviewing_day2_combined_dir = fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day2');
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');
    
    % Getting hippocampal and entorhinal anatomical ROIs of this subject from freesurfer parcellation
    curr_subject_functional_directory = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'fmriprep');
    curr_subject_functional_directory = fullfile(fullfile(curr_subject_functional_directory, subj_data_dirs{subj}), 'ses-day1/func/');
    curr_parcellation_file = dir(fullfile(curr_subject_functional_directory, '*aparcaseg*.nii'));
    curr_parcellation = reslice_data(fullfile(curr_subject_functional_directory, curr_parcellation_file(1).name), betas_file1, 0);
    left_HC = curr_parcellation == 17; right_HC = curr_parcellation == 53; bilat_HC = left_HC + right_HC;
    all_ROIs = {lRSC, rRSC, lPPA, rPPA, lOPA, rOPA, left_HC, right_HC};
    
    temp_ROIs_jrd = cell(size(all_ROIs));
    % Reading the 100 most active voxels in JRD in each hemisphere and roi, and combining bilaterally
    if exist(fullfile(curr_subj_jrd_combined_dir, 'spmT_0001.nii'),'file')
        % Getting the t-values for all objects and averaging them
        all_tvalues_ROI_jrd = nan(length(all_ROIs{1}(:)), num_betas);
        t_values_files_current = dir(fullfile(curr_subj_jrd_combined_dir, 'spmT*.nii'));
        for b = 1:num_betas
            curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_jrd_combined_dir, t_values_files_current(b).name)));
            all_tvalues_ROI_jrd(:,b) = curr_tvalues_image(:);
        end
        jrd_contrast = nanmean(all_tvalues_ROI_jrd, 2);
        jrd_contrast = reshape(jrd_contrast, size(all_ROIs{1}));
        % Defining the ROIs based on the top 100 voxels
        for i=1:length(all_ROIs)
            current_threshold = min(maxk(jrd_contrast(all_ROIs{i} == 1), 100));       % Choosing the 100 voxels with the maximal t value inside the ROI
            temp_ROIs_jrd{i} = all_ROIs{i} & (jrd_contrast >= current_threshold);     % masking the first ROIs with the individual localizer images
        end
    else
        for i=1:length(all_ROIs)
            temp_ROIs_jrd{i} = nan(size(all_ROIs{1}));
        end
    end
    % Reading the 100 most active voxels in object-viewing day2 in each hemisphere and roi, and combining bilaterally
    temp_ROIs_objview = cell(size(all_ROIs));
    if exist(fullfile(curr_subj_objectviewing_day2_combined_dir, 'spmT_0001.nii'),'file')
        % Getting the t-values for all objects and averaging them
        all_tvalues_ROI_objview = nan(length(all_ROIs{1}(:)), num_betas);
        t_values_files_current = dir(fullfile(curr_subj_objectviewing_day2_combined_dir, 'spmT*.nii'));
        for b = 1:num_betas
            curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_objectviewing_day2_combined_dir, t_values_files_current(b).name)));
            all_tvalues_ROI_objview(:,b) = curr_tvalues_image(:);
        end
        objectviewing_contrast = nanmean(all_tvalues_ROI_objview, 2);
        objectviewing_contrast = reshape(objectviewing_contrast, size(all_ROIs{1}));
        % Defining the ROIs based on the top 100 voxels
        for i=1:length(all_ROIs)
            current_threshold = min(maxk(objectviewing_contrast(all_ROIs{i} == 1), 100));       % Choosing the 100 voxels with the maximal t value inside the ROI
            temp_ROIs_objview{i} = all_ROIs{i} & (objectviewing_contrast >= current_threshold);     % masking the first ROIs with the individual localizer images
        end        
    else
        for i=1:length(all_ROIs)
            temp_ROIs_objview{i} = nan(size(all_ROIs{1}));
        end
    end
    
    % Saving the per-subject ROIs
    for i = 1:length(all_ROI_names)
        % Saving the per-subject ROI
        save_mat_to_nifti(betas_file1, temp_ROIs_jrd{i}, fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_activation_jrd_individual.nii']));
        save_mat_to_nifti(betas_file1, temp_ROIs_objview{i}, fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_activation_objview_individual.nii']));
    end
end
% % Combining the bilateral activation-based ROIs and saving
new_ROI_names = {'RSC','OPA','PPA','HC'};
for subj = 1:length(subj_data_dirs)
    disp(subj)

    % Defining current subject directories
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');

    for i = 1:length(new_ROI_names)
        curr_l_ROI_jrd = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, ['l', new_ROI_names{i}, '_activation_jrd_individual.nii'])));
        curr_r_ROI_jrd = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, ['r', new_ROI_names{i}, '_activation_jrd_individual.nii'])));
        curr_l_ROI_objview = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, ['l', new_ROI_names{i}, '_activation_objview_individual.nii'])));
        curr_r_ROI_objview = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, ['r', new_ROI_names{i}, '_activation_objview_individual.nii'])));

        curr_bilat_ROI_jrd = curr_l_ROI_jrd + curr_r_ROI_jrd;
        curr_bilat_ROI_objview = curr_l_ROI_objview + curr_r_ROI_objview;

        save_mat_to_nifti(betas_file1, curr_bilat_ROI_jrd, fullfile(curr_subj_localizer_analysis_dir, ['bilat_', new_ROI_names{i}, '_activation_jrd_individual.nii']));
        save_mat_to_nifti(betas_file1, curr_bilat_ROI_objview, fullfile(curr_subj_localizer_analysis_dir, ['bilat_', new_ROI_names{i}, '_activation_objview_individual.nii']));
    end
end

