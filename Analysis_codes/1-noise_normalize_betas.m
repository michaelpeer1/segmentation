% Perform noise normalization on the data
% The actual noise Normalization uses code from the Decoding Toolbox (TDT - Hebart, Gorgen, Haynes 2015)

betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

%% Go over subjects
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    % Object viewing analysis dirs
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objectviewing_analysis_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run4')};
    num_runs_objectviewing = length(curr_subj_objectviewing_analysis_dirs);
    curr_subj_objectviewing_combined_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day2')};
    % JRD analysis dirs
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    num_runs_jrd = length(dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat')));
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');
    % Functional localizer analysis dirs
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    
    %     all_analysis_dirs = [curr_subj_objectviewing_analysis_dirs, curr_subj_objectviewing_combined_dirs, curr_subj_jrd_analysis_dirs, curr_subj_jrd_combined_dir];
    all_analysis_dirs = [curr_subj_objectviewing_analysis_dirs, curr_subj_jrd_analysis_dirs];
    
    % Go over analysis directories
    for d = 1:length(all_analysis_dirs)
        curr_analysis_dir = all_analysis_dirs{d};
        
        % Read SPM file
        if exist(fullfile(curr_analysis_dir, 'SPM.mat'),'file')
            SPM = load(fullfile(curr_analysis_dir, 'SPM.mat'));
            SPM = SPM.SPM;
            
            % Change filenames to correct path
            for f = 1:length(SPM.xY.VY)
                curr_filename = SPM.xY.VY(f).fname;
                loc_start_subj_filename = strfind(curr_filename, 'sub-seg');
                curr_new_filename = [data_parent_dir, curr_filename(loc_start_subj_filename-1:end)];
                SPM.xY.VY(f).fname = curr_new_filename;
                SPM.xY.VY(f).private.dat.fname = curr_new_filename;
            end
            % Change SPM directory to correct path
            curr_filename = SPM.swd;
            loc_start_subj_filename = strfind(curr_filename, 'sub-seg');
            curr_new_filename = [data_parent_dir, curr_filename(loc_start_subj_filename-1:end)];
            SPM.swd = curr_new_filename;
            
            % Use the SPM_write_residuals function to save the GLM residuals
            residuals_filenames = spm_write_residuals(SPM,0);
            num_timepoints = length(residuals_filenames);
            
            % Getting this subject's ROIs and finding all of their coordinates
            all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC'};
            all_ROIs_combined = [];
            for i = 1:length(all_ROI_names)
                all_ROIs_combined = cat(4, all_ROIs_combined, spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii']))));
            end
            all_ROIs_combined = nanmean(all_ROIs_combined, 4);
            all_ROIs_coords = find(all_ROIs_combined(:) > 0);
            num_ROI_voxels = length(all_ROIs_coords);
            [all_ROI_coordx,all_ROI_coordy,all_ROI_coordz] = ind2sub(size(all_ROIs_combined), all_ROIs_coords);
            
            % Read residuals only from relevant voxels (using individual ROIs combination)
            all_residuals_ROI = nan(num_timepoints, num_ROI_voxels);
            for t = 1:num_timepoints
                all_residuals_ROI(t,:) = spm_sample_vol(spm_vol(fullfile(curr_analysis_dir, residuals_filenames(t).fname)), all_ROI_coordx, all_ROI_coordy, all_ROI_coordz, 0);
            end
            % Read the original betas
            beta_files = dir(fullfile(curr_analysis_dir,'beta*.nii'));
            betas_original = nan(num_conditions, num_ROI_voxels);
            for b = 1:num_conditions
                betas_original(b,:) = spm_sample_vol(spm_vol(fullfile(curr_analysis_dir, beta_files(b).name)), all_ROI_coordx, all_ROI_coordy, all_ROI_coordz, 0);
            end
            
            % Removing voxels with NaN beta values
            locs_notnan = find(~isnan(sum(betas_original)));
            all_residuals_ROI = all_residuals_ROI(:, locs_notnan);
            betas_original = betas_original(:, locs_notnan);
            
            % Compute the covariance matrix and shrink it
            % Get the covariance matrix of the error residuals, with regularization of the matrix through optimal shrikage(Ledoit & Wolf (2004));
            % using the Decoding Toolbox (TDT) function (Hebart, Gorgen, Haynes 2015), since the RSA Toolbox function has a bug (seehttps://github.com/rsagroup/rsatoolbox_matlab/issues/35)
            regularized_covmat = covshrink_lw2(all_residuals_ROI);
            
            % Normalize the betas by the residuals (error) covariance matrix
            [E,D] = eig(regularized_covmat);
            normalized_betas = betas_original*E*diag(1./real(sqrt(diag(D))))*E';
            % Save to a new image
            for b = 1:num_conditions
                new_image = nan(size(all_ROIs_combined));
                new_image(all_ROIs_coords(locs_notnan)) = normalized_betas(b,:);
                numstr_b = num2str(b); if b<10, numstr_b = ['0', numstr_b]; end
                save_mat_to_nifti(betas_file1, new_image, fullfile(curr_analysis_dir, ['normalized_beta', numstr_b, '.nii']));
                clear new_image;
            end
            
            % Delete the residual files from the analysis folder
            for t = 1:num_timepoints
                curr_redisual_filename = fullfile(curr_analysis_dir, residuals_filenames(t).fname);
                delete(curr_redisual_filename);
            end
            % Clear variables
            clear SPM; clear all_ROIs_combined; clear all_residuals_ROI; clear regularized_covmat; clear betas_original; clear normalized_betas; clear E; clear D;
        end
    end
end


%% JRD single-trial betas - get residuals and do multivariate noise normalization
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    % JRD analysis dirs
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    num_runs_jrd = length(dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat')));
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');
    % Functional localizer analysis dirs
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs(subj).name);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);

    % Defining general SPM GLM parameters
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
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;      % Saving the residuals for the multivariate noise normalization
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    jrd_data_files = dir(fullfile(curr_subj_session_2_dir, 'sm_*jrd*preproc_bold.nii'));
    % Getting this subject's ROIs and finding all of their coordinates
    all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC'};
    all_ROIs_combined = [];
    for roi = 1:length(all_ROI_names)
        all_ROIs_combined = cat(4, all_ROIs_combined, spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{roi}, '_individual.nii']))));
    end
    all_ROIs_combined = nanmean(all_ROIs_combined, 4);
    all_ROIs_coords = find(all_ROIs_combined(:) > 0);
    num_ROI_voxels = length(all_ROIs_coords);
    [all_ROI_coordx,all_ROI_coordy,all_ROI_coordz] = ind2sub(size(all_ROIs_combined), all_ROIs_coords);
    % Go over JRD runs
    for r = 1:length(curr_subj_jrd_analysis_dirs)
        curr_subj_jrd_singletrials_analysis_dir = fullfile(curr_subj_jrd_analysis_dirs{r}, 'single_trials');
        curr_subj_jrd_singletrials_temp_dir = fullfile(curr_subj_jrd_singletrials_analysis_dir, 'temp');
        % Reading the design file for this run
        design_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, '*Segmentation*.mat'));
        if ~isempty(design_file)
            design_file = load(fullfile(design_file(1).folder, design_file(1).name));
            [all_onsets, indices] = sort(cell2mat(design_file.onsets));                                     % Sorting by the order of appearance
            all_durations = cell2mat(design_file.durations); all_durations = all_durations(indices);        % Sorting by the order of the onsets sort
            % Defining general SPM batch parameters
            matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_subj_jrd_singletrials_temp_dir};    % The output directory
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(curr_subj_jrd_singletrials_temp_dir, 'SPM.mat')};     % The resulting design file from the first stage, to be estimated
            % Entering the functional data filenames
            n = nifti(fullfile(jrd_data_files(r).folder, jrd_data_files(r).name)); num_volumes = size(n.dat, 4); clear n;
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cell(num_volumes, 1);
            for v = 1:num_volumes
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans{v} = [fullfile(jrd_data_files(r).folder, jrd_data_files(r).name), ',', num2str(v)];
            end
            % Running the estimation for each trial separately
            for i = 1:length(all_onsets)
                if exist(fullfile(curr_subj_jrd_singletrials_temp_dir, 'SPM.mat')), delete(fullfile(curr_subj_jrd_singletrials_temp_dir, 'SPM.mat')); end
                % Defining all of the regressors - one for the specific trial, the other for all the rest
                onsets = {}; durations = {};
                onsets{1} = all_onsets(i);
                onsets{2} = all_onsets; onsets{2}(i) = [];    % A regressor for all the other trials
                durations{1} = all_durations(i);     % A regressor only with this specific trial
                durations{2} = all_durations; durations{2}(i) = [];
                names = {'current_trial','all_other_trials'};
                % Saving the design file
                temp_design_file_name = fullfile(curr_subj_jrd_singletrials_temp_dir, 'temp_design.mat');
                save(temp_design_file_name, 'names','onsets','durations');
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {temp_design_file_name};      % Defining the design file in the SPM matlab batch
                % Running the design creation and estimation script
                spm_jobman('run', matlabbatch);
                
                % Multivariate noise normalization
                residuals_filenames = dir(fullfile(curr_subj_jrd_singletrials_temp_dir, 'res_*.nii'));
                num_timepoints = length(residuals_filenames);
                % Read residuals only from relevant voxels (using individual ROIs combination)
                all_residuals_ROI = nan(num_timepoints, num_ROI_voxels);
                for t = 1:num_timepoints
                    all_residuals_ROI(t,:) = spm_sample_vol(spm_vol(fullfile(curr_subj_jrd_singletrials_temp_dir, residuals_filenames(t).name)), all_ROI_coordx, all_ROI_coordy, all_ROI_coordz, 0);
                end
                % Read the original betas
                beta_file = fullfile(curr_subj_jrd_singletrials_temp_dir,'beta_0001.nii');
                betas_original = spm_sample_vol(spm_vol(beta_file), all_ROI_coordx, all_ROI_coordy, all_ROI_coordz, 0);
                % Removing voxels with NaN beta values
                locs_notnan = find(~isnan(betas_original));
                betas_original = betas_original(locs_notnan);
                all_residuals_ROI = all_residuals_ROI(:,locs_notnan);
                % Compute the covariance matrix and shrink it
                % Get the covariance matrix of the error residuals, with regularization of the matrix through optimal shrikage(Ledoit & Wolf (2004));
                % using the Decoding Toolbox (TDT) function (Hebart, Gorgen, Haynes 2015), since the RSA Toolbox function has a bug (seehttps://github.com/rsagroup/rsatoolbox_matlab/issues/35)
                regularized_covmat = covshrink_lw2(all_residuals_ROI);
                % Normalize the betas by the residuals (error) covariance matrix
                [E,D] = eig(regularized_covmat);
                normalized_betas = betas_original'*E*diag(1./real(sqrt(diag(D))))*E';
                % Save to a new image
                new_image = nan(size(all_ROIs_combined));
                new_image(all_ROIs_coords(locs_notnan)) = normalized_betas;
                save_mat_to_nifti(betas_file1, new_image, fullfile(curr_subj_jrd_singletrials_temp_dir, ['normalized_beta_0001.nii']));
                % Clear variables
                clear new_image; clear all_residuals_ROI; clear regularized_covmat; clear betas_original; clear normalized_betas; clear E; clear D;
                % Copying the resulting normalized beta1 (the current trial estimate) to the parent directory
                copyfile(fullfile(curr_subj_jrd_singletrials_temp_dir, 'normalized_beta_0001.nii'), fullfile(curr_subj_jrd_singletrials_analysis_dir, ['norm_singletrial_beta_time_',num2str(onsets{1},'%03.0f'),'.nii']));
            end
            % Delete the residual files from the analysis folder
            for t = 1:num_timepoints
                curr_redisual_filename = fullfile(curr_subj_jrd_singletrials_temp_dir, residuals_filenames(t).name);
                delete(curr_redisual_filename);
            end
        end        
    end
end

