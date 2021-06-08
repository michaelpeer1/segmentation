% Calculating the tSNR across scans (tSNR equals the mean of the fMRI signal, divided
% by its standard deviation - Murphy et al., NeuroImage 2007

betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC'};
num_ROIs = length(all_ROI_names);


all_subjs_ROIs_tSNR = nan(length(subj_data_dirs), num_ROIs, 7);
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    curr_subj_parent_dir = fullfile(data_parent_dir, subj_data_dirs(subj).name);
    curr_subj_session_1_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day1/func']);
    curr_subj_session_2_dir = fullfile(curr_subj_parent_dir, ['fmriprep/sub-seg', curr_subj_parent_dir(end-1:end), '/ses-day2/func']);
    curr_subj_objview_func_files_day1 = dir(fullfile(curr_subj_session_1_dir, 'sm*object_viewing*preproc_bold.nii'));
    curr_subj_objview_func_files_day2 = dir(fullfile(curr_subj_session_2_dir, 'sm*object_viewing*preproc_bold.nii'));
    curr_subj_objview_func_files_day1 = fullfile(curr_subj_session_1_dir, {curr_subj_objview_func_files_day1.name});
    curr_subj_objview_func_files_day2 = fullfile(curr_subj_session_2_dir, {curr_subj_objview_func_files_day2.name});
    curr_subj_jrd_func_files = dir(fullfile(curr_subj_session_2_dir, 'sm*jrd*preproc_bold.nii'));
    curr_subj_jrd_func_files = fullfile(curr_subj_session_2_dir, {curr_subj_jrd_func_files.name});
    curr_subj_all_func_files = [curr_subj_objview_func_files_day1, curr_subj_objview_func_files_day2, curr_subj_jrd_func_files];
    
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    
    %% Loading the subject specific ROIs
    all_ROIs_individual = {};
    for i = 1:length(all_ROI_names)
        all_ROIs_individual{i} = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
    end
    
    %% Compute tSNR
    for f = 1:length(curr_subj_all_func_files)
        % Read the functional data of each run
        temp_file = spm_read_vols(spm_vol(curr_subj_all_func_files{f}));
        sz = size(temp_file);
        % reshaping to voxels X time
        temp_file = reshape(temp_file, [sz(1)*sz(2)*sz(3), sz(4)]);
        
        for roi = 1:num_ROIs        % Going over ROIs
            curr_ROI = find(all_ROIs_individual{roi}(:)==1);
            
            % Compute average signal in ROI across time
            curr_ROI_signal = nanmean(temp_file(curr_ROI, :));
            
            % Compute tSNR for that run in each ROI
            curr_mean = mean(curr_ROI_signal);
            curr_std = std(curr_ROI_signal);
            
            all_subjs_ROIs_tSNR(subj, roi, f) = curr_mean / curr_std;
        end
        clear temp_file;
    end
end


% Average tSNR across subjects

t = nanmean(all_subjs_ROIs_tSNR,3);
tt = [nanmean(t(:,[1,5]),2), nanmean(t(:,[2,6]),2), nanmean(t(:,[3,7]),2), t(:,4),t(:,8)];      % RSC, PPA, OPA, HC, ERC
tmean = nanmean(tt);
