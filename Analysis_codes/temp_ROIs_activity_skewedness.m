% Test activation strength across ROIs
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data';
names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'framewise_displacement', 'csf', 'white_matter'};      % Confounds to use
num_confounds = length(names);
% Reading the subjects file, with the path of all logfiles for each subject
[~, ~, subjects_file_data]= xlsread(fullfile(data_parent_dir, 'fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx')); 
subjects_file_data = subjects_file_data(4:end, :);
for i = 1:size(subjects_file_data,1), if isnan(subjects_file_data{i,1}), subjects_file_data = subjects_file_data(1:i-1,:); break; end,end       % Getting rid of rows without subjects data
subj_data_dirs = subjects_file_data(:,29);
num_betas = 16;

betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
lRSC = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/lRSC.nii', betas_file1, 0);
rRSC = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/rRSC.nii', betas_file1, 0);
lPPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/lPPA.nii', betas_file1, 0);
rPPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/rPPA.nii', betas_file1, 0);
lOPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/lTOS.nii', betas_file1, 0);
rOPA = reslice_data('/Users/mpeer/Dropbox/Michael_templates/func_localizer_parcels_Julian_2012/rTOS.nii', betas_file1, 0);
bilat_RSC = lRSC + rRSC; bilat_PPA = lPPA + rPPA; bilat_OPA = lOPA + rOPA;

all_ROI_names = {'bilat_RSC', 'bilat_PPA', 'bilat_OPA', 'bilat_HC', 'bilat_EC'};
temp_ROIs_jrd = cell(1,length(all_ROI_names));
temp_ROIs_objview = cell(1,length(all_ROI_names));
for subj = 1:length(subj_data_dirs)
    disp(subj)
    
    % Defining current subject directories
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs{subj}), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    %     curr_subj_jrd_univariate_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_univariate');
    %     curr_subj_objectviewing_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing_adaptation');
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
    left_EC = curr_parcellation == 1006; right_EC = curr_parcellation == 2006; bilat_EC = left_EC + right_EC;
    all_ROIs = {bilat_RSC, bilat_PPA, bilat_OPA, bilat_HC, bilat_EC};

    % Reading the voxels activity in the JRD task
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
        for i=1:length(all_ROIs)
            curr_locations = find(all_ROIs{i}(:) > 0);
            temp_ROIs_jrd{i}{subj} = jrd_contrast(curr_locations);
        end
    else
        for i=1:length(all_ROIs)
            temp_ROIs_jrd{i}{subj} = nan;
        end
    end
    % Reading the voxels activity in object-viewing day2
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
        for i=1:length(all_ROIs)
            curr_locations = find(all_ROIs{i}(:) > 0);
            temp_ROIs_objview{i}{subj} = objectviewing_contrast(curr_locations);
        end
    else
        for i=1:length(all_ROIs)
            temp_ROIs_objview{i}{subj} = nan;
        end
    end
end


% Plot all subjects' ROI activity
% edges = [-5:0.2:5];
% a=[];b=[];aa=[];bb=[];
% for i=1:length(temp_ROIs_jrd{4})
%     a(i,:) = histcounts(temp_ROIs_jrd{4}{i},edges);
%     b(i,:) = histcounts(temp_ROIs_objview{4}{i},edges);
%     aa(i,:) = histcounts(temp_ROIs_jrd{4}{i},20);
%     bb(i,:) = histcounts(temp_ROIs_objview{4}{i},20);
% end
roi = 3;
figure;subplot(6,4,1);
for i=1:length(temp_ROIs_jrd{roi})
    subplot(6,4,i);
    hist(temp_ROIs_jrd{roi}{i},20);
end
figure;subplot(6,4,1);
for i=1:length(temp_ROIs_objview{roi})
    subplot(6,4,i);
    hist(temp_ROIs_objview{roi}{i},20);
end


% Calculate skewness of the distribution, based on the Fisher-Pearson coefficient of skewness
all_skewness_jrd = []; all_skewness_objview = [];
for roi = 1:length(temp_ROIs_jrd)
    for i = 1:length(temp_ROIs_jrd{roi})
        current_distribution = temp_ROIs_jrd{roi}{i};
%         all_skewness_jrd(roi,i) = nansum((current_distribution - nanmean(current_distribution)).^3) / (length(current_distribution) * nanstd(current_distribution)^3);
        all_nonparam_skewness_jrd(roi,i) = (nanmean(current_distribution) - nanmedian(current_distribution)) / nanstd(current_distribution);
        current_distribution = temp_ROIs_objview{roi}{i};
%         all_skewness_objview(roi,i) = nansum((current_distribution - nanmean(current_distribution)).^3) / (length(current_distribution) * nanstd(current_distribution)^3);
        all_nonparam_skewness_objview(roi,i) = (nanmean(current_distribution) - nanmedian(current_distribution)) / nanstd(current_distribution);        
    end
end
