
betas_file1 = 'C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\1 - Segmentation\Segmentation_data\sub-seg01\Analysis\Analysis_objectviewing\Run1\beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = 'C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\1 - Segmentation\Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

all_ROI_names = {'bilat_RSC', 'bilat_PPA', 'bilat_OPA'};
num_ROIs = length(all_ROI_names);


num_overlapping_vox_loc_jrd = [];
num_overlapping_vox_loc_objview = [];
num_overlapping_vox_jrd_objview = [];
num_vox_loc = [];
num_vox_jrd = [];
num_vox_objview = [];

for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    % Object viewing analysis dirs
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    
    %% Loading the subject specific ROIs
    for i = 1:length(all_ROI_names)
        curr_ROI_localizer = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
        curr_ROI_activation_jrd = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_activation_jrd_individual.nii'])));
        curr_ROI_activation_objview = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_activation_objview_individual.nii'])));
        num_overlapping_vox_loc_jrd(subj,i) = length(find(curr_ROI_localizer(:)>0 & curr_ROI_activation_jrd(:)>0));
        num_overlapping_vox_loc_objview(subj,i) = length(find(curr_ROI_localizer(:)>0 & curr_ROI_activation_objview(:)>0));
        num_overlapping_vox_jrd_objview(subj,i) = length(find(curr_ROI_activation_jrd(:)>0 & curr_ROI_activation_objview(:)>0));
        num_vox_loc(subj,i) = length(find(curr_ROI_localizer(:)>0));
        num_vox_jrd(subj,i) = length(find(curr_ROI_activation_jrd(:)>0));
        num_vox_objview(subj,i) = length(find(curr_ROI_activation_objview(:)>0));
        
    end
end
