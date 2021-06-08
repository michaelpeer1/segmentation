% Defining general variables
betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC', 'lHCpost','rHCpost','lHCant','rHCant'};
num_ROIs = length(all_ROI_names);



%% Creating the different similarity matrices
% Define locations of objects and create the distance matrix
% (in the paper's figures, the same locations that are here 1-16 are 5,6,...,15,16,1,2,3,4)
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
mat_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
mat_distances = mat_distances / max(mat_distances(:));    % Normalizing to range 0-1
% Converting all matrices from dissimilarity to similarity matrices
mat_distances = 1 - mat_distances;


% Overlay - with distance
locations_new_overlay = locations_all_objects;
locations_new_overlay(5:12,2) = locations_new_overlay(5:12,2) + 75;
mat_overlay = pdist2(locations_new_overlay, locations_new_overlay);         % Distances between all objects
mat_overlay = mat_overlay / max(mat_overlay(:));    % Normalizing to range 0-1
mat_overlay = 1 - mat_overlay;

% Lower triangle elements only
mat_nondiag_distance = mat_distances(find(tril(ones(num_conditions),-1)));
mat_nondiag_overlay = mat_overlay(find(tril(ones(num_conditions),-1)));


for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    % Object viewing analysis dirs
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objectviewing_analysis_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2')};
    num_runs_objectviewing = length(curr_subj_objectviewing_analysis_dirs);
    % Functional localizer analysis dirs
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    
    
    %% Loading the subject specific ROIs
    all_ROIs_individual = {};
    for i = 1:length(all_ROI_names)
        all_ROIs_individual{i} = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
    end
    
    
    for roi = 1:num_ROIs        % Going over ROIs
        curr_ROI = find(all_ROIs_individual{roi}(:)==1);
        if ~isempty(curr_ROI)
            % Getting the ROI xyz coordinates
            [curr_ROI_coordx,curr_ROI_coordy,curr_ROI_coordz] = ind2sub(size(all_ROIs_individual{roi}),curr_ROI);
            
            %% Object viewing data
            all_betas_ROI_objview = nan(length(curr_ROI), num_conditions, num_runs_objectviewing);
            all_betas_ROI_objview_meanremoved = nan(length(curr_ROI), num_conditions, num_runs_objectviewing);
            for run = 1:num_runs_objectviewing
                % Get the already noise-normalized betas (obtained with the temp_noise_normalize_betas script)
                norm_beta_files_current = dir(fullfile(curr_subj_objectviewing_analysis_dirs{run}, 'normalized_beta*.nii'));
                for b = 1:num_conditions
                    curr_norm_beta_image = spm_sample_vol(spm_vol(fullfile(curr_subj_objectviewing_analysis_dirs{run}, norm_beta_files_current(b).name)), curr_ROI_coordx, curr_ROI_coordy, curr_ROI_coordz, 0);
                    all_betas_ROI_objview(:,b,run) = curr_norm_beta_image;
                end
                
                % Removing the cocktail mean from each day (mean activity at each voxel across patterns)
                mean_run_pattern = nanmean(all_betas_ROI_objview(:,:,run), 2);
                for b = 1:num_conditions
                    all_betas_ROI_objview_meanremoved(:,b,run) = all_betas_ROI_objview(:,b,run) - mean_run_pattern;
                end
            end
            
            % Getting the average pattern across runs, and calculating the correlation between these average patterns
            all_betas_day1 = nanmean(all_betas_ROI_objview_meanremoved(:,:,1:2), 3);
            corr_ROI_objview_day1 = fisherz(corr(all_betas_day1, all_betas_day1, 'rows', 'complete'));            
            % Non-diagonal elements
            corr_ROI_objview_day1_nondiag = corr_ROI_objview_day1(find(tril(ones(num_conditions),-1)));       % Taking only the lower triangle elements as a vector
            all_corrs_objview_day1(:,subj,roi) = corr_ROI_objview_day1_nondiag;
            % Correlation to different matrices
            corr_objview_d2minusd1_distances(subj, roi) = corr(corr_ROI_objview_day1_nondiag, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_overlay(subj, roi) = corr(corr_ROI_objview_day1_nondiag, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
            
            
            % Same thing for the patterns without the mean removed
            all_betas_day1 = nanmean(all_betas_ROI_objview(:,:,1:2), 3);
            corr_ROI_objview_day1 = fisherz(corr(all_betas_day1, all_betas_day1, 'rows', 'complete'));            
            % Non-diagonal elements
            corr_ROI_objview_day1_nondiag = corr_ROI_objview_day1(find(tril(ones(num_conditions),-1)));       % Taking only the lower triangle elements as a vector
            all_corrs_objview_day1_nomeanremoval(:,subj,roi) = corr_ROI_objview_day1_nondiag;
            % Correlation to different matrices
            corr_objview_d2minusd1_distances_nomeanremoval(subj, roi) = corr(corr_ROI_objview_day1_nondiag, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_overlay_nomeanremoval(subj, roi) = corr(corr_ROI_objview_day1_nondiag, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');

        end
    end
end
