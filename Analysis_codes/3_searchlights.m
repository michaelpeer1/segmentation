% Searchlights
betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;


%% Defining dissimilarity matrices
% Define locations of objects and create the distance matrix
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
mat_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
mat_distances = mat_distances / max(mat_distances(:));    % Normalizing to range 0-1
% Create the segment identity similarity matrix
mat_segments = ones(num_conditions);
mat_segments(1:4,1:4) = 0; mat_segments(13:16,13:16) = 0; mat_segments(5:12,5:12) = 0; mat_segments(13:16,1:4) = 0; mat_segments(1:4,13:16) = 0;
% Create the orthogonal segment identity similarity matrix
mat_orth_segments = ones(num_conditions);
mat_orth_segments(1:8,1:8) = 0; mat_orth_segments(9:16,9:16) = 0;
% Create the quadrant identity similarity matrix
mat_quadrants = ones(num_conditions);
mat_quadrants(1:4,1:4) = 0; mat_quadrants(5:8,5:8) = 0; mat_quadrants(9:12,9:12) = 0; mat_quadrants(13:16,13:16) = 0;
% Converting all matrices from dissimilarity to similarity matrices
mat_distances = 1 - mat_distances; mat_segments = 1 - mat_segments; mat_quadrants = 1 - mat_quadrants; mat_orth_segments = 1 - mat_orth_segments;

% Schema preserving - with distance
locations_new_schema = locations_all_objects;
for i = 1:4
    locations_new_schema(i, :) = locations_new_schema(i + num_conditions/2, :);
end
for i = 13:16
    locations_new_schema(i, :) = locations_new_schema(i - num_conditions/2, :);
end
mat_schema_preserving_w_dist = pdist2(locations_new_schema, locations_new_schema);         % Distances between all objects
mat_schema_preserving_w_dist = mat_schema_preserving_w_dist / max(mat_schema_preserving_w_dist(:));    % Normalizing to range 0-1
mat_schema_preserving_w_dist = 1 - mat_schema_preserving_w_dist;
% Flipping - with distance
locations_new_flipping = locations_all_objects;
locations_new_flipping(:,2) = abs(locations_new_flipping(:,2));
mat_flipping_w_dist = pdist2(locations_new_flipping, locations_new_flipping);         % Distances between all objects
mat_flipping_w_dist = mat_flipping_w_dist / max(mat_flipping_w_dist(:));    % Normalizing to range 0-1
mat_flipping_w_dist = 1 - mat_flipping_w_dist;
% Overlay - with distance
locations_new_overlay = locations_all_objects;
locations_new_overlay(5:12,2) = locations_new_overlay(5:12,2) + 75;
mat_overlay_w_dist = pdist2(locations_new_overlay, locations_new_overlay);         % Distances between all objects
mat_overlay_w_dist = mat_overlay_w_dist / max(mat_overlay_w_dist(:));    % Normalizing to range 0-1
mat_overlay_w_dist = 1 - mat_overlay_w_dist;
% Organization on Y axis (along river direction)
mat_y_axis = pdist2(locations_all_objects(:,1), locations_all_objects(:,1));
mat_y_axis = mat_y_axis / max(mat_y_axis(:));    % Normalizing to range 0-1
mat_y_axis = 1 - mat_y_axis;

% A matrix coding if it is an adjacent quadrant within and between segment
mat_within_seg_adj_quad = zeros(num_conditions);
mat_within_seg_adj_quad(1:4,13:16) = 1; mat_within_seg_adj_quad(5:8,9:12) = 1; mat_within_seg_adj_quad(9:12,5:8) = 1; mat_within_seg_adj_quad(13:16,1:4) = 1;
mat_between_seg_adj_quad = zeros(num_conditions);
mat_between_seg_adj_quad(1:4,5:8) = 1; mat_between_seg_adj_quad(5:8,1:4) = 1; mat_between_seg_adj_quad(9:12,13:16) = 1; mat_between_seg_adj_quad(13:16,9:12) = 1;
% Distance within segment, adjacent quadrants only
mat_dist_within_segment_adj_quadrants = mat_distances;
mat_dist_within_segment_adj_quadrants(mat_segments==0) = nan;
mat_dist_within_segment_adj_quadrants(mat_quadrants==1) = nan;
% Distance between segments, adjacent quadrants only
mat_diagonal_quadrants = zeros(16); mat_diagonal_quadrants(1:4,9:12) = 1; mat_diagonal_quadrants(9:12, 1:4) = 1; mat_diagonal_quadrants(5:8,13:16) = 1; mat_diagonal_quadrants(13:16,5:8) = 1;
mat_dist_between_segments_adj_quadrants = mat_distances;
mat_dist_between_segments_adj_quadrants(mat_segments==1) = nan;
mat_dist_between_segments_adj_quadrants(mat_diagonal_quadrants==1) = nan;

% Distance within segment only
mat_dist_within_segment = mat_distances;
mat_dist_within_segment(mat_segments==0) = nan;
% Distance between segments only
mat_dist_between_segments = mat_distances;
mat_dist_between_segments(mat_segments==1) = nan;

% A matrix coding if it is an adjacent quadrant within and between segment
mat_within_seg_adj_quad = zeros(num_conditions);
mat_within_seg_adj_quad(1:4,13:16) = 1; mat_within_seg_adj_quad(5:8,9:12) = 1; mat_within_seg_adj_quad(9:12,5:8) = 1; mat_within_seg_adj_quad(13:16,1:4) = 1;
mat_between_seg_adj_quad = zeros(num_conditions);
mat_between_seg_adj_quad(1:4,5:8) = 1; mat_between_seg_adj_quad(5:8,1:4) = 1; mat_between_seg_adj_quad(9:12,13:16) = 1; mat_between_seg_adj_quad(13:16,9:12) = 1;

% A matrix coding a segment grouping index - pattern similarity within segment larger than between segments (adjacent quadrants only)
mat_segment_grouping = nan(num_conditions);
mat_segment_grouping(mat_within_seg_adj_quad == 1) = 1;
mat_segment_grouping(mat_between_seg_adj_quad == 1) = -1;


% Getting the nondiagonal elements
mat_nondiag_distance = mat_distances(find(tril(ones(num_conditions),-1)));
mat_nondiag_segments = mat_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_orthsegments = mat_orth_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_quadrants = mat_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_flipping = mat_flipping_w_dist(find(tril(ones(num_conditions),-1)));
mat_nondiag_overlay = mat_overlay_w_dist(find(tril(ones(num_conditions),-1)));
mat_nondiag_schema = mat_schema_preserving_w_dist(find(tril(ones(num_conditions),-1)));
mat_nondiag_y_axis = mat_y_axis(find(tril(ones(num_conditions),-1)));
mat_nondiag_within_seg_adj_quad = mat_within_seg_adj_quad(find(tril(ones(num_conditions),-1)));
mat_nondiag_between_seg_adj_quad = mat_between_seg_adj_quad(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_within_segment_adj_quadrants = mat_dist_within_segment_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_between_segments_adj_quadrants = mat_dist_between_segments_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_within_segment = mat_dist_within_segment(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_between_segments = mat_dist_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_segment_grouping = mat_segment_grouping(find(tril(ones(num_conditions),-1)));




%% HC and whole-brain and all-ROIs searchlight in each subject

searchlight_radius = 3;

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
    
    % Getting the ROIs of this subject from freesurfer parcellation
    curr_subject_functional_directory = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'fmriprep');
    curr_subject_functional_directory = fullfile(fullfile(curr_subject_functional_directory, subj_data_dirs(subj).name), 'ses-day1/func/');
    % Hippocampus ROI
    curr_parcellation_file = dir(fullfile(curr_subject_functional_directory, '*aparcaseg*.nii'));
    curr_parcellation = reslice_data(fullfile(curr_subject_functional_directory, curr_parcellation_file(1).name), betas_file1, 0);
    hippocampus_ROI = (curr_parcellation == 17 | curr_parcellation == 53);
    % Whole-brain ROI
    curr_parcellation_file = dir(fullfile(curr_subject_functional_directory, '*-aseg_dseg.nii'));
    curr_parcellation = reslice_data(fullfile(curr_subject_functional_directory, curr_parcellation_file(1).name), betas_file1, 0);
    wholebrain_ROI = (curr_parcellation == 3 | curr_parcellation == 17 | curr_parcellation == 18 | curr_parcellation == 42 | curr_parcellation == 53 | curr_parcellation == 54);       % Cortex, hippocampus, amygdala
    % Getting the Julian et al. group ROIs and combining to a large ROI of all + hippocampus
    lRSC = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/lRSC.nii', betas_file1, 0);
    rRSC = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/rRSC.nii', betas_file1, 0);
    lPPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/lPPA.nii', betas_file1, 0);
    rPPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/rPPA.nii', betas_file1, 0);
    lOPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/lTOS.nii', betas_file1, 0);
    rOPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/rTOS.nii', betas_file1, 0);
    combinedROIs = lRSC+rRSC+lPPA+rPPA+lOPA+rOPA+hippocampus_ROI;
    combinedROIs = combinedROIs > 0;
    
    all_searchlight_ROIs = {wholebrain_ROI};
    all_searchlight_ROI_names = {'wholebrain'};
    size_overall = size(curr_parcellation,1)*size(curr_parcellation,2)*size(curr_parcellation,3);
    
    
    for roi = 1:length(all_searchlight_ROIs)
        curr_ROI = all_searchlight_ROIs{roi};
        curr_ROI_voxels = find(curr_ROI);   % Getting all of the indices of voxels in the ROI
        
        all_tvalues_ROI_objectviewing = nan(size_overall, num_conditions, length(curr_subj_objectviewing_combined_dirs));
        for r = 1:length(curr_subj_objectviewing_combined_dirs)
            t_values_files_current = dir(fullfile(curr_subj_objectviewing_combined_dirs{r}, 'spmT*.nii'));
            for b = 1:num_conditions
                curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_objectviewing_combined_dirs{r}, t_values_files_current(b).name)));
                all_tvalues_ROI_objectviewing(:,b,r) = curr_tvalues_image(:);
            end
        end
        % Removing the cocktail mean from each day (mean activity at each voxel across patterns)
        for r = 1:length(curr_subj_objectviewing_combined_dirs)
            mean_day_pattern = nanmean(all_tvalues_ROI_objectviewing(:,:,r), 2);
            for b = 1:num_conditions
                all_tvalues_ROI_objectviewing(:,b,r) = all_tvalues_ROI_objectviewing(:,b,r) - mean_day_pattern;
            end
        end
        if num_runs_jrd > 0
            all_tvalues_ROI_jrd = nan(size_overall, num_conditions);
            t_values_files_current = dir(fullfile(curr_subj_jrd_combined_dir, 'spmT*.nii'));
            for b = 1:num_conditions
                curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_jrd_combined_dir, t_values_files_current(b).name)));
                all_tvalues_ROI_jrd(:,b) = curr_tvalues_image(:);
            end
            % Removing the cocktail mean (mean activity at each voxel across patterns)
            mean_jrd_pattern = nanmean(all_tvalues_ROI_jrd, 2);
            for b = 1:num_conditions
                all_tvalues_ROI_jrd(:,b) = all_tvalues_ROI_jrd(:,b) - mean_jrd_pattern;
            end
        end
        
        % Running the searchlight
        ROI_searchlight_results_objectviewing_distances = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_segments = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_quadrants = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_schema = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_flipping = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_overlay = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_yaxis = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_distances = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_segments = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_quadrants = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_schema = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_flipping = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_overlay = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_yaxis = nan(size(curr_ROI));
        ROI_searchlight_partialvar_objview_day2minusday1_distances = nan(size(curr_ROI));
        ROI_searchlight_partialvar_objview_day2minusday1_segments = nan(size(curr_ROI));
        ROI_searchlight_partialvar_objview_day2minusday1_quadrants = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_grouping = nan(size(curr_ROI));
        ROI_searchlight_results_objectviewing_remapping = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_grouping = nan(size(curr_ROI));
        ROI_searchlight_results_jrd_remapping = nan(size(curr_ROI));
        
        for i = 1:length(curr_ROI_voxels)
            % Defining a sphere around the voxel according to the searchlight radius
            [curr_vox_coord_x, curr_vox_coord_y, curr_vox_coord_z] = ind2sub(size(curr_ROI), curr_ROI_voxels(i));        % Getting the 3D coordinate of the current voxel
            [x y z] = ind2sub(size(curr_ROI),1:numel(curr_ROI));      % Getting the 3D coordinates of all voxels
            M = [x', y', z'] - repmat([curr_vox_coord_x, curr_vox_coord_y, curr_vox_coord_z], numel(curr_ROI), 1);     % The distance along each dimension of each voxel from the current one
            D = (sqrt(sum(M.^2,2)) <= searchlight_radius);        % Identifying only the locations whose sum of squared distances from the voxels are equal to or smaller than the searchlight radius
            D = reshape(D, size(curr_ROI));        % Reshaping the sphere locations to the original matrix size
            
            % Masking with the ROI to get only voxels that are in it and in the searchlight radius
            curr_searchlight_voxels = find(D ~= 0 & curr_ROI ~= 0);
            
            % Calculating the result in the current searchlight sphere - objectviewing day2 minus day1 correlation, correlation to distances and segments matrices
            patterns_day1 = nanmean(all_tvalues_ROI_objectviewing(curr_searchlight_voxels,:,1), 3);
            curr_corr_ROI_objectviewing_day1 = fisherz(corr(patterns_day1, patterns_day1, 'rows', 'complete'));                    % Correlation between day 1 patterns, averaged across runs
            curr_corr_ROI_objectviewing_day1_nondiag = curr_corr_ROI_objectviewing_day1(find(tril(ones(num_conditions),-1)));       % Taking only the lower triangle elements as a vector
            patterns_day2 = nanmean(all_tvalues_ROI_objectviewing(curr_searchlight_voxels,:,2), 3);
            curr_corr_ROI_objectviewing_day2 = fisherz(corr(patterns_day2, patterns_day2, 'rows', 'complete'));                    % Correlation between day 2 patterns, averaged across runs
            curr_corr_ROI_objectviewing_day2_nondiag = curr_corr_ROI_objectviewing_day2(find(tril(ones(num_conditions),-1)));       % Taking only the lower triangle elements as a vector
            corr_ROI_objview_d2minusd1_nondiag = curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag;
            ROI_searchlight_results_objectviewing_distances(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            ROI_searchlight_results_objectviewing_segments(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
            ROI_searchlight_results_objectviewing_quadrants(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            ROI_searchlight_results_objectviewing_schema(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_schema, 'Type', 'Spearman', 'rows', 'complete');
            ROI_searchlight_results_objectviewing_flipping(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_flipping, 'Type', 'Spearman', 'rows', 'complete');
            ROI_searchlight_results_objectviewing_overlay(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
            ROI_searchlight_results_objectviewing_yaxis(curr_ROI_voxels(i)) = corr((curr_corr_ROI_objectviewing_day2_nondiag - curr_corr_ROI_objectviewing_day1_nondiag), mat_nondiag_y_axis, 'Type', 'Spearman', 'rows', 'complete');
            % Segment grouping
            ROI_searchlight_results_objectviewing_grouping(curr_ROI_voxels(i)) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
            % Segment remapping
            ROI_searchlight_results_objectviewing_remapping(curr_ROI_voxels(i)) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete') - corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
            
            
            if num_runs_jrd > 0
                % Getting patterns and correlating to different matrices
                corr_ROI_jrd_nondiag = fisherz(corr(all_tvalues_ROI_jrd(curr_searchlight_voxels,:), all_tvalues_ROI_jrd(curr_searchlight_voxels,:), 'rows', 'complete'));
                corr_ROI_jrd_nondiag = corr_ROI_jrd_nondiag(find(tril(ones(num_conditions),-1)));   % Taking only the lower triangle elements as a vector
                ROI_searchlight_results_jrd_distances(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
                ROI_searchlight_results_jrd_segments(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
                ROI_searchlight_results_jrd_quadrants(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_quadrants, 'Type', 'Spearman', 'rows', 'complete');
                ROI_searchlight_results_jrd_schema(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_schema, 'Type', 'Spearman', 'rows', 'complete');
                ROI_searchlight_results_jrd_flipping(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_flipping, 'Type', 'Spearman', 'rows', 'complete');
                ROI_searchlight_results_jrd_overlay(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
                ROI_searchlight_results_jrd_yaxis(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_y_axis, 'Type', 'Spearman', 'rows', 'complete');
                % Segment grouping
                ROI_searchlight_results_jrd_grouping(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
                % Segment remapping
                ROI_searchlight_results_jrd_remapping(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete') - corr(corr_ROI_jrd_nondiag(:), mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
            end
        end
        
        % Saving the searchlight results
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_distances, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_distances.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_segments, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_segments.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_quadrants, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_quadrants.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_schema, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_schema.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_flipping, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_flipping.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_overlay, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_overlay.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_yaxis, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_yaxis.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_grouping, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_grouping.nii']));
        save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_objectviewing_remapping, fullfile(curr_subj_objectviewing_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_objectviewing_day2minusday1_remapping.nii']));        
        if num_runs_jrd > 0
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_distances, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_distances.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_segments, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_segments.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_quadrants, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_quadrants.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_schema, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_schema.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_flipping, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_flipping.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_overlay, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_overlay.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_yaxis, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_yaxis.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_grouping, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_grouping.nii']));
            save_mat_to_nifti(fullfile(t_values_files_current(1).folder, t_values_files_current(1).name), ROI_searchlight_results_jrd_remapping, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_', all_searchlight_ROI_names{roi} ,'_avgruns_searchlight_jrd_remapping.nii']));
        end
    end
end




%% Smoothing individual subjects' searchlight results

spm('defaults', 'FMRI');
for subj = 1:length(subj_data_dirs)
    disp(subj)
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    
    % Object viewing runs
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    %     curr_searchlight_files_objectviewing = dir(fullfile(curr_subj_objectviewing_analysis_dir, '*searchlight*.nii'));
    curr_searchlight_files_objectviewing = dir(fullfile(curr_subj_objectviewing_analysis_dir, '*avgruns_searchlight*.nii'));
    curr_searchlight_files_objectviewing = fullfile(curr_subj_objectviewing_analysis_dir, {curr_searchlight_files_objectviewing.name});
    % JRD runs
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    %     curr_searchlight_files_jrd = dir(fullfile(curr_subj_jrd_analysis_dir, '*searchlight*.nii'));
    curr_searchlight_files_jrd = dir(fullfile(curr_subj_jrd_analysis_dir, '*avgruns_searchlight*.nii'));
    curr_searchlight_files_jrd = fullfile(curr_subj_jrd_analysis_dir, {curr_searchlight_files_jrd.name});
    % Combining JRD and objectviewing files
    curr_searchlight_files = [curr_searchlight_files_objectviewing, curr_searchlight_files_jrd];
    
    % Removing smoothed files from the files to smooth
    curr_searchlight_files(cellfun(@(x) ~isempty(strfind(x,'sm_')), curr_searchlight_files)) = [];
    
    % Defining the spm batch parameters for smoothing
    matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'sm_';
    matlabbatch{1}.spm.spatial.smooth.data = curr_searchlight_files';     % The actual data
    
    % Running the smoothing
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch;
    
    % Masking with individual subject brain mask, to remove non-brain signals resulting from smoothing
    curr_subject_functional_directory = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'fmriprep');
    curr_subject_functional_directory = fullfile(fullfile(curr_subject_functional_directory, subj_data_dirs(subj).name), 'ses-day1/func/');
    curr_brainmask_file = dir(fullfile(curr_subject_functional_directory, '*brain_mask*.nii'));
    curr_brainmask = reslice_data(fullfile(curr_subject_functional_directory, curr_brainmask_file(1).name), betas_file1, 0);
    % Getting the smoothed files
    curr_searchlight_files_objectviewing = dir(fullfile(curr_subj_objectviewing_analysis_dir, 'sm_*avgruns_searchlight*.nii'));
    curr_searchlight_files_objectviewing = fullfile(curr_subj_objectviewing_analysis_dir, {curr_searchlight_files_objectviewing.name});
    curr_searchlight_files_jrd = dir(fullfile(curr_subj_jrd_analysis_dir, 'sm_*avgruns_searchlight*.nii'));
    curr_searchlight_files_jrd = fullfile(curr_subj_jrd_analysis_dir, {curr_searchlight_files_jrd.name});
    curr_searchlight_files = [curr_searchlight_files_objectviewing, curr_searchlight_files_jrd];
    % Reading each file, masking and saving
    for i = 1:length(curr_searchlight_files)
        curr_file = curr_searchlight_files{i};
        a = spm_read_vols(spm_vol(curr_file));
        a(curr_brainmask==0) = 0;
        save_mat_to_nifti(curr_file, a, curr_file);
    end
    
end




%% TFCE group analysis results using cosmomvpa
group_analysis_dir = fullfile(data_parent_dir, 'Analysis_group');
num_iterations = 10000;     % Number of iterations for the permutation testing

RSA_filenames = {'_avgruns_searchlight_objectviewing_day2minusday1_distances.nii',...
    '_avgruns_searchlight_objectviewing_day2minusday1_overlay.nii',...
    '_avgruns_searchlight_objectviewing_day2minusday1_grouping.nii',...
    '_avgruns_searchlight_objectviewing_day2minusday1_remapping.nii',...
    '_avgruns_searchlight_jrd_distances.nii',...
    '_avgruns_searchlight_jrd_overlay.nii',...
    '_avgruns_searchlight_jrd_grouping.nii',...
    '_avgruns_searchlight_jrd_remapping.nii',...
    };
%     '_avgruns_searchlight_objectviewing_day2minusday1_schema.nii',...
%     '_avgruns_searchlight_objectviewing_day2minusday1_flipping.nii',...
%     '_avgruns_searchlight_objectviewing_day2minusday1_yaxis.nii',...
%     '_avgruns_searchlight_jrd_schema.nii',...
%     '_avgruns_searchlight_jrd_flipping.nii',...
%     '_avgruns_searchlight_jrd_yaxis.nii',...


% Defining ROIs for analysis
% Whole-brain ROI from the AAL3 atlas
wholebrain_ROI = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/AAL3v1_1mm.nii', betas_file1, 0); 
wholebrain_ROI = wholebrain_ROI > 0;
% Combined ROIs (RSC, OPA, PPA from Julian 2012 NeuroImage and hippocampus from the AAL3 atlas)
lRSC = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/lRSC.nii', betas_file1, 0);
rRSC = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/rRSC.nii', betas_file1, 0);
lPPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/lPPA.nii', betas_file1, 0);
rPPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/rPPA.nii', betas_file1, 0);
lOPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/lTOS.nii', betas_file1, 0);
rOPA = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/func_localizer_parcels_Julian_2012/rTOS.nii', betas_file1, 0);
% MTL ROI (hippocampus + entorhinal) - combination of all subjects' individual ROIs
% hippocampus_ROI = reslice_data('/Users/mpeer/Dropbox/Michael_scripts/Templates/AAL3_hipp.nii', betas_file1, 0);
for subj=1:length(subj_data_dirs)
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    curr_HC = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, 'bilat_HC_individual.nii')));
    curr_EC = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, 'bilat_EC_individual.nii')));
    hc = hc+curr_HC; ec = ec+curr_EC;
end
MTL_ROI = hc + ec; MTL_ROI = MTL_ROI > 0;
% Combining the ROIs
% combinedROIs = lRSC+rRSC+lPPA+rPPA+lOPA+rOPA+hippocampus_ROI;
combinedROIs = lRSC+rRSC+lPPA+rPPA+lOPA+rOPA+MTL_ROI;
combinedROIs = combinedROIs > 0;
% Saving the ROIs
save_mat_to_nifti(betas_file1, wholebrain_ROI, fullfile(group_analysis_dir, 'wholebrain_ROI.nii'));
% save_mat_to_nifti(betas_file1, hippocampus_ROI, fullfile(group_analysis_dir, 'hippocampus_ROI.nii'));
save_mat_to_nifti(betas_file1, MTL_ROI, fullfile(group_analysis_dir, 'MTL_ROI.nii'));
save_mat_to_nifti(betas_file1, ec>0, fullfile(group_analysis_dir, 'EC_individuals_combined_ROI.nii'));
save_mat_to_nifti(betas_file1, hc>0, fullfile(group_analysis_dir, 'HC_individuals_combined_ROI.nii'));
save_mat_to_nifti(betas_file1, MTL_ROI, fullfile(group_analysis_dir, 'MTL_ROI.nii'));
save_mat_to_nifti(betas_file1, combinedROIs, fullfile(group_analysis_dir, 'combinedROIs_new.nii'));

% all_searchlight_mask_files = {fullfile(group_analysis_dir, 'wholebrain_ROI.nii'), fullfile(group_analysis_dir, 'hippocampus_ROI.nii'), fullfile(group_analysis_dir, 'combinedROIs.nii')};
% all_searchlight_mask_names = {'_wholebrain_', '_HC_', '_combinedROIs_'};
all_searchlight_mask_files = {fullfile(group_analysis_dir, 'wholebrain_ROI.nii'), fullfile(group_analysis_dir, 'MTL_ROI.nii'), fullfile(group_analysis_dir, 'combinedROIs_new.nii')};
all_searchlight_mask_names = {'_wholebrain_', '_MTL_', '_combinedROIs_new_'};

% Running the group analysis separately for different dissimilarity matrix results, in each ROI
for i = 1:length(RSA_filenames)
    RSA_filename = RSA_filenames{i};
    disp(['Current file: ', RSA_filename])
    
    for roi = 1:length(all_searchlight_mask_files)
        % Reading the RSA analysis results files from all subjects
        all_RSA_results = cell(length(subj_data_dirs), 1);
        for subj = 1:length(subj_data_dirs)
            curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
            if ~isempty(strfind(RSA_filename, 'objectviewing')) | ~isempty(strfind(RSA_filename, 'objview'))
                curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
                curr_filename = fullfile(curr_subj_objectviewing_analysis_dir, ['sm_', subj_data_dirs(subj).name, '_wholebrain', RSA_filename]);
                all_RSA_results{subj} = cosmo_fmri_dataset(curr_filename, 'mask', all_searchlight_mask_files{roi});
            elseif ~isempty(strfind(RSA_filename, 'jrd'))
                curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
                curr_filename = fullfile(curr_subj_jrd_analysis_dir, ['sm_', subj_data_dirs(subj).name, '_wholebrain', RSA_filename]);
                if exist(curr_filename, 'file')
                    all_RSA_results{subj} = cosmo_fmri_dataset(curr_filename, 'mask', all_searchlight_mask_files{roi});
                end
            end
        end
        
        % Organizing the across-subjects dataset
        all_RSA_results(cellfun(@isempty, all_RSA_results)) = [];
        [idxs, ds_cell_common] = cosmo_mask_dim_intersect(all_RSA_results);     % Choosing only voxels overlapping across all subjects
        ds_group = cosmo_stack(ds_cell_common, 1, 'drop_nonunique');
        ds_group.sa.targets = ones(size(ds_group.samples, 1),1);
        ds_group.sa.chunks = [1:size(ds_group.samples, 1)]';
        % Removing uninformative locations (NaN values, etc.)
        ds_group = cosmo_remove_useless_data(ds_group);
        % Removing locations not existing in half the subjects
        s = sum(~isnan(ds_group.samples) & ds_group.samples~=0);
        f = find(s > size(ds_group.samples,1) / 2);
        ds_group = cosmo_slice(ds_group, f, 2);
        % Running the group analysis
        h0_mean = 0;    % The null hypothesis - no correlation
        cluster_nbrhood = cosmo_cluster_neighborhood(ds_group);
        stat_map = cosmo_montecarlo_cluster_stat(ds_group, cluster_nbrhood, 'niter', num_iterations, 'h0_mean', h0_mean);
        % saving the result as a nifti file
        cosmo_map2fmri(stat_map, fullfile(group_analysis_dir, ['group_analysis_TFCE_sm_', num2str(num_iterations), '_iters', all_searchlight_mask_names{roi}, RSA_filename]));
    end
end

% DISPLAY IN CONNECTOME VIEWER - self threshold, 1.645 (p=0.05); palette red-yellow, fixed, min 1.64, max 2.33


%% Getting the average maps of individual whole-brain searchlights, and uncorrected t-test results (one-sample one-tailed t-test across subjects)
group_analysis_dir = fullfile(data_parent_dir, 'Analysis_group');
searchlight_averages_dir = fullfile(group_analysis_dir, 'Searchlight_map_averages');

RSA_filenames = {'_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_distances.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_segments.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_quadrants.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_schema.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_flipping.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_overlay.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_yaxis.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_grouping.nii',...
    '_wholebrain_avgruns_searchlight_objectviewing_day2minusday1_remapping.nii',...
    '_wholebrain_avgruns_searchlight_jrd_distances.nii',...
    '_wholebrain_avgruns_searchlight_jrd_segments.nii',...
    '_wholebrain_avgruns_searchlight_jrd_quadrants.nii',...
    '_wholebrain_avgruns_searchlight_jrd_schema.nii',...
    '_wholebrain_avgruns_searchlight_jrd_flipping.nii',...
    '_wholebrain_avgruns_searchlight_jrd_overlay.nii',...
    '_wholebrain_avgruns_searchlight_jrd_yaxis.nii',...
    '_wholebrain_avgruns_searchlight_jrd_grouping.nii',...
    '_wholebrain_avgruns_searchlight_jrd_remapping.nii',...
};

for i = 1:length(RSA_filenames)
    RSA_filename = RSA_filenames{i};
    disp(['Current file: ', RSA_filename])
    
    % Reading subject 1's searchlight result
    subj = 1;
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    if ~isempty(strfind(RSA_filename, 'objectviewing'))
        curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
        all_RSA_results = spm_read_vols(spm_vol(fullfile(curr_subj_objectviewing_analysis_dir, ['sm_', subj_data_dirs(subj).name, RSA_filename])));
    elseif ~isempty(strfind(RSA_filename, 'jrd'))
        curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
        if exist(fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, RSA_filename]), 'file')
            all_RSA_results = spm_read_vols(spm_vol(fullfile(curr_subj_jrd_analysis_dir, ['sm_', subj_data_dirs(subj).name, RSA_filename])));
        end
    end
    
    % Reading the RSA analysis results files from all subjects
    for subj = 2:length(subj_data_dirs)
        curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
        if ~isempty(strfind(RSA_filename, 'objectviewing'))
            curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
            all_RSA_results = cat(4, all_RSA_results, spm_read_vols(spm_vol(fullfile(curr_subj_objectviewing_analysis_dir, ['sm_', subj_data_dirs(subj).name, RSA_filename]))));
        elseif ~isempty(strfind(RSA_filename, 'jrd'))
            curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
            if exist(fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, RSA_filename]), 'file')
                all_RSA_results = cat(4, all_RSA_results, spm_read_vols(spm_vol(fullfile(curr_subj_jrd_analysis_dir, ['sm_', subj_data_dirs(subj).name, RSA_filename]))));
            end
        end
    end
    
    % Average of all subjects' RSA searchlight correlation values
    curr_average_map = nanmean(all_RSA_results, 4);
    save_mat_to_nifti(betas_file1, curr_average_map, fullfile(searchlight_averages_dir, ['average_map', RSA_filename]));
    
    % T-test across subjects
    size_data = size(all_RSA_results);
    all_RSA_results = reshape(all_RSA_results, [size_data(1)*size_data(2)*size_data(3), size_data(4)])';
    [~,p,~,ttest_all_subjects] = ttest(all_RSA_results, [],'tail','right');
    tvalues = ttest_all_subjects.tstat;
    tvalues = reshape(tvalues, size_data(1:3));
    save_mat_to_nifti(betas_file1, tvalues, fullfile(searchlight_averages_dir, ['ttest', RSA_filename]));
    
end

% SIGNIFICANT VALUE OF ONE-TAILED T-TEST FOR P=0.05 WITH 22 DFs - 1.717; FOR P=0.01 - 2.508