% Code to analyze the JRD data using single-trial beta estimates instead of an average pattern across the run

betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;
num_trials_per_run = 48;

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_EC'};
% all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC'};
%     all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview'};
%     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC','rHC','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview'};
%     all_ROI_names = {'bilat_HC', 'bilat_EC', 'bilat_RSC', 'bilat_PPA', 'bilat_OPA', 'bilat_LOC'};
%     all_ROI_names = {'bilat_HC_activation_jrd', 'bilat_HC_activation_objview', 'bilat_RSC', 'bilat_PPA', 'bilat_OPA', 'bilat_LOC'};
%     all_ROI_names = {'lRSC', 'rRSC', 'lPPA', 'rPPA', 'lOPA', 'rOPA', 'lLOC', 'rLOC','lHC_activation_jrd','rHC_activation_jrd','lHC_activation_objview','rHC_activation_objview'};
% %     all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC_activation_jrd','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview', 'bilat_HC_activation_objview'};
% %     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC_activation_jrd','rHC_activation_jrd','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview', 'lHC_activation_objview','rHC_activation_objview'};
num_ROIs = length(all_ROI_names);


% Locations of objects and their distances
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
mat_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
mat_distances = mat_distances / max(mat_distances(:));    % Normalizing to range 0-1
% Create the segment identity similarity matrix
mat_segments = ones(num_conditions);
mat_segments(1:4,1:4) = 0; mat_segments(13:16,13:16) = 0; mat_segments(5:12,5:12) = 0; mat_segments(13:16,1:4) = 0; mat_segments(1:4,13:16) = 0;
% Create the quadrant identity similarity matrix
mat_quadrants = ones(num_conditions);
mat_quadrants(1:4,1:4) = 0; mat_quadrants(5:8,5:8) = 0; mat_quadrants(9:12,9:12) = 0; mat_quadrants(13:16,13:16) = 0;
% Create the orthogonal segment identity similarity matrix
mat_orth_segments = ones(num_conditions);
mat_orth_segments(1:8,1:8) = 0; mat_orth_segments(9:16,9:16) = 0;
% Converting all matrices from dissimilarity to similarity matrices
mat_distances = 1 - mat_distances; mat_segments = 1 - mat_segments; mat_quadrants = 1 - mat_quadrants; mat_orth_segments = 1 - mat_orth_segments;


% Overlay - with distance
locations_new_overlay = locations_all_objects;
locations_new_overlay(5:12,2) = locations_new_overlay(5:12,2) + 75;
mat_overlay = pdist2(locations_new_overlay, locations_new_overlay);         % Distances between all objects
mat_overlay = mat_overlay / max(mat_overlay(:));    % Normalizing to range 0-1
mat_overlay = 1 - mat_overlay;

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
mat_grouping = nan(num_conditions);
mat_grouping(mat_within_seg_adj_quad == 1) = 1;
mat_grouping(mat_between_seg_adj_quad == 1) = -1;


% Defining variables for MDS
allsubjs_pattern_dist = cell(size(all_ROI_names));
allsubjs_pattern_dist_withinsegtrials = cell(size(all_ROI_names));
allsubjs_pattern_dist_betweensegtrials = cell(size(all_ROI_names));
for i=1:length(all_ROI_names)
    allsubjs_pattern_dist{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
    allsubjs_pattern_dist_withinsegtrials{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
    allsubjs_pattern_dist_betweensegtrials{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
end



% Defining variables
corr_object1dist = nan(length(subj_data_dirs), num_ROIs);
corr_object1dist_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_object1dist_grouping = nan(length(subj_data_dirs), num_ROIs);
corr_object1dist_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_object1dist_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_object1dist_remapping = nan(length(subj_data_dirs), num_ROIs);
corr_object1_identity = nan(length(subj_data_dirs), num_ROIs);
corr_segments = nan(length(subj_data_dirs), num_ROIs);
corr_orth_segments = nan(length(subj_data_dirs), num_ROIs);
corr_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_viewdirdiff = nan(length(subj_data_dirs), num_ROIs);
corr_viewdirdiff_withinquadrant = nan(length(subj_data_dirs), num_ROIs);
corr_viewdirdiff_withinsegment_adjquadrant = nan(length(subj_data_dirs), num_ROIs);
corr_viewdirdiff_betweensegments_adjquadrant = nan(length(subj_data_dirs), num_ROIs);
corr_viewdirdiff_diagonalquadrant = nan(length(subj_data_dirs), num_ROIs);
corr_viewdirdiff_withinsegment = nan(length(subj_data_dirs), num_ROIs);
corr_viewdir_identity = nan(length(subj_data_dirs), num_ROIs);
corr_viewdir_identity_bin = nan(length(subj_data_dirs), num_ROIs);
corr_withinsegtrials_object1dist = nan(length(subj_data_dirs), num_ROIs);
corr_withinsegtrials_object1dist_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_withinsegtrials_object1dist_grouping = nan(length(subj_data_dirs), num_ROIs);
corr_withinsegtrials_object1dist_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_withinsegtrials_object1dist_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_withinsegtrials_object1dist_remapping = nan(length(subj_data_dirs), num_ROIs);
corr_betweensegtrials_object1dist = nan(length(subj_data_dirs), num_ROIs);
corr_betweensegtrials_object1dist_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_betweensegtrials_object1dist_grouping = nan(length(subj_data_dirs), num_ROIs);
corr_betweensegtrials_object1dist_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_betweensegtrials_object1dist_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_betweensegtrials_object1dist_remapping = nan(length(subj_data_dirs), num_ROIs);


for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    curr_subj_jrd_singletrials_analysis_dirs = fullfile(curr_subj_jrd_analysis_dirs, 'single_trials');
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    num_runs_jrd = length(dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat')));
    
    
    % Creating a matrix representing the runs (within run is nan, across runs 1)
    runs_mat = ones(num_trials_per_run * num_runs_jrd);
    for i = 1:num_runs_jrd
        runs_mat((i-1)*num_trials_per_run+1:i*num_trials_per_run, (i-1)*num_trials_per_run+1:i*num_trials_per_run) = nan;
    end
    
    
    %% Loading the subject specific ROIs and the JRD data
    % Reading the ROIs
    all_ROIs_individual = cell(length(all_ROI_names));
    for i = 1:length(all_ROI_names)
        all_ROIs_individual{i} = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
    end
    num_ROIs_individual = length(all_ROIs_individual);
    
    % Reading all betas from this subject's JRD folders
    all_betas = [];
    for r = 1:num_runs_jrd
        temp_files = dir(fullfile(curr_subj_jrd_singletrials_analysis_dirs{r}, 'singletrial_beta*.nii'));
        if ~isempty(temp_files)
            for i = 1:length(temp_files)
                all_betas = cat(4, all_betas, spm_read_vols(spm_vol(fullfile(temp_files(i).folder, temp_files(i).name))));
            end
        end
    end
    all_betas = reshape(all_betas, [size(all_betas, 1)*size(all_betas, 2)*size(all_betas, 3), size(all_betas, 4)]);
    
    
    if ~isempty(all_betas)        % JRD runs exist
        %% Getting the information on objects and directions in each trial and creating dissimilarity matrices
        all_curr_subj_info = [];
        all_object1_start = [];
        all_object3_target = [];
        all_viewing_direction = [];
        all_quadrant1_start = [];
        all_quadrant2_target = [];
        for r = 1:num_runs_jrd
            curr_info_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, 'all_info*.mat'));
            curr_info = load(fullfile(curr_info_file.folder, curr_info_file.name));
            all_curr_subj_info = [all_curr_subj_info;curr_info.all_info];
            all_object1_start = [all_object1_start; all_curr_subj_info(r).object1_unity_location(1:num_trials_per_run) + 1];
            all_object3_target = [all_object3_target; all_curr_subj_info(r).object3_unity_location(1:num_trials_per_run) + 1];
            all_viewing_direction = [all_viewing_direction; all_curr_subj_info(r).viewing_direction(1:num_trials_per_run)'];
            all_quadrant1_start = [all_quadrant1_start; all_curr_subj_info(r).quadrant_start(1:num_trials_per_run)];
            all_quadrant2_target = [all_quadrant2_target; all_curr_subj_info(r).quadrant_target(1:num_trials_per_run)];
        end
        % Figuring out the segment for each trial
        all_segment1_start = (all_quadrant1_start == 2 | all_quadrant1_start == 3);
        all_segment2_target = (all_quadrant2_target == 2 | all_quadrant2_target == 3);
        all_orth_segment1_start = (all_quadrant1_start == 3 | all_quadrant1_start == 4);
        all_orth_segment2 = (all_quadrant2_target == 3 | all_quadrant2_target == 4);
        % Viewing direction of each trial - lowering the resolution to bins
        all_viewing_direction_bins = nan(size(all_viewing_direction));
        all_viewing_direction_bins(all_viewing_direction>-45 & all_viewing_direction<=45) = 1;      % North - toward mountains
        all_viewing_direction_bins(all_viewing_direction>45 & all_viewing_direction<=135) = 2;      % West
        all_viewing_direction_bins(all_viewing_direction>135 | all_viewing_direction<=-135) = 3;    % South
        all_viewing_direction_bins(all_viewing_direction>-135 & all_viewing_direction<=-45) = 4;    % East
        
        
        % Calculating matrices for all trials differences in distance, viewing direction, etc.
        curr_distances_all_object1 = nan(length(all_object1_start));
        curr_distances_all_object1_overlay = nan(length(all_object1_start));
        curr_distances_all_object1_grouping = nan(length(all_object1_start));
        curr_distances_all_object1_within_segment = nan(length(all_object1_start));
        curr_distances_all_object1_between_segments = nan(length(all_object1_start));
        curr_same_object1 = nan(length(all_object1_start));
        curr_same_segment = nan(length(all_object1_start));
        curr_same_orth_segment = nan(length(all_object1_start));
        curr_same_quadrant = nan(length(all_object1_start));
        curr_viewdir_diff_all_objects = nan(length(all_object1_start));
        curr_same_viewdir = nan(length(all_object1_start));
        curr_same_viewdir_bin = nan(length(all_object1_start));
        
        for i = 1:length(all_object1_start)
            for j = 1:length(all_object1_start)
                if i ~= j
                    % Object distance coding and segment / quadrant coding
                    if all_object1_start(i) ~= all_object1_start(j)     % If this is not the same starting object
                        % Distance between objects in different trials
                        curr_distances_all_object1(i,j) = mat_distances(all_object1_start(i), all_object1_start(j));
                        curr_distances_all_object1_overlay(i,j) = mat_overlay(all_object1_start(i), all_object1_start(j));
                        curr_distances_all_object1_grouping(i,j) = mat_grouping(all_object1_start(i), all_object1_start(j));
                        curr_distances_all_object1_within_segment(i,j) = mat_dist_within_segment(all_object1_start(i), all_object1_start(j));
                        curr_distances_all_object1_between_segments(i,j) = mat_dist_between_segments(all_object1_start(i), all_object1_start(j));
                        % Same / different segment or quadrant between trials
                        curr_same_segment(i,j) = (all_segment1_start(i) == all_segment1_start(j));
                        curr_same_orth_segment(i,j) = (all_orth_segment1_start(i) == all_orth_segment1_start(j));
                        curr_same_quadrant(i,j) = (all_quadrant1_start(i) == all_quadrant1_start(j));
                        curr_same_object1(i,j) = 0;
                    else
                        curr_same_object1(i,j) = 1;
                    end
                    
                    % Viewing direction differences between trials
                    % Calculating the angle between the viewing directions
                    curr_angle_diff = all_viewing_direction(i) - all_viewing_direction(j);
                    if abs(curr_angle_diff) > 180, curr_angle_diff = curr_angle_diff - 360*sign(curr_angle_diff); end
                    if curr_angle_diff < 0, curr_angle_diff = 360 + curr_angle_diff; end
                    if curr_angle_diff > 180, curr_angle_diff = 360 - curr_angle_diff; end     % Choosing the smaller of the two possible angles
                    curr_viewdir_diff_all_objects(i,j) = curr_angle_diff;
                    % Is it the same viewing direction, or viewing direction bin
                    curr_same_viewdir(i,j) = (curr_angle_diff == 0);
                    curr_same_viewdir_bin(i,j) = (all_viewing_direction_bins(i) == all_viewing_direction_bins(j));
                end
            end
        end
        
        % Converting the viewing direction matrix to similarity matrix (similar angles -> high value, dissimilar angles - low values)
        curr_viewdir_diff_all_objects = curr_viewdir_diff_all_objects / max(curr_viewdir_diff_all_objects(:));
        curr_viewdir_diff_all_objects = 1 - curr_viewdir_diff_all_objects;
        
        % Direction coding only inside quadrant, segment, etc.
        curr_viewdir_diff_withinquadrant = curr_viewdir_diff_all_objects; curr_viewdir_diff_withinquadrant(curr_same_quadrant == 0) = nan;
        curr_viewdir_diff_withinsegment_adjquadrant = curr_viewdir_diff_all_objects; curr_viewdir_diff_withinsegment_adjquadrant(curr_same_quadrant ~= 0 | curr_same_segment == 0) = nan;
        curr_viewdir_diff_betweensegments_adjquadrant = curr_viewdir_diff_all_objects; curr_viewdir_diff_betweensegments_adjquadrant(curr_same_quadrant ~= 0 | curr_same_orth_segment == 0) = nan;
        curr_viewdir_diff_diagonalquadrant = curr_viewdir_diff_all_objects; curr_viewdir_diff_diagonalquadrant(curr_same_orth_segment == 1 | curr_same_segment == 1 | isnan(curr_same_quadrant)) = nan;
        curr_viewdir_diff_withinsegment = curr_viewdir_diff_all_objects; curr_viewdir_diff_withinsegment(curr_same_segment == 0) = nan;
        
        % Getting the non-diagonal elements of each matrix
        curr_distances_all_object1_nondiag = curr_distances_all_object1(find(tril(ones(length(all_object1_start)),-1)));
        curr_distances_all_object1_overlay_nondiag = curr_distances_all_object1_overlay(find(tril(ones(length(all_object1_start)),-1)));
        curr_distances_all_object1_grouping_nondiag = curr_distances_all_object1_grouping(find(tril(ones(length(all_object1_start)),-1)));        curr_angles_viewdir_diff_all_objects_nondiag = curr_viewdir_diff_all_objects(find(tril(ones(length(all_object1_start)),-1)));
        curr_distances_all_object1_within_segment_nondiag = curr_distances_all_object1_within_segment(find(tril(ones(length(all_object1_start)),-1)));
        curr_distances_all_object1_between_segments_nondiag = curr_distances_all_object1_between_segments(find(tril(ones(length(all_object1_start)),-1)));
        curr_same_object1_nondiag = curr_same_object1(find(tril(ones(length(all_object1_start)),-1)));
        curr_same_segment_nondiag = curr_same_segment(find(tril(ones(length(all_object1_start)),-1)));
        curr_same_orth_segment_nondiag = curr_same_orth_segment(find(tril(ones(length(all_object1_start)),-1)));
        curr_same_quadrant_nondiag = curr_same_quadrant(find(tril(ones(length(all_object1_start)),-1)));
        curr_viewdir_diff_all_objects_nondiag = curr_viewdir_diff_all_objects(find(tril(ones(length(all_object1_start)),-1)));
        curr_viewdir_diff_betweensegments_adjquadrant_nondiag = curr_viewdir_diff_betweensegments_adjquadrant(find(tril(ones(length(all_object1_start)),-1)));
        curr_viewdir_diff_withinquadrant_nondiag = curr_viewdir_diff_withinquadrant(find(tril(ones(length(all_object1_start)),-1)));
        curr_viewdir_diff_withinsegment_adjquadrant_nondiag = curr_viewdir_diff_withinsegment_adjquadrant(find(tril(ones(length(all_object1_start)),-1)));
        curr_viewdir_diff_diagonalquadrant_nondiag = curr_viewdir_diff_diagonalquadrant(find(tril(ones(length(all_object1_start)),-1)));
        curr_viewdir_diff_withinsegment_nondiag = curr_viewdir_diff_withinsegment(find(tril(ones(length(all_object1_start)),-1)));
        curr_same_viewdir_nondiag = curr_same_viewdir(find(tril(ones(length(all_object1_start)),-1)));
        curr_same_viewdir_bin_nondiag = curr_same_viewdir_bin(find(tril(ones(length(all_object1_start)),-1)));
        runs_mat_nondiag = runs_mat(find(tril(ones(length(all_object1_start)),-1)));

%         % Getting only the between-run elements (Mumford, Davis & Poldrack, NeuroImage 2014)
%         curr_distances_all_object1_nondiag = curr_distances_all_object1_nondiag(runs_mat_nondiag==1);
%         curr_distances_all_object1_overlay_nondiag = curr_distances_all_object1_overlay_nondiag(runs_mat_nondiag==1);
%         curr_distances_all_object1_grouping_nondiag = curr_distances_all_object1_grouping_nondiag(runs_mat_nondiag==1);
%         curr_distances_all_object1_within_segment_nondiag = curr_distances_all_object1_within_segment_nondiag(runs_mat_nondiag==1);
%         curr_distances_all_object1_between_segments_nondiag = curr_distances_all_object1_between_segments_nondiag(runs_mat_nondiag==1);
%         curr_same_object1_nondiag = curr_same_object1_nondiag(runs_mat_nondiag==1);
%         curr_same_segment_nondiag = curr_same_segment_nondiag(runs_mat_nondiag==1);
%         curr_same_orth_segment_nondiag = curr_same_orth_segment_nondiag(runs_mat_nondiag==1);
%         curr_same_quadrant_nondiag = curr_same_quadrant_nondiag(runs_mat_nondiag==1);
%         curr_viewdir_diff_all_objects_nondiag = curr_viewdir_diff_all_objects_nondiag(runs_mat_nondiag==1);
%         curr_viewdir_diff_betweensegments_adjquadrant_nondiag = curr_viewdir_diff_betweensegments_adjquadrant_nondiag(runs_mat_nondiag==1);
%         curr_viewdir_diff_withinquadrant_nondiag = curr_viewdir_diff_withinquadrant_nondiag(runs_mat_nondiag==1);
%         curr_viewdir_diff_withinsegment_adjquadrant_nondiag = curr_viewdir_diff_withinsegment_adjquadrant_nondiag(runs_mat_nondiag==1);
%         curr_viewdir_diff_diagonalquadrant_nondiag = curr_viewdir_diff_diagonalquadrant_nondiag(runs_mat_nondiag==1);
%         curr_viewdir_diff_withinsegment_nondiag = curr_viewdir_diff_withinsegment_nondiag(runs_mat_nondiag==1);
%         curr_same_viewdir_nondiag = curr_same_viewdir_nondiag(runs_mat_nondiag==1);
%         curr_same_viewdir_bin_nondiag = curr_same_viewdir_bin_nondiag(runs_mat_nondiag==1);
        
        
        
        
        %         %% Reading the data - average pattern in ROI for each object and direction
        %
        %         % Getting the average pattern for each starting location and viewing direction (16*4 = 64 patterns)
        %         average_patterns = nan(size(all_betas, 1), size(unique_o1_viewdir, 1));
        %         for i = 1:size(unique_o1_viewdir, 1)
        %             average_patterns(:,i) = nanmean(all_betas(:, find(all_o1_viewdir(:,1)==unique_o1_viewdir(i,1) & all_o1_viewdir(:,2)==unique_o1_viewdir(i,2))), 2);
        %         end
        %         % Removing the mean pattern (mean of the 48 average patterns)
        %         mean_pattern_all = nanmean(average_patterns, 2);
        %         for i = 1:length(unique_o1_viewdir), average_patterns(:,i) = average_patterns(:,i) - mean_pattern_all; end
        %
        %
        %         % Getting the average pattern for each starting location and viewing direction - FOR WITHIN SEGMENT TRIALS ONLY
        %         average_patterns_samesegquestions = nan(size(all_betas, 1), size(unique_o1_viewdir, 1));
        %         for i = 1:size(unique_o1_viewdir, 1)
        %             average_patterns_samesegquestions(:,i) = nanmean(all_betas(:, find(all_o1_viewdir(:,1)==unique_o1_viewdir(i,1) & all_o1_viewdir(:,2)==unique_o1_viewdir(i,2) & all_is_same_segment == 1)), 2);
        %         end
        %         % Removing the mean pattern (mean of the 48 average patterns)
        %         mean_pattern_all = nanmean(average_patterns_samesegquestions, 2);
        %         for i = 1:length(unique_o1_viewdir), average_patterns_samesegquestions(:,i) = average_patterns_samesegquestions(:,i) - mean_pattern_all; end
        
        
        
        %% MVPA analysis in each ROI
        for roi = 1:num_ROIs_individual        % Going over ROIs
            curr_ROI = find(all_ROIs_individual{roi}(:)==1);
            if ~isempty(curr_ROI)
                % Reading the ROI data and calculating correlations
                curr_betas = all_betas(curr_ROI, :);
                % Removing the mean pattern from each run
                for r = 1:num_runs_jrd
                    temp_range = (r-1)*num_trials_per_run+1:r*num_trials_per_run;
                    mean_pattern_all = nanmean(curr_betas(:, temp_range), 2);
                    for i = 1:length(temp_range), curr_betas(:,temp_range(i)) = curr_betas(:,temp_range(i)) - mean_pattern_all; end
                end
                % Computing the correlation matrix in this ROI between all betas, and getting its non-diagonal elements
                curr_corr_betas = fisherz(corr(curr_betas, 'rows', 'complete'));
                curr_corr_betas_nondiag = curr_corr_betas(find(tril(ones(size(curr_corr_betas)),-1)));
%                 % Getting only the between-run elements
%                 curr_corr_betas_nondiag = curr_corr_betas_nondiag(runs_mat_nondiag == 1);
                
                % Reading the ROI data and calculating correlations - SAME SEGMENT TRIALS ONLY
                curr_betas_withinsegtrials = curr_betas;
                curr_betas_betweensegtrials = curr_betas;
                % Removing trials where the question was between segments
                curr_betas_withinsegtrials(:, all_segment1_start ~= all_segment2_target) = nan;
                curr_betas_betweensegtrials(:, all_segment1_start == all_segment2_target) = nan;
                % Removing the mean pattern from each run
                for r=1:num_runs_jrd
                    temp_range = (r-1)*num_trials_per_run+1:r*num_trials_per_run;
                    mean_pattern_all = nanmean(curr_betas_withinsegtrials(:, temp_range), 2);
                    for i = 1:length(temp_range), curr_betas_withinsegtrials(:,temp_range(i)) = curr_betas_withinsegtrials(:,temp_range(i)) - mean_pattern_all; end
                    mean_pattern_all = nanmean(curr_betas_betweensegtrials(:, temp_range), 2);
                    for i = 1:length(temp_range), curr_betas_betweensegtrials(:,temp_range(i)) = curr_betas_betweensegtrials(:,temp_range(i)) - mean_pattern_all; end
                end
                %                 mean_pattern_all = nanmean(curr_betas_withinsegtrials, 2);
                %                 for i = 1:length(all_object1), curr_betas_withinsegtrials(:,i) = curr_betas_withinsegtrials(:,i) - mean_pattern_all; end
                % Computing correlation between all betas
                curr_corr_betas_withinsegtrials = nan(size(curr_betas_withinsegtrials, 2));
                curr_corr_betas_betweensegtrials = nan(size(curr_betas_betweensegtrials, 2));
                for i=1:size(curr_betas_withinsegtrials, 2)
                    for j=1:size(curr_betas_withinsegtrials, 2)
                        curr_corr_betas_withinsegtrials(i,j) = fisherz(corr(curr_betas_withinsegtrials(:,i),curr_betas_withinsegtrials(:,j),'rows','complete'));
                    end
                end
                for i=1:size(curr_betas_betweensegtrials, 2)
                    for j=1:size(curr_betas_betweensegtrials, 2)
                        curr_corr_betas_betweensegtrials(i,j) = fisherz(corr(curr_betas_betweensegtrials(:,i),curr_betas_betweensegtrials(:,j),'rows','complete'));
                    end
                end
                curr_corr_betas_withinsegtrials_nondiag = curr_corr_betas_withinsegtrials(find(tril(ones(size(curr_corr_betas_withinsegtrials)),-1)));
                curr_corr_betas_betweensegtrials_nondiag = curr_corr_betas_betweensegtrials(find(tril(ones(size(curr_corr_betas_betweensegtrials)),-1)));
%                 % Getting only the between-run elements
%                 curr_corr_betas_withinsegtrials_nondiag = curr_corr_betas_withinsegtrials_nondiag(runs_mat_nondiag == 1);
%                 curr_corr_betas_betweensegtrials_nondiag = curr_corr_betas_betweensegtrials_nondiag(runs_mat_nondiag == 1);
                
                
                % Correlation of pattern dissimilarity matrix to the different model matrices
                corr_object1dist(subj,roi) = corr(curr_corr_betas_nondiag, curr_distances_all_object1_nondiag, 'type','Spearman','rows','complete');
                corr_object1dist_overlay(subj,roi) = corr(curr_corr_betas_nondiag, curr_distances_all_object1_overlay_nondiag, 'type','Spearman','rows','complete');
                corr_object1dist_grouping(subj, roi) = corr(curr_corr_betas_nondiag, curr_distances_all_object1_grouping_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_object1dist_within_segment(subj, roi) = corr(curr_corr_betas_nondiag, curr_distances_all_object1_within_segment_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_object1dist_between_segments(subj, roi) = corr(curr_corr_betas_nondiag, curr_distances_all_object1_between_segments_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_object1dist_remapping(subj, roi) = corr_object1dist_within_segment(subj, roi) - corr_object1dist_between_segments(subj, roi);
                corr_object1_identity(subj,roi) = corr(curr_corr_betas_nondiag, curr_same_object1_nondiag, 'type','Spearman','rows','complete');
                corr_segments(subj,roi) = corr(curr_corr_betas_nondiag, curr_same_segment_nondiag, 'type','Spearman','rows','complete');
                corr_orth_segments(subj,roi) = corr(curr_corr_betas_nondiag, curr_same_orth_segment_nondiag, 'type','Spearman','rows','complete');
                corr_quadrants(subj,roi) = corr(curr_corr_betas_nondiag, curr_same_quadrant_nondiag, 'type','Spearman','rows','complete');
                corr_viewdirdiff(subj,roi) = -1 * corr(curr_corr_betas_nondiag, curr_viewdir_diff_all_objects_nondiag, 'type','Spearman','rows','complete');
                corr_viewdirdiff_withinquadrant(subj,roi) = corr(curr_corr_betas_nondiag, curr_viewdir_diff_withinquadrant_nondiag, 'type','Spearman','rows','complete');
                corr_viewdirdiff_withinsegment_adjquadrant(subj,roi) = corr(curr_corr_betas_nondiag, curr_viewdir_diff_withinsegment_adjquadrant_nondiag, 'type','Spearman','rows','complete');
                corr_viewdirdiff_betweensegments_adjquadrant(subj,roi) = corr(curr_corr_betas_nondiag, curr_viewdir_diff_betweensegments_adjquadrant_nondiag, 'type','Spearman','rows','complete');
                corr_viewdirdiff_diagonalquadrant(subj,roi) = corr(curr_corr_betas_nondiag, curr_viewdir_diff_diagonalquadrant_nondiag, 'type','Spearman','rows','complete');
                corr_viewdirdiff_withinsegment(subj,roi) = corr(curr_corr_betas_nondiag, curr_viewdir_diff_withinsegment_nondiag, 'type','Spearman','rows','complete');
                corr_viewdir_identity(subj,roi) = corr(curr_corr_betas_nondiag, curr_same_viewdir_nondiag, 'type','Spearman','rows','complete');
                corr_viewdir_identity_bin(subj,roi) = corr(curr_corr_betas_nondiag, curr_same_viewdir_bin_nondiag, 'type','Spearman','rows','complete');
                % Correlations - within-segment trials only
                corr_withinsegtrials_object1dist(subj,roi) = corr(curr_corr_betas_withinsegtrials_nondiag, curr_distances_all_object1_nondiag, 'type','Spearman','rows','complete');
                corr_withinsegtrials_object1dist_overlay(subj,roi) = corr(curr_corr_betas_withinsegtrials_nondiag, curr_distances_all_object1_overlay_nondiag, 'type','Spearman','rows','complete');
                corr_withinsegtrials_object1dist_grouping(subj, roi) = corr(curr_corr_betas_withinsegtrials_nondiag, curr_distances_all_object1_grouping_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_withinsegtrials_object1dist_within_segment(subj, roi) = corr(curr_corr_betas_withinsegtrials_nondiag, curr_distances_all_object1_within_segment_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_withinsegtrials_object1dist_between_segments(subj, roi) = corr(curr_corr_betas_withinsegtrials_nondiag, curr_distances_all_object1_between_segments_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_withinsegtrials_object1dist_remapping(subj, roi) = corr_withinsegtrials_object1dist_within_segment(subj, roi) - corr_withinsegtrials_object1dist_between_segments(subj, roi);
                corr_betweensegtrials_object1dist(subj,roi) = corr(curr_corr_betas_betweensegtrials_nondiag, curr_distances_all_object1_nondiag, 'type','Spearman','rows','complete');
                corr_betweensegtrials_object1dist_overlay(subj,roi) = corr(curr_corr_betas_betweensegtrials_nondiag, curr_distances_all_object1_overlay_nondiag, 'type','Spearman','rows','complete');
                corr_betweensegtrials_object1dist_grouping(subj, roi) = corr(curr_corr_betas_betweensegtrials_nondiag, curr_distances_all_object1_grouping_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_betweensegtrials_object1dist_within_segment(subj, roi) = corr(curr_corr_betas_betweensegtrials_nondiag, curr_distances_all_object1_within_segment_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_betweensegtrials_object1dist_between_segments(subj, roi) = corr(curr_corr_betas_betweensegtrials_nondiag, curr_distances_all_object1_between_segments_nondiag, 'Type', 'Spearman', 'rows', 'complete');
                corr_betweensegtrials_object1dist_remapping(subj, roi) = corr_betweensegtrials_object1dist_within_segment(subj, roi) - corr_betweensegtrials_object1dist_between_segments(subj, roi);
                
                
                % Getting the patterns for MDS
                % Averaging per object
                curr_patterns = []; curr_patterns_withinsegtrials = []; curr_patterns_betweensegtrials = [];
                for i = 1:num_conditions
                    curr_patterns(:,i) = nanmean(curr_betas(:,all_object1_start==i),2);
                    curr_patterns_withinsegtrials(:,i) = nanmean(curr_betas_withinsegtrials(:,all_object1_start==i),2);
                    curr_patterns_betweensegtrials(:,i) = nanmean(curr_betas_betweensegtrials(:,all_object1_start==i),2);
                end
                pattern_dist = 1 - fisherz(corr(curr_patterns, 'rows', 'complete'));
                pattern_dist_withinsegtrials = 1 - fisherz(corr(curr_patterns_withinsegtrials, 'rows', 'complete'));
                pattern_dist_betweensegtrials = 1 - fisherz(corr(curr_patterns_betweensegtrials, 'rows', 'complete'));
                allsubjs_pattern_dist{roi}(:,:,subj) = pattern_dist;
                allsubjs_pattern_dist_withinsegtrials{roi}(:,:,subj) = pattern_dist_withinsegtrials;
                allsubjs_pattern_dist_betweensegtrials{roi}(:,:,subj) = pattern_dist_betweensegtrials;
            end
        end
        
        
        
        %         %% Searchlight
        %         searchlight_radius = 3;
        %
        %         % Getting the grey-matter mask from freesurfer parcellation
        %         curr_subject_functional_directory = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'fmriprep');
        %         curr_subject_functional_directory = fullfile(fullfile(curr_subject_functional_directory, subj_data_dirs(subj).name), 'ses-day1/func/');
        %         curr_parcellation_file = dir(fullfile(curr_subject_functional_directory, '*-aseg_dseg.nii'));
        %         curr_parcellation = reslice_data(fullfile(curr_subject_functional_directory, curr_parcellation_file(1).name), betas_file1, 0);
        %         wholebrain_ROI = (curr_parcellation == 3 | curr_parcellation == 17 | curr_parcellation == 18 | curr_parcellation == 42 | curr_parcellation == 53 | curr_parcellation == 54);       % Cortex, hippocampus, amygdala
        %         curr_ROI = wholebrain_ROI;
        %         curr_ROI_voxels = find(curr_ROI);   % Getting all of the indices of voxels in the ROI
        %
        %         % Extracting the beta values from the grey-matter mask
        %         all_betas_ROI_jrd = average_patterns;
        %
        %         ROI_searchlight_results_jrd_distances = nan(size(curr_ROI));
        %         ROI_searchlight_results_jrd_segments = nan(size(curr_ROI));
        %         ROI_searchlight_results_jrd_quadrants = nan(size(curr_ROI));
        %         ROI_searchlight_results_jrd_viewdirections = nan(size(curr_ROI));
        %         for i = 1:length(curr_ROI_voxels)
        %             % Defining a sphere around the voxel according to the searchlight radius
        %             [curr_vox_coord_x, curr_vox_coord_y, curr_vox_coord_z] = ind2sub(size(curr_ROI), curr_ROI_voxels(i));        % Getting the 3D coordinate of the current voxel
        %             [x y z] = ind2sub(size(curr_ROI),1:numel(curr_ROI));      % Getting the 3D coordinates of all voxels
        %             M = [x', y', z'] - repmat([curr_vox_coord_x, curr_vox_coord_y, curr_vox_coord_z], numel(curr_ROI), 1);     % The distance along each dimension of each voxel from the current one
        %             D = (sqrt(sum(M.^2,2)) <= searchlight_radius);        % Identifying only the locations whose sum of squared distances from the voxels are equal to or smaller than the searchlight radius
        %             D = reshape(D, size(curr_ROI));        % Reshaping the sphere locations to the original matrix size
        %
        %             % Masking with the ROI to get only voxels that are in it and in the searchlight radius
        %             curr_searchlight_voxels = find(D ~= 0 & curr_ROI ~= 0);
        %
        %             % calculating correlation between patterns in the current ROI sphere across runs
        %             corr_ROI_jrd = fisherz(corr(all_betas_ROI_jrd(curr_searchlight_voxels,:), all_betas_ROI_jrd(curr_searchlight_voxels,:)));      % Correlation between runs i and j
        %
        %             if ~isempty(corr_ROI_jrd)
        %                 % Getting patterns and correlating to different matrices
        %                 corr_ROI_jrd_nondiag = corr_ROI_jrd(find(tril(ones(size(all_betas_ROI_jrd, 2)),-1)));   % Taking only the lower triangle elements as a vector
        %                 curr_dist_temp = curr_distances_all_object1(find(tril(ones(size(all_betas_ROI_jrd, 2)),-1)));
        %                 curr_seg_temp = curr_same_segment(find(tril(ones(size(all_betas_ROI_jrd, 2)),-1)));
        %                 curr_quad_temp = curr_same_quadrant(find(tril(ones(size(all_betas_ROI_jrd, 2)),-1)));
        %                 curr_view_temp = curr_angles_viewdirections_diff_all_objects(find(tril(ones(size(all_betas_ROI_jrd, 2)),-1)));
        %                 ROI_searchlight_results_jrd_distances(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(~isnan(curr_dist_temp)), curr_dist_temp(~isnan(curr_dist_temp)), 'Type', 'Spearman');
        %                 ROI_searchlight_results_jrd_segments(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(~isnan(curr_seg_temp)), curr_seg_temp(~isnan(curr_seg_temp)), 'Type', 'Spearman');
        %                 ROI_searchlight_results_jrd_quadrants(curr_ROI_voxels(i)) = corr(corr_ROI_jrd_nondiag(~isnan(curr_quad_temp)), curr_quad_temp(~isnan(curr_quad_temp)), 'Type', 'Spearman');
        %                 ROI_searchlight_results_jrd_viewdirections(curr_ROI_voxels(i)) = -1 * corr(corr_ROI_jrd_nondiag(~isnan(curr_view_temp)), curr_view_temp(~isnan(curr_view_temp)), 'Type', 'Spearman');
        %             end
        %         end
        %
        %         % Saving the searchlight results
        %         save_mat_to_nifti(betas_file1, ROI_searchlight_results_jrd_distances, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_searchlight_jrd_singletrials_avgruns_distances.nii']));
        %         save_mat_to_nifti(betas_file1, ROI_searchlight_results_jrd_segments, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_searchlight_jrd_singletrials_avgruns_segments.nii']));
        %         save_mat_to_nifti(betas_file1, ROI_searchlight_results_jrd_quadrants, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_searchlight_jrd_singletrials_avgruns_quadrants.nii']));
        %         save_mat_to_nifti(betas_file1, ROI_searchlight_results_jrd_viewdirections, fullfile(curr_subj_jrd_analysis_dir, [subj_data_dirs(subj).name, '_searchlight_jrd_singletrials_avgruns_viewdirections.nii']));
        
    end
end



curr_corr = corr_object1dist; [~,p_corr_object1dist] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_object1dist_overlay; [~,p_corr_object1dist_overlay] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_object1dist_grouping; [~,p_corr_object1dist_grouping] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_object1dist_remapping; [~,p_corr_object1dist_remapping] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_object1_identity; [~,p_corr_object1_identity] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_segments; [~,p_corr_segments] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_orth_segments; [~,p_corr_orth_segments] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_quadrants; [~,p_corr_quadrants] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdirdiff; [~,p_corr_viewdirdiff] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdirdiff_withinquadrant; [~,p_corr_viewdirdiff_withinquadrant] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdirdiff_withinsegment_adjquadrant; [~,p_corr_viewdirdiff_withinsegment_adjquadrant] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdirdiff_betweensegments_adjquadrant; [~,p_corr_viewdirdiff_betweensegments_adjquadrant] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdirdiff_diagonalquadrant; [~,p_corr_viewdirdiff_diagonalquadrant] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdirdiff_withinsegment; [~,p_corr_viewdirdiff_withinsegment] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdir_identity; [~,p_corr_viewdir_identity] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_viewdir_identity_bin; [~,p_corr_viedir_bin_identity] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');

curr_corr = corr_withinsegtrials_object1dist; [~,p_corr_withinsegtrials_object1dist] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_withinsegtrials_object1dist_overlay; [~,p_corr_withinsegtrials_object1dist_overlay] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_withinsegtrials_object1dist_grouping; [~,p_corr_withinsegtrials_object1dist_grouping] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_withinsegtrials_object1dist_remapping; [~,p_corr_withinsegtrials_object1dist_remapping] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_betweensegtrials_object1dist; [~,p_corr_betweensegtrials_object1dist] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_betweensegtrials_object1dist_overlay; [~,p_corr_betweensegtrials_object1dist_overlay] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_betweensegtrials_object1dist_grouping; [~,p_corr_betweensegtrials_object1dist_grouping] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');
curr_corr = corr_betweensegtrials_object1dist_remapping; [~,p_corr_betweensegtrials_object1dist_remapping] = ttest(curr_corr(sum(curr_corr, 2)~=0, :), [], 'tail', 'right');


%% MDS

% Summarizing pattern distances
for roi=1:length(all_ROI_names)
    allsubjs_mean_pattern_dist{roi} = nanmean(allsubjs_pattern_dist{roi},3);
    allsubjs_mean_pattern_dist_withinsegtrials{roi} = nanmean(allsubjs_pattern_dist_withinsegtrials{roi},3);
    allsubjs_mean_pattern_dist_betweensegtrials{roi} = nanmean(allsubjs_pattern_dist_betweensegtrials{roi},3);
end
% MDS
curr_mat = allsubjs_mean_pattern_dist_withinsegtrials;
roi = 1;

curr_mat{roi}(find(eye(size(curr_mat{roi})))) = 0;      % Zeroing the diagonal
curr_mat{roi} = curr_mat{roi} - min(curr_mat{roi}(:));     % Normalizing so all values will be positive
curr_mat{roi}(find(eye(size(curr_mat{roi})))) = 0;      % Zeroing the diagonal
curr_mat{roi} = curr_mat{roi} / max(curr_mat{roi}(:));  % Normalizing to range 0-1

Y = cmdscale(curr_mat{roi},2);      % Multidimensional scaling

% Plotting the MDS, with gradient descent to fit it to actual distance locations as much as possible
nums={}; for i=1:16, nums{i}=num2str(i);end
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
l1 = locations_all_objects;
l1 = l1/150;
i = -90; R = [cosd(i) -sind(i); sind(i) cosd(i)]; l1 = l1*R; % Rotating to correct environment orientation
% Implementing a numerical gradient descent
param = {[0:360], [-0.02:0.01:0.02], [-0.02:0.01:0.02], [1,1.1], [-1 1]};    % Rotation degree, x-translation, y-translation, scaling, flipping
max_iters = 100;
max_dist_improvement = 0.0001;
curr_dist_improvement = 100;
iter = 1;
while iter <= max_iters && curr_dist_improvement > max_dist_improvement
    % Calculating the current distance
    curr_dist = pdist2(l1, Y);
    curr_dist = sum(diag(curr_dist).^2);    % sum of squared distances
    % Going over parameters to find the best ones
    all_dists = [];
    for deg = 1:length(param{1})     % Degree of rotation
        for xtrans = 1:length(param{2})    % X translation
            for ytrans = 1:length(param{3})  % Y translation
                for sc = 1:length(param{4})    % Scaling
                    for rot = 1:length(param{5})    % Flipping
                        curr_deg = param{1}(deg); curr_xtrans = param{2}(xtrans); curr_ytrans = param{3}(ytrans); curr_sc = param{4}(sc); curr_rot = param{5}(rot);
                        Y_new = Y;
                        Y_new(:,1) = curr_rot * Y_new(:,1);
                        R = [cosd(curr_deg) -sind(curr_deg); sind(curr_deg) cosd(curr_deg)];
                        Y_new = Y_new*R;
                        Y_new(:,1) = Y_new(:,1) + curr_xtrans;
                        Y_new(:,2) = Y_new(:,2) + curr_ytrans;
                        Y_new = Y_new * curr_sc;
                        dist_temp = pdist2(l1, Y_new);
                        all_dists(deg,xtrans,ytrans,sc,rot) = sum(diag(dist_temp).^2);
                    end
                end
            end
        end
    end
    [D,I] = min(all_dists(:));
    [I1,I2,I3,I4,I5] = ind2sub(size(all_dists),I);
    curr_deg = param{1}(I1); curr_xtrans = param{2}(I2); curr_ytrans = param{3}(I3); curr_sc = param{4}(I4); curr_rot = param{5}(I5);
    Y_new = Y;
    Y_new(:,1) = curr_rot * Y_new(:,1);
    R = [cosd(curr_deg) -sind(curr_deg); sind(curr_deg) cosd(curr_deg)];
    Y_new = Y_new*R;
    Y_new(:,1) = Y_new(:,1) + curr_xtrans;
    Y_new(:,2) = Y_new(:,2) + curr_ytrans;
    Y_new = Y_new * curr_sc;
    Y = Y_new;
    curr_dist_improvement = curr_dist - D;
    disp([iter, curr_dist_improvement])
    iter = iter + 1;
end
% figure;plot(Y(:,1),Y(:,2),'*'), text(Y(:,1),Y(:,2),nums), title(['roi ' all_ROI_names{roi}]);
% hold on; scatter(l1(:,1),l1(:,2),'r'), text(l1(:,1),l1(:,2),nums,'r')
% diffs = Y-l1; quiver(l1(:,1),l1(:,2),diffs(:,1),diffs(:,2),0)
% figure;hold on;plot(Y(1:4,1),Y(1:4,2),'.g',Y(5:8,1),Y(5:8,2),'.b',Y(9:12,1),Y(9:12,2),'.r',Y(13:16,1),Y(13:16,2),'.y','MarkerSize',30), text(Y(:,1),Y(:,2),nums,'FontSize',12), title(['roi ' all_ROI_names{roi}]);
figure;hold on;
scatter(Y(1:2,1),Y(1:2,2),200,'g','p','filled')
scatter(Y(3:4,1),Y(3:4,2),100,'g','s','filled')
scatter(Y([5,8],1),Y([5,8],2),200,'b','p','filled')
scatter(Y([6,7],1),Y([6,7],2),100,'b','s','filled')
scatter(Y([11,12],1),Y([11,12],2),200,'r','p','filled')
scatter(Y([9,10],1),Y([9,10],2),100,'r','s','filled')
scatter(Y([14,15],1),Y([14,15],2),200,'y','p','filled')
scatter(Y([13,16],1),Y([13,16],2),100,'y','s','filled')
text(Y(:,1)+0.01,Y(:,2),nums,'FontSize',12), title(['roi ' all_ROI_names{roi}])
