% Code to analyze the JRD data using single-trial beta estimates instead of an average pattern across the run

betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_EC'};
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

% Schema preserving (termed "ROTATION" in the paper) - with distance
locations_new_schema = locations_all_objects;
for i = 1:4
    locations_new_schema(i, :) = locations_new_schema(i + num_conditions/2, :);
end
for i = 13:16
    locations_new_schema(i, :) = locations_new_schema(i - num_conditions/2, :);
end
mat_schema = pdist2(locations_new_schema, locations_new_schema);         % Distances between all objects
mat_schema = mat_schema / max(mat_schema(:));    % Normalizing to range 0-1
mat_schema = 1 - mat_schema;

% Flipping (termed "MIRRORING" in the paper) - with distance
locations_new_flipping = locations_all_objects;
locations_new_flipping(:,2) = abs(locations_new_flipping(:,2));
mat_flipping = pdist2(locations_new_flipping, locations_new_flipping);         % Distances between all objects
mat_flipping = mat_flipping / max(mat_flipping(:));    % Normalizing to range 0-1
mat_flipping = 1 - mat_flipping;

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


% Quadrants overlay model (while maintaining the original N/S/E/W directions)
locations_new_quadrants_overlay = locations_all_objects;
locations_new_quadrants_overlay(5:12,2) = locations_new_quadrants_overlay(5:12,2) + 75;
locations_new_quadrants_overlay(9:16,1) = locations_new_quadrants_overlay(9:16,1) + 75;
mat_quadrants_overlay = pdist2(locations_new_quadrants_overlay, locations_new_quadrants_overlay);         % Distances between all objects
mat_quadrants_overlay = mat_quadrants_overlay / max(mat_quadrants_overlay(:));    % Normalizing to range 0-1
mat_quadrants_overlay = 1 - mat_quadrants_overlay;
% Quadrants overlay model (while maintaining corner and river locations)
locations_new_quadrants_overlay_flipping = locations_all_objects;
locations_new_quadrants_overlay_flipping(:,2) = abs(locations_new_quadrants_overlay_flipping(:,2));
locations_new_quadrants_overlay_flipping(:,1) = abs(locations_new_quadrants_overlay_flipping(:,1));
mat_quadrants_overlay_flipping = pdist2(locations_new_quadrants_overlay_flipping, locations_new_quadrants_overlay_flipping);         % Distances between all objects
mat_quadrants_overlay_flipping = mat_quadrants_overlay_flipping / max(mat_quadrants_overlay_flipping(:));    % Normalizing to range 0-1
mat_quadrants_overlay_flipping = 1 - mat_quadrants_overlay_flipping;



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
    
    
    if num_runs_jrd>0
        %% Getting the information on objects and directions in each trial and creating dissimilarity matrices
        all_curr_subj_info = [];
        all_object1 = [];
        all_object3 = [];
        all_viewing_direction = [];
        all_quadrant1 = [];
        all_quadrant2 = [];
        all_obj13_dist=[];
        for r = 1:num_runs_jrd
            curr_info_file = dir(fullfile(curr_subj_jrd_analysis_dirs{r}, 'all_info*.mat'));
            curr_info = load(fullfile(curr_info_file.folder, curr_info_file.name));
            all_curr_subj_info = [all_curr_subj_info;curr_info.all_info];
            all_object1 = [all_object1; all_curr_subj_info(r).object1_unity_location(1:48) + 1];
            all_object3 = [all_object3; all_curr_subj_info(r).object3_unity_location(1:48) + 1];
            all_viewing_direction = [all_viewing_direction; all_curr_subj_info(r).viewing_direction(1:48)'];
            all_quadrant1 = [all_quadrant1; all_curr_subj_info(r).quadrant_start(1:48)];
            all_quadrant2 = [all_quadrant2; all_curr_subj_info(r).quadrant_target(1:48)];
        end
        for i=1:length(all_object1)
            all_obj13_dist(i) = mat_distances(all_object1(i), all_object3(i));
        end
        % Figuring out the segment for each trial
        all_segment1 = (all_quadrant1 == 2 | all_quadrant1 == 3);
        all_segment2 = (all_quadrant2 == 2 | all_quadrant2 == 3);
        all_orth_segment1 = (all_quadrant1 == 3 | all_quadrant1 == 4);
        all_orth_segment2 = (all_quadrant2 == 3 | all_quadrant2 == 4);
        % Viewing direction of each trial - lowering the resolution to bins
        all_viewing_direction_bins = nan(size(all_viewing_direction));
        all_viewing_direction_bins(all_viewing_direction>-45 & all_viewing_direction<=45) = 1;      % North - toward mountains
        all_viewing_direction_bins(all_viewing_direction>45 & all_viewing_direction<=135) = 2;      % West
        all_viewing_direction_bins(all_viewing_direction>135 | all_viewing_direction<=-135) = 3;    % South
        all_viewing_direction_bins(all_viewing_direction>-135 & all_viewing_direction<=-45) = 4;    % East
        
        
        % Calculating matrices for all trials differences in distance, viewing direction, etc.
        curr_distances_all_object1 = nan(length(all_object1));
        curr_distances_all_object1_overlay = nan(length(all_object1));
        curr_distances_all_object1_grouping = nan(length(all_object1));
        curr_distances_all_object1_within_segment = nan(length(all_object1));
        curr_distances_all_object1_between_segments = nan(length(all_object1));
        curr_same_object1 = nan(length(all_object1));
        curr_same_segment = nan(length(all_object1));
        curr_same_orth_segment = nan(length(all_object1));
        curr_same_quadrant = nan(length(all_object1));
        curr_viewdir_diff_all_objects = nan(length(all_object1));
        curr_same_viewdir = nan(length(all_object1));
        curr_same_viewdir_bin = nan(length(all_object1));
        
        for i = 1:length(all_object1)
            for j = 1:length(all_object1)
                if i ~= j
                    % Object distance coding and segment / quadrant coding
                    %                 if all_object1(i) ~= all_object1(j)     % If this is not the same starting object
                    % Distance between objects in different trials
                    curr_distances_all_object1(i,j) = mat_distances(all_object1(i), all_object1(j));
                    curr_distances_all_object1_overlay(i,j) = mat_overlay(all_object1(i), all_object1(j));
                    curr_distances_all_object1_grouping(i,j) = mat_grouping(all_object1(i), all_object1(j));
                    curr_distances_all_object1_within_segment(i,j) = mat_dist_within_segment(all_object1(i), all_object1(j));
                    curr_distances_all_object1_between_segments(i,j) = mat_dist_between_segments(all_object1(i), all_object1(j));
                    % Same / different segment or quadrant between trials
                    curr_same_segment(i,j) = (all_segment1(i) == all_segment1(j));
                    curr_same_orth_segment(i,j) = (all_orth_segment1(i) == all_orth_segment1(j));
                    curr_same_quadrant(i,j) = (all_quadrant1(i) == all_quadrant1(j));
                    curr_same_object1(i,j) = 0;
                    %                     else
                    %                         curr_same_object1(i,j) = 1;
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
            %             end
        end
        
        % Converting the viewing direction matrix to similarity matrix (similar angles -> high value, dissimilar angles - low values)
        curr_viewdir_diff_all_objects = curr_viewdir_diff_all_objects / max(curr_viewdir_diff_all_objects(:));
        curr_viewdir_diff_all_objects = 1 - curr_viewdir_diff_all_objects;
        
        curr_diff_obj13_dist = [];
        for i=1:length(all_object1)
            for j=1:length(all_object1)
                curr_diff_obj13_dist(i,j) = abs(all_obj13_dist(i) - all_obj13_dist(j));
            end
        end
        
        % Getting the non-diagonal elements of each matrix
        curr_distances_all_object1_nondiag = curr_distances_all_object1(find(tril(ones(length(all_object1)),-1)));
        curr_diff_obj13_dist_nondiag = curr_diff_obj13_dist(find(tril(ones(length(all_object1)),-1)));
        curr_viewdir_diff_nondiag = curr_viewdir_diff_all_objects(find(tril(ones(length(all_object1)),-1)));
        
        corr_diff_dists(subj) = corr(curr_diff_obj13_dist_nondiag, curr_distances_all_object1_nondiag,'rows','complete');
        corr_dist_viewdir(subj) = corr(curr_viewdir_diff_nondiag, curr_distances_all_object1_nondiag,'rows','complete');
                
    end
end


%% Creating an average viewing direction per object matrix
average_viewdir_per_location = nan(num_conditions,1);
for i=1:num_conditions
    average_viewdir_per_location(i) = mean(all_viewing_direction(all_object1 == i));
end
mat_diff_average_viewdir = [];
for i=1:num_conditions
    for j=1:num_conditions
        curr_angle_diff = average_viewdir_per_location(i) - average_viewdir_per_location(j);
        if abs(curr_angle_diff) > 180, curr_angle_diff = curr_angle_diff - 360*sign(curr_angle_diff); end
        if curr_angle_diff < 0, curr_angle_diff = 360 + curr_angle_diff; end
        if curr_angle_diff > 180, curr_angle_diff = 360 - curr_angle_diff; end     % Choosing the smaller of the two possible angles
        mat_diff_average_viewdir(i,j) = curr_angle_diff;
    end
end
mat_diff_average_viewdir = mat_diff_average_viewdir / max(mat_diff_average_viewdir(:));
mat_diff_average_viewdir = 1-mat_diff_average_viewdir;
