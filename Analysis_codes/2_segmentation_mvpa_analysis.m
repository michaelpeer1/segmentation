% Defining general variables
betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview','bilat_EC'};
%     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC','rHC','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview','lEC','rEC'};
% all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview'};
%     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC','rHC','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview'};
% %     all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC_activation_jrd','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview', 'bilat_HC_activation_objview'};
% %     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC_activation_jrd','rHC_activation_jrd','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview', 'lHC_activation_objview','rHC_activation_objview'};
%     all_ROI_names = {'bilat_HC_activation_jrd', 'bilat_HC_activation_objview', 'bilat_RSC', 'bilat_PPA', 'bilat_OPA', 'bilat_LOC'};
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

% Locations of diagonal quadrants
mat_diagonal_quadrants = zeros(num_conditions); mat_diagonal_quadrants(1:4,9:12) = 1; mat_diagonal_quadrants(9:12, 1:4) = 1; mat_diagonal_quadrants(5:8,13:16) = 1; mat_diagonal_quadrants(13:16,5:8) = 1;

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
% Schema preserving - with distance between items in different segments only
mat_schema_between_segments = mat_schema;
mat_schema_between_segments(mat_segments==1) = nan;
% Flipping - between segments adjacent quadrants only
mat_schema_between_segments_adj_quadrants = mat_schema;
mat_schema_between_segments_adj_quadrants(mat_segments==1) = nan;
mat_schema_between_segments_adj_quadrants(mat_diagonal_quadrants==1) = nan;
% Schema preserving - y axis only model
mat_schema_y_axis_only = pdist2(locations_new_schema(:,1), locations_new_schema(:,1));
mat_schema_y_axis_only = mat_schema_y_axis_only / max(mat_schema_y_axis_only(:));    % Normalizing to range 0-1
mat_schema_y_axis_only = 1 - mat_schema_y_axis_only;



% Flipping (termed "MIRRORING" in the paper) - with distance
locations_new_flipping = locations_all_objects;
locations_new_flipping(:,2) = abs(locations_new_flipping(:,2));
mat_flipping = pdist2(locations_new_flipping, locations_new_flipping);         % Distances between all objects
mat_flipping = mat_flipping / max(mat_flipping(:));    % Normalizing to range 0-1
mat_flipping = 1 - mat_flipping;
% Flipping - with distance between items in different segments only
mat_flipping_between_segments = mat_flipping;
mat_flipping_between_segments(mat_segments==1) = nan;
% Flipping - between segments adjacent quadrants only
mat_flipping_between_segments_adj_quadrants = mat_flipping;
mat_flipping_between_segments_adj_quadrants(mat_segments==1) = nan;
mat_flipping_between_segments_adj_quadrants(mat_diagonal_quadrants==1) = nan;
% Flipping - x axis only model
mat_flipping_x_axis_only = pdist2(locations_new_flipping(:,2), locations_new_flipping(:,2));
mat_flipping_x_axis_only = mat_flipping_x_axis_only / max(mat_flipping_x_axis_only(:));    % Normalizing to range 0-1
mat_flipping_x_axis_only = 1 - mat_flipping_x_axis_only;
% Creating the opposite matrix - flipping orthogonal to the river
locations_new_flipping_orthogonal = locations_all_objects;
locations_new_flipping_orthogonal(:,1) = abs(locations_new_flipping_orthogonal(:,1));
mat_orth_flipping = pdist2(locations_new_flipping_orthogonal, locations_new_flipping_orthogonal);         % Distances between all objects
mat_orth_flipping = mat_orth_flipping / max(mat_orth_flipping(:));    % Normalizing to range 0-1
mat_orth_flipping = 1 - mat_orth_flipping;

% Overlay - with distance
locations_new_overlay = locations_all_objects;
locations_new_overlay(5:12,2) = locations_new_overlay(5:12,2) + 75;
mat_overlay = pdist2(locations_new_overlay, locations_new_overlay);         % Distances between all objects
mat_overlay = mat_overlay / max(mat_overlay(:));    % Normalizing to range 0-1
mat_overlay = 1 - mat_overlay;
% Overlay - with distance between items in different segments only
mat_overlay_between_segments = mat_overlay;
mat_overlay_between_segments(mat_segments==1) = nan;
% Overlay - between segments adjacent quadrants only
mat_overlay_between_segments_adj_quadrants = mat_overlay;
mat_overlay_between_segments_adj_quadrants(mat_segments==1) = nan;
mat_overlay_between_segments_adj_quadrants(mat_diagonal_quadrants==1) = nan;
% Overlay - x axis only model
mat_overlay_x_axis_only = pdist2(locations_new_overlay(:,2), locations_new_overlay(:,2));
mat_overlay_x_axis_only = mat_overlay_x_axis_only / max(mat_overlay_x_axis_only(:));    % Normalizing to range 0-1
mat_overlay_x_axis_only = 1 - mat_overlay_x_axis_only;
% Creating the opposite matrix - overlay orthogonal to the river
locations_new_overlay_orthogonal = locations_all_objects;
locations_new_overlay_orthogonal(9:16,1) = locations_new_overlay_orthogonal(9:16,1) + 75;
mat_orth_overlay = pdist2(locations_new_overlay_orthogonal, locations_new_overlay_orthogonal);         % Distances between all objects
mat_orth_overlay = mat_orth_overlay / max(mat_orth_overlay(:));    % Normalizing to range 0-1
mat_orth_overlay = 1 - mat_orth_overlay;

% Organization on X axis (between segments direction)
mat_x_axis = pdist2(locations_all_objects(:,2), locations_all_objects(:,2));
mat_x_axis = mat_x_axis / max(mat_x_axis(:));    % Normalizing to range 0-1
mat_x_axis = 1 - mat_x_axis;
% Organization on Y axis (along river direction)
mat_y_axis = pdist2(locations_all_objects(:,1), locations_all_objects(:,1));
mat_y_axis = mat_y_axis / max(mat_y_axis(:));    % Normalizing to range 0-1
mat_y_axis = 1 - mat_y_axis;
% Organization on X and Y axis within segment and within orthogonal segment only
mat_x_axis_within_segment = mat_x_axis;
mat_x_axis_within_segment(mat_segments==0) = nan;
mat_x_axis_between_segments = mat_x_axis;
mat_x_axis_between_segments(mat_segments==1) = nan;
mat_y_axis_within_segment = mat_y_axis;
mat_y_axis_within_segment(mat_segments==0) = nan;
mat_y_axis_between_segments = mat_y_axis;
mat_y_axis_between_segments(mat_segments==1) = nan;

% Distance within segment only
mat_dist_within_segment = mat_distances;
mat_dist_within_segment(mat_segments==0) = nan;
% Distance between segments only
mat_dist_between_segments = mat_distances;
mat_dist_between_segments(mat_segments==1) = nan;
% Distance within segment, adjacent quadrants only
mat_dist_within_segment_adj_quadrants = mat_distances;
mat_dist_within_segment_adj_quadrants(mat_segments==0) = nan;
mat_dist_within_segment_adj_quadrants(mat_quadrants==1) = nan;
% Distance between segments, adjacent quadrants only
mat_dist_between_segments_adj_quadrants = mat_distances;
mat_dist_between_segments_adj_quadrants(mat_segments==1) = nan;
mat_dist_between_segments_adj_quadrants(mat_diagonal_quadrants==1) = nan;
% Distance between objects in diagonal quadrants only
mat_dist_diagonal_quadrants = nan(num_conditions);
mat_dist_diagonal_quadrants(mat_diagonal_quadrants==1) = mat_distances(mat_diagonal_quadrants==1);
% Distances within quadrant only
mat_dist_within_quadrant = mat_distances;
mat_dist_within_quadrant(mat_quadrants==0) = nan;

% Distances between quadrants (generalizing within quadrant)
mat_quadrants_dist = mat_quadrants;
mat_quadrants_dist(1:4,[5:8,13:16]) = 0.5; mat_quadrants_dist([5:8,13:16],1:4) = 0.5;
mat_quadrants_dist(5:8,9:12) = 0.5; mat_quadrants_dist(9:12,5:8) = 0.5;
mat_quadrants_dist(9:12,13:16) = 0.5; mat_quadrants_dist(13:16,9:12) = 0.5;

% Distance from river
distances_from_river = abs(locations_all_objects(:,2));
mat_distance_from_river = pdist2(distances_from_river, distances_from_river);
mat_distance_from_river = mat_distance_from_river / max(mat_distance_from_river(:));    % Normalizing to range 0-1
mat_distance_from_river = 1 - mat_distance_from_river;

% A matrix coding if it is an adjacent quadrant within and between segment
mat_within_seg_adj_quad = zeros(num_conditions);
mat_within_seg_adj_quad(1:4,13:16) = 1; mat_within_seg_adj_quad(5:8,9:12) = 1; mat_within_seg_adj_quad(9:12,5:8) = 1; mat_within_seg_adj_quad(13:16,1:4) = 1;
mat_between_seg_adj_quad = zeros(num_conditions);
mat_between_seg_adj_quad(1:4,5:8) = 1; mat_between_seg_adj_quad(5:8,1:4) = 1; mat_between_seg_adj_quad(9:12,13:16) = 1; mat_between_seg_adj_quad(13:16,9:12) = 1;

% A matrix coding a segment grouping index - pattern similarity within segment larger than between segments (adjacent quadrants only)
mat_segment_grouping = nan(num_conditions);
mat_segment_grouping(mat_within_seg_adj_quad == 1) = 1;
mat_segment_grouping(mat_between_seg_adj_quad == 1) = -1;
% % A matrix coding a segment grouping index - all regions - by calculating the expected difference under a Euclidean assumption
% grouping_newmeasure_factor = nanmean(mat_distances(mat_segments==1 & eye(num_conditions)~=1)) / nanmean(mat_distances(mat_segments==0 & eye(num_conditions)~=1));
% mat_segment_grouping_newmeasure = nan(num_conditions);
% mat_segment_grouping_newmeasure(mat_segments == 1) = 1;
% mat_segment_grouping_newmeasure(mat_segments == 0) = -1 * grouping_newmeasure_factor;


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


% Lower triangle elements only
mat_nondiag_distance = mat_distances(find(tril(ones(num_conditions),-1)));
mat_nondiag_segments = mat_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_orthsegments = mat_orth_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_quadrants = mat_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_schema = mat_schema(find(tril(ones(num_conditions),-1)));
mat_nondiag_schema_btwn_segments = mat_schema_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_schema_between_seg_adj_quad = mat_schema_between_segments_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_schema_yaxis_only = mat_schema_y_axis_only(find(tril(ones(num_conditions),-1)));
mat_nondiag_flipping = mat_flipping(find(tril(ones(num_conditions),-1)));
mat_nondiag_flipping_btwn_segments = mat_flipping_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_flipping_between_seg_adj_quad = mat_flipping_between_segments_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_flipping_xaxis_only = mat_flipping_x_axis_only(find(tril(ones(num_conditions),-1)));
mat_nondiag_orth_flipping = mat_orth_flipping(find(tril(ones(num_conditions),-1)));
mat_nondiag_overlay = mat_overlay(find(tril(ones(num_conditions),-1)));
mat_nondiag_overlay_btwn_segments = mat_overlay_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_overlay_between_seg_adj_quad = mat_overlay_between_segments_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_overlay_xaxis_only = mat_overlay_x_axis_only(find(tril(ones(num_conditions),-1)));
mat_nondiag_orth_overlay = mat_orth_overlay(find(tril(ones(num_conditions),-1)));
mat_nondiag_y_axis = mat_y_axis(find(tril(ones(num_conditions),-1)));
mat_nondiag_y_axis_between = mat_y_axis_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_y_axis_within = mat_y_axis_within_segment(find(tril(ones(num_conditions),-1)));
mat_nondiag_x_axis = mat_x_axis(find(tril(ones(num_conditions),-1)));
mat_nondiag_x_axis_between = mat_x_axis_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_x_axis_within = mat_x_axis_within_segment(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_within_segment = mat_dist_within_segment(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_between_segments = mat_dist_between_segments(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_within_segment_adj_quadrants = mat_dist_within_segment_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_between_segments_adj_quadrants = mat_dist_between_segments_adj_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_diagonal_quadrants = mat_dist_diagonal_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_within_quadrant = mat_dist_within_quadrant(find(tril(ones(num_conditions),-1)));
mat_nondiag_dist_diagonal_quadrants = mat_dist_diagonal_quadrants(find(tril(ones(num_conditions),-1)));
mat_nondiag_distance_from_river = mat_distance_from_river(find(tril(ones(num_conditions),-1)));
mat_nondiag_quadrants_dist = mat_quadrants_dist(find(tril(ones(num_conditions),-1)));
mat_nondiag_within_seg_adj_quad = mat_within_seg_adj_quad(find(tril(ones(num_conditions),-1)));
mat_nondiag_between_seg_adj_quad = mat_between_seg_adj_quad(find(tril(ones(num_conditions),-1)));
mat_nondiag_segment_grouping = mat_segment_grouping(find(tril(ones(num_conditions),-1)));
% mat_nondiag_segment_grouping_newmeasure = mat_segment_grouping_newmeasure(find(tril(ones(num_conditions),-1)));
mat_nondiag_quadrants_overlay = mat_quadrants_overlay(find(tril(ones(num_conditions),-1)));
mat_nondiag_quadrants_overlay_flipping = mat_quadrants_overlay_flipping(find(tril(ones(num_conditions),-1)));

% Creating a template that balances the diagonal minus nondiagonal elements, to calculate their weighted subraction
diag_nondiag_template = (eye(num_conditions) - 1 / num_conditions) / (num_conditions - 1);



%% RSA in ROIs

% Initialize variables
means_diag_minus_nondiag_objview_halfruns = nan(length(subj_data_dirs), num_ROIs);

corr_objview_d2minusd1_distances = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_orth_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_schema = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_schema_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_schema_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_schema_yaxis_only = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_flipping = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_flipping_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_flipping_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_flipping_xaxis_only = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_orth_flipping = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_overlay_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_overlay_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_overlay_xaxis_only = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_orth_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_distance_from_river = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_x_axis = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_y_axis = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_dist_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_dist_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_dist_within_quadrant = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_dist_within_segment_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_dist_between_segments_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_dist_diagonal_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_quadrants_dist = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_estimated_dist = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_x_axis_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_y_axis_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_x_axis_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_y_axis_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_quadrants_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_objview_d2minusd1_quadrants_overlay_flipping = nan(length(subj_data_dirs), num_ROIs);

corr_objview_day2_distances = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_orth_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_coviewing = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_schema = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_schema_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_schema_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_schema_y_axis_only = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_flipping = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_flipping_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_flipping_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_flipping_x_axis_only = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_orth_flipping = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_overlay_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_overlay_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_overlay_x_axis_only = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_orth_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_distance_from_river = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_x_axis = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_y_axis = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_dist_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_dist_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_dist_within_quadrant = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_dist_within_segment_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_dist_between_segments_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_dist_diagonal_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_quadrants_dist = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_estimated_dist = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_x_axis_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_y_axis_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_x_axis_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_y_axis_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_quadrants_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_objview_day2_quadrants_overlay_flipping = nan(length(subj_data_dirs), num_ROIs);

corr_jrd_distances = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_segments = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_orth_segments = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_schema = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_schema_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_schema_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_schema_y_axis_only = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_flipping = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_flipping_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_flipping_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_flipping_x_axis_only = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_orth_flipping = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_overlay_btwn = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_overlay_btwn_seg_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_overlay_x_axis_only = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_orth_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_distance_from_river = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_x_axis = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_y_axis = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_dist_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_dist_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_dist_within_quadrant = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_dist_within_segment_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_dist_between_segments_adj_quad = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_dist_diagonal_quadrants = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_quadrants_dist = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_estimated_dist = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_x_axis_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_y_axis_within_segment = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_x_axis_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_y_axis_between_segments = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_quadrants_overlay = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_quadrants_overlay_flipping = nan(length(subj_data_dirs), num_ROIs);

corr_jrd_coviewing = nan(length(subj_data_dirs), num_ROIs);
corr_jrd_resp_similarity = nan(length(subj_data_dirs), num_ROIs);

partialcorr_objview_d2minusd1_distances = nan(length(subj_data_dirs), num_ROIs);
partialcorr_objview_d2minusd1_segments = nan(length(subj_data_dirs), num_ROIs);
partialcorr_objview_d2minusd1_orthsegments = nan(length(subj_data_dirs), num_ROIs);
partialcorr_objview_d2minusd1_quadrants = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_distances = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_distances_cov_resp = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_overlay_cov_resp = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_flipping_cov_resp = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_yaxis_cov_resp = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_segments = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_orthsegments = nan(length(subj_data_dirs), num_ROIs);
partialcorr_jrd_quadrants = nan(length(subj_data_dirs), num_ROIs);

segment_grouping_objview_d2minusd1 = nan(length(subj_data_dirs), num_ROIs);
segment_grouping_newmeasure_objview_d2minusd1 = nan(length(subj_data_dirs), num_ROIs);
segment_remapping_objview_d2minusd1_full = nan(length(subj_data_dirs), num_ROIs);
segment_remapping_objview_d2minusd1_adjquads = nan(length(subj_data_dirs), num_ROIs);
segment_grouping_jrd = nan(length(subj_data_dirs), num_ROIs);
segment_grouping_newmeasure_jrd = nan(length(subj_data_dirs), num_ROIs);
segment_remapping_jrd_full = nan(length(subj_data_dirs), num_ROIs);
segment_remapping_jrd_adjquads = nan(length(subj_data_dirs), num_ROIs);

all_neural_corr_ROI_jrd_nondiag = nan(length(mat_nondiag_distance), num_ROIs, length(subj_data_dirs));
all_neural_corr_objview_day2_nondiag = nan(length(mat_nondiag_distance), num_ROIs, length(subj_data_dirs));
all_neural_corr_objview_d2minusd1_nondiag = nan(length(mat_nondiag_distance), num_ROIs, length(subj_data_dirs));

% Calculation RSA in each subject and ROI
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
    
    
    %% Loading the subject specific ROIs
    all_ROIs_individual = {};
    for i = 1:length(all_ROI_names)
        all_ROIs_individual{i} = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
    end
    
    
    %% Reading the individual subject's estimated distance matrix
    mat_estimated_dist = load(fullfile(curr_subj_jrd_analysis_dir, 'estimated_distances_mat.mat'));
    mat_estimated_dist = mat_estimated_dist.dist_matrix;
    % Converting to a maximum of 1
    mat_estimated_dist = mat_estimated_dist / max(mat_estimated_dist(:));
    % Converting to similarity matrix
    mat_estimated_dist = 1 - mat_estimated_dist;
    mat_nondiag_estimated_dist = mat_estimated_dist(find(tril(ones(num_conditions),-1)));
    
    
    %% MVPA analysis in each ROI
    for roi = 1:num_ROIs        % Going over ROIs
        curr_ROI = find(all_ROIs_individual{roi}(:)==1);
        if ~isempty(curr_ROI)
            % Reading object viewing data in the ROI
            all_tvalues_ROI_objview = nan(length(curr_ROI), num_conditions, length(curr_subj_objectviewing_combined_dirs));
            for r = 1:length(curr_subj_objectviewing_combined_dirs)
                t_values_files_current = dir(fullfile(curr_subj_objectviewing_combined_dirs{r}, 'spmT*.nii'));
                for b = 1:num_conditions
                    curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_objectviewing_combined_dirs{r}, t_values_files_current(b).name)));
                    all_tvalues_ROI_objview(:,b,r) = curr_tvalues_image(curr_ROI);
                end
            end
            % Removing the cocktail mean from each day (mean activity at each voxel across patterns)
            for r = 1:length(curr_subj_objectviewing_combined_dirs)
                mean_day_pattern = nanmean(all_tvalues_ROI_objview(:,:,r), 2);
                for b = 1:num_conditions
                    all_tvalues_ROI_objview(:,b,r) = all_tvalues_ROI_objview(:,b,r) - mean_day_pattern;
                end
            end
            
            % Object viewing runs - calculating correlation between day1 patterns and day2 patterns
            curr_corr = fisherz(corr(all_tvalues_ROI_objview(:,:,1), all_tvalues_ROI_objview(:,:,2), 'rows', 'complete'));                 % Correlation between first and second halves of runs
            means_diag_minus_nondiag_objview_halfruns(subj, roi) = sum(sum(diag_nondiag_template .* curr_corr));                          % Diagonal minus nondiagonal elements
            
            % Object viewing runs - calculating correlation to distance and segment matrices
            corr_ROI_objview_day1 = fisherz(corr(all_tvalues_ROI_objview(:,:,1), all_tvalues_ROI_objview(:,:,1), 'rows', 'complete'));                    % Correlation between day 1 patterns, averaged across runs
            corr_ROI_objview_day2 = fisherz(corr(all_tvalues_ROI_objview(:,:,2), all_tvalues_ROI_objview(:,:,2),'rows', 'complete'));                    % Correlation between day 2 patterns, averaged across runs
            % Non-diagonal elements
            corr_ROI_objview_day1_nondiag = corr_ROI_objview_day1(find(tril(ones(num_conditions),-1)));       % Taking only the lower triangle elements as a vector
            corr_ROI_objview_day2_nondiag = corr_ROI_objview_day2(find(tril(ones(num_conditions),-1)));       % Taking only the lower triangle elements as a vector
            corr_ROI_objview_d2minusd1_nondiag = (corr_ROI_objview_day2_nondiag - corr_ROI_objview_day1_nondiag);            % Day2 minus day1
            
            % Calculating correlation to different matrices
            % Day 2
            corr_objview_day2_distances(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_segments(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_quadrants(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_orth_segments(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_orthsegments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_schema(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_schema, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_schema_btwn(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_schema_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_schema_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_schema_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_schema_y_axis_only(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_schema_yaxis_only, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_flipping(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_flipping, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_flipping_btwn(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_flipping_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_flipping_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_flipping_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_flipping_x_axis_only(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_flipping_xaxis_only, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_orth_flipping(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_orth_flipping, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_overlay(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_overlay_btwn(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_overlay_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_overlay_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_overlay_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_overlay_x_axis_only(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_overlay_xaxis_only, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_orth_overlay(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_orth_overlay, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_dist_within_segment(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_dist_between_segments(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_dist_within_quadrant(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_dist_within_quadrant, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_dist_within_segment_adj_quad(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_dist_within_segment_adj_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_dist_between_segments_adj_quad(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_dist_between_segments_adj_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_dist_diagonal_quadrants(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_dist_diagonal_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_x_axis(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_x_axis, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_x_axis_between_segments(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_x_axis_between, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_x_axis_within_segment(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_x_axis_within, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_y_axis(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_y_axis, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_y_axis_between_segments(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_y_axis_between, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_y_axis_within_segment(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_y_axis_within, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_distance_from_river(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_distance_from_river, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_quadrants_dist(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_quadrants_dist, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_estimated_dist(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_estimated_dist, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_quadrants_overlay(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_quadrants_overlay, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_day2_quadrants_overlay_flipping(subj, roi) = corr(corr_ROI_objview_day2_nondiag, mat_nondiag_quadrants_overlay_flipping, 'Type', 'Spearman', 'rows', 'complete');
            
            % Day 2 minus day 1 (change in corrrelation between days)
            corr_objview_d2minusd1_distances(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_segments(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_quadrants(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_orth_segments(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_orthsegments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_schema(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_schema, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_schema_btwn(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_schema_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_schema_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_schema_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_schema_y_axis_only(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_schema_yaxis_only, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_flipping(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_flipping, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_flipping_btwn(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_flipping_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_flipping_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_flipping_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_flipping_x_axis_only(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_flipping_xaxis_only, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_orth_flipping(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_orth_flipping, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_overlay(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_overlay_btwn(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_overlay_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_overlay_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_overlay_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_overlay_x_axis_only(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_overlay_xaxis_only, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_orth_overlay(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_orth_overlay, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_dist_within_segment(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_dist_between_segments(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_dist_within_quadrant(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_within_quadrant, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_dist_within_segment_adj_quad(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_within_segment_adj_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_dist_between_segments_adj_quad(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_between_segments_adj_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_dist_diagonal_quadrants(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_dist_diagonal_quadrants, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_x_axis(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_x_axis, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_x_axis_between_segments(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_x_axis_between, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_x_axis_within_segment(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_x_axis_within, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_y_axis(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_y_axis, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_y_axis_between_segments(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_y_axis_between, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_y_axis_within_segment(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_y_axis_within, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_distance_from_river(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_distance_from_river, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_quadrants_dist(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_quadrants_dist, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_estimated_dist(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_estimated_dist, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_quadrants_overlay(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_quadrants_overlay, 'Type', 'Spearman', 'rows', 'complete');
            corr_objview_d2minusd1_quadrants_overlay_flipping(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_quadrants_overlay_flipping, 'Type', 'Spearman', 'rows', 'complete');
            
            % Object viewing runs - partial correlation - day 2 minus day 1 correlations
            partialcorr_objview_d2minusd1_distances(subj, roi) = partialcorr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_distance, mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
            partialcorr_objview_d2minusd1_segments(subj, roi) = partialcorr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_segments, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            partialcorr_objview_d2minusd1_orthsegments(subj, roi) = partialcorr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_orthsegments, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            partialcorr_objview_d2minusd1_quadrants(subj, roi) = partialcorr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_quadrants, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            partialcorr_objview_d2minusd1_distances_overlay(subj, roi) = partialcorr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_distance, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
            partialcorr_objview_d2minusd1_overlay_distances(subj, roi) = partialcorr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_overlay, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
            
            % Calculating the segment grouping measure (higher neural similarity within segment compared to between)
            segment_grouping_objview_d2minusd1(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
%             segment_grouping_newmeasure_objview_d2minusd1(subj, roi) = corr(corr_ROI_objview_d2minusd1_nondiag, mat_nondiag_segment_grouping_newmeasure, 'Type', 'Spearman', 'rows', 'complete');
            grouping_newmeasure_factor = nanmean(mat_nondiag_distance(mat_nondiag_segments==1))/nanmean(mat_nondiag_distance(mat_nondiag_segments==0));
            segment_grouping_newmeasure_objview_d2minusd1(subj, roi) = grouping_newmeasure_factor * nanmean(corr_ROI_objview_d2minusd1_nondiag(mat_nondiag_segments==1)) - nanmean(corr_ROI_objview_d2minusd1_nondiag(mat_nondiag_segments==0));
            % Calculating the segment remapping measure (higher correlation to distance matrix within segment vs. between segments)
            segment_remapping_objview_d2minusd1_full(subj, roi) = corr_objview_d2minusd1_dist_within_segment(subj, roi) - corr_objview_d2minusd1_dist_between_segments(subj, roi);
            segment_remapping_objview_d2minusd1_adjquads(subj, roi) = corr_objview_d2minusd1_dist_within_segment_adj_quad(subj, roi) - corr_objview_d2minusd1_dist_between_segments_adj_quad(subj, roi);
            
            % Saving the actual neural correlations
            all_neural_corr_objview_day2_nondiag(:,roi,subj) = corr_ROI_objview_day2_nondiag;
            all_neural_corr_objview_d2minusd1_nondiag(:,roi,subj) = corr_ROI_objview_d2minusd1_nondiag;

            
            % JRD
            if num_runs_jrd > 0
                % Extracting the t-values at each ROI - JRD
                all_tvalues_ROI_jrd = nan(length(curr_ROI), num_conditions);
                t_values_files_current = dir(fullfile(curr_subj_jrd_combined_dir, 'spmT*.nii'));
                for b = 1:num_conditions
                    curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_jrd_combined_dir, t_values_files_current(b).name)));
                    all_tvalues_ROI_jrd(:,b) = curr_tvalues_image(curr_ROI);
                end
                % Removing the cocktail mean (mean activity at each voxel across patterns)
                mean_jrd_pattern = nanmean(all_tvalues_ROI_jrd, 2);
                for b = 1:num_conditions
                    all_tvalues_ROI_jrd(:,b) = all_tvalues_ROI_jrd(:,b) - mean_jrd_pattern;
                end
                
                % JRD - Calculating correlation to different matrices
                corr_ROI_jrd = fisherz(corr(all_tvalues_ROI_jrd, all_tvalues_ROI_jrd, 'rows', 'complete'));
                if ~isempty(corr_ROI_jrd)
                    % Getting the non-diagonal elements only of the across-run-combinations-average correlation matrix, averaged across upper and lower triangles
                    corr_ROI_jrd_nondiag = corr_ROI_jrd(find(tril(ones(num_conditions),-1)));
                    
                    % Reading this subject's co-viewing stats and responses similarity stats (created with the behavioral analysis script - analyze_psychopy_all_subjects.m)
                    curr_coviewing_mat_file = fullfile(curr_subj_jrd_analysis_dir, 'object_coviewing_JRD.mat');
                    if exist(curr_coviewing_mat_file, 'file')
                        coviewing_mat = load(fullfile(curr_subj_jrd_analysis_dir, 'object_coviewing_JRD.mat'));
                        coviewing_mat = coviewing_mat.current_subject_coviewing_mat;
                    end
                    curr_resp_similarity_mat_file = fullfile(curr_subj_jrd_analysis_dir, 'responses_similarity_JRD.mat');
                    if exist(curr_resp_similarity_mat_file, 'file')
                        resp_similarity_mat = load(fullfile(curr_subj_jrd_analysis_dir, 'responses_similarity_JRD.mat'));
                        resp_similarity_mat = 1 - resp_similarity_mat.current_responses_diff;       % Converting to similarity instead of dissimilarity
                    end
                    mat_nondiag_coviewing = coviewing_mat(find(tril(ones(num_conditions),-1)));
                    mat_nondiag_resp_similarity = resp_similarity_mat(find(tril(ones(num_conditions),-1)));
                    
                    % Correlations to different matrices
                    corr_jrd_distances(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_segments(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_quadrants(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_quadrants, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_orth_segments(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_orthsegments, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_schema(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_schema, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_schema_btwn(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_schema_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_schema_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_schema_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_schema_y_axis_only(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_schema_yaxis_only, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_flipping(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_flipping, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_flipping_btwn(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_flipping_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_flipping_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_flipping_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_flipping_x_axis_only(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_flipping_xaxis_only, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_orth_flipping(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_orth_flipping, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_overlay(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_overlay_btwn(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_overlay_btwn_segments, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_overlay_btwn_seg_adj_quad(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_overlay_between_seg_adj_quad, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_overlay_x_axis_only(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_overlay_xaxis_only, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_orth_overlay(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_orth_overlay, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_dist_within_segment(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_dist_between_segments(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_dist_within_quadrant(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_dist_within_quadrant, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_dist_within_segment_adj_quad(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_dist_within_segment_adj_quadrants, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_dist_between_segments_adj_quad(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_dist_between_segments_adj_quadrants, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_dist_diagonal_quadrants(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_dist_diagonal_quadrants, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_x_axis(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_x_axis, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_x_axis_between_segments(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_x_axis_between, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_x_axis_within_segment(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_x_axis_within, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_y_axis(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_y_axis, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_y_axis_between_segments(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_y_axis_between, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_y_axis_within_segment(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_y_axis_within, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_distance_from_river(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_distance_from_river, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_quadrants_dist(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_quadrants_dist, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_estimated_dist(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_estimated_dist, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_quadrants_overlay(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_quadrants_overlay, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_quadrants_overlay_flipping(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_quadrants_overlay_flipping, 'Type', 'Spearman', 'rows', 'complete');
                    
                    corr_jrd_coviewing(subj, roi) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_coviewing, 'Type', 'Spearman', 'rows', 'complete');
                    corr_jrd_resp_similarity(subj, roi) = corr(corr_ROI_jrd_nondiag(:), mat_nondiag_resp_similarity, 'Type', 'Spearman', 'rows', 'complete');
                    
                    % Partial correlation
                    partialcorr_jrd_distances(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_distance, mat_nondiag_segments, 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_distances_cov_resp(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_distance, [mat_nondiag_coviewing, mat_nondiag_resp_similarity], 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_overlay_cov_resp(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_overlay, [mat_nondiag_coviewing, mat_nondiag_resp_similarity], 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_flipping_cov_resp(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_flipping, [mat_nondiag_coviewing, mat_nondiag_resp_similarity], 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_yaxis_cov_resp(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_y_axis, [mat_nondiag_coviewing, mat_nondiag_resp_similarity], 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_segments(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_segments, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_orthsegments(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_orthsegments, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_quadrants(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_quadrants, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_distances_overlay(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_distance, mat_nondiag_overlay, 'Type', 'Spearman', 'rows', 'complete');
                    partialcorr_jrd_overlay_distances(subj, roi) = partialcorr(corr_ROI_jrd_nondiag, mat_nondiag_overlay, mat_nondiag_distance, 'Type', 'Spearman', 'rows', 'complete');
                    
                    % Calculating the segment grouping measure (higher neural similarity within segment compared to between)
                    segment_grouping_jrd(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
                    %                     segment_grouping_newmeasure_jrd(subj, roi) = corr(corr_ROI_jrd_nondiag, mat_nondiag_segment_grouping_newmeasure, 'Type', 'Spearman', 'rows', 'complete');
                    grouping_newmeasure_factor = nanmean(mat_nondiag_distance(mat_nondiag_segments==1))/nanmean(mat_nondiag_distance(mat_nondiag_segments==0));
                    segment_grouping_newmeasure_jrd(subj, roi) = grouping_newmeasure_factor * nanmean(corr_ROI_jrd_nondiag(mat_nondiag_segments==1)) - nanmean(corr_ROI_jrd_nondiag(mat_nondiag_segments==0));
                    % Calculating the segment remapping measure (higher correlation to distance matrix within segment vs. between segments)
                    segment_remapping_jrd_full(subj, roi) = corr_jrd_dist_within_segment(subj, roi) - corr_jrd_dist_between_segments(subj, roi);
                    segment_remapping_jrd_adjquads(subj, roi) = corr_jrd_dist_within_segment_adj_quad(subj, roi) - corr_jrd_dist_between_segments_adj_quad(subj, roi);
                    
                    % Saving the actual neural correlation
                    all_neural_corr_ROI_jrd_nondiag(:,roi,subj) = corr_ROI_jrd_nondiag;
                end
            end
        end
    end
end


%% Summarizing results

n = length(all_ROI_names);

[~, pvals_means_diag_nondiag_objview_halfruns] = ttest(means_diag_minus_nondiag_objview_halfruns,[],'tail','right');
mm1_means_pv_diag_nondiag = [nanmean(means_diag_minus_nondiag_objview_halfruns); pvals_means_diag_nondiag_objview_halfruns];

[~, pvals_segment_grouping_jrd] = ttest(segment_grouping_jrd,[],'tail','right');
[~, pvals_segment_grouping_objview] = ttest(segment_grouping_objview_d2minusd1,[],'tail','right');
fdr_pvals_segment_grouping_jrd = nan(size(pvals_segment_grouping_jrd)); fdr_pvals_segment_grouping_jrd([1:4,8]) = mafdr(pvals_segment_grouping_jrd([1:4,8]),'BHFDR','true');
fdr_pvals_segment_grouping_objview = nan(size(pvals_segment_grouping_objview)); fdr_pvals_segment_grouping_objview([4:8]) = mafdr(pvals_segment_grouping_objview([4:8]),'BHFDR','true');
mm2_means_pv_segment_grouping_jrd = [nanmean(segment_grouping_jrd); fdr_pvals_segment_grouping_jrd; pvals_segment_grouping_jrd];
mm3_means_pv_segment_grouping_objview = [nanmean(segment_grouping_objview_d2minusd1); fdr_pvals_segment_grouping_objview; pvals_segment_grouping_objview];

[~, pvals_segment_remapping_jrd] = ttest(segment_remapping_jrd_full,[],'tail','right');
[~, pvals_segment_remapping_objview] = ttest(segment_remapping_objview_d2minusd1_full,[],'tail','right');
fdr_pvals_segment_remapping_jrd = nan(size(pvals_segment_remapping_jrd)); fdr_pvals_segment_remapping_jrd([1:4,8]) = mafdr(pvals_segment_remapping_jrd([1:4,8]),'BHFDR','true');
fdr_pvals_segment_remapping_objview = nan(size(pvals_segment_remapping_objview)); fdr_pvals_segment_remapping_objview([4:8]) = mafdr(pvals_segment_remapping_objview([4:8]),'BHFDR','true');
mm4_means_pv_segment_remapping_jrd = [nanmean(segment_remapping_jrd_full); fdr_pvals_segment_remapping_jrd; pvals_segment_remapping_jrd];
mm5_means_pv_segment_remapping_objview = [nanmean(segment_remapping_objview_d2minusd1_full); fdr_pvals_segment_remapping_objview; pvals_segment_remapping_objview];



[~, pvals_partialcorr_jrd_distances] = ttest(partialcorr_jrd_distances,[]);
[~, pvals_partialcorr_jrd_segments] = ttest(partialcorr_jrd_segments,[]);
[~, pvals_partialcorr_jrd_quadrants] = ttest(partialcorr_jrd_quadrants,[]);
[~, pvals_partialcorr_jrd_orthsegments] = ttest(partialcorr_jrd_orthsegments,[]);
[~, pvals_partialcorr_objview_d2minusd1_distances] = ttest(partialcorr_objview_d2minusd1_distances,[]);
[~, pvals_partialcorr_objview_d2minusd1_segments] = ttest(partialcorr_objview_d2minusd1_segments,[]);
[~, pvals_partialcorr_objview_d2minusd1_quadrants] = ttest(partialcorr_objview_d2minusd1_quadrants,[]);
[~, pvals_partialcorr_objview_d2minusd1_orthsegments] = ttest(partialcorr_objview_d2minusd1_orthsegments,[]);
pv_partialcorr1_jrd = [pvals_partialcorr_jrd_distances; pvals_partialcorr_jrd_segments; pvals_partialcorr_jrd_quadrants; pvals_partialcorr_jrd_orthsegments];
pv_partialcorr2_objview = [pvals_partialcorr_objview_d2minusd1_distances; pvals_partialcorr_objview_d2minusd1_segments; pvals_partialcorr_objview_d2minusd1_quadrants; pvals_partialcorr_objview_d2minusd1_orthsegments];
mm6_means_partialcorr_jrd = [nanmean(partialcorr_jrd_distances); nanmean(partialcorr_jrd_segments); nanmean(partialcorr_jrd_quadrants); nanmean(partialcorr_jrd_orthsegments)];
mm7_means_partialcorr_objview = [nanmean(partialcorr_objview_d2minusd1_distances); nanmean(partialcorr_objview_d2minusd1_segments); nanmean(partialcorr_objview_d2minusd1_quadrants); nanmean(partialcorr_objview_d2minusd1_orthsegments)];

means1_jrd = [nanmean(corr_jrd_distances);...
    nanmean(corr_jrd_dist_within_quadrant);nanmean(corr_jrd_dist_within_segment_adj_quad);nanmean(corr_jrd_dist_between_segments_adj_quad);nanmean(corr_jrd_dist_diagonal_quadrants);...
    nan(1,n);nanmean(corr_jrd_quadrants_dist);nanmean(corr_jrd_estimated_dist);...
    nan(1,n);nan(1,n);nanmean(corr_jrd_schema);nanmean(corr_jrd_schema_btwn);nanmean(corr_jrd_schema_y_axis_only);...
    nanmean(corr_jrd_flipping);nanmean(corr_jrd_flipping_btwn);nanmean(corr_jrd_flipping_x_axis_only);...
    nanmean(corr_jrd_overlay);nanmean(corr_jrd_overlay_btwn);nanmean(corr_jrd_overlay_x_axis_only);...
    nan(1,n);nan(1,n);nanmean(corr_jrd_segments);nanmean(corr_jrd_quadrants);nanmean(corr_jrd_orth_segments);...
    nan(1,n);nanmean(corr_jrd_y_axis);nanmean(corr_jrd_x_axis);
    nanmean(corr_jrd_y_axis_within_segment);nanmean(corr_jrd_x_axis_within_segment);nanmean(corr_jrd_y_axis_between_segments);nanmean(corr_jrd_x_axis_between_segments);...
    nan(1,n);nanmean(corr_jrd_dist_within_segment);nanmean(corr_jrd_dist_between_segments);
    nan(1,n);nanmean(corr_jrd_schema_btwn_seg_adj_quad);nanmean(corr_jrd_flipping_btwn_seg_adj_quad);nanmean(corr_jrd_overlay_btwn_seg_adj_quad)];
[~,p1] = ttest(corr_jrd_distances,[],'tail','right');
[~,p2] = ttest(corr_jrd_dist_within_quadrant,[],'tail','right');
[~,p3] = ttest(corr_jrd_dist_within_segment_adj_quad,[],'tail','right');
[~,p4] = ttest(corr_jrd_dist_between_segments_adj_quad,[],'tail','right');
[~,p5] = ttest(corr_jrd_dist_diagonal_quadrants,[],'tail','right');
[~,p6] = ttest(corr_jrd_quadrants_dist,[],'tail','right');
[~,p7] = ttest(corr_jrd_estimated_dist,[],'tail','right');
[~,p8] = ttest(corr_jrd_schema,[],'tail','right');
[~,p9] = ttest(corr_jrd_schema_btwn,[],'tail','right');
[~,p9b] = ttest(corr_jrd_schema_y_axis_only,[],'tail','right');
[~,p10] = ttest(corr_jrd_flipping,[],'tail','right');
[~,p11] = ttest(corr_jrd_flipping_btwn,[],'tail','right');
[~,p12] = ttest(corr_jrd_flipping_x_axis_only,[],'tail','right');
[~,p13] = ttest(corr_jrd_overlay,[],'tail','right');
[~,p14] = ttest(corr_jrd_overlay_btwn,[],'tail','right');
[~,p15] = ttest(corr_jrd_overlay_x_axis_only,[],'tail','right');
[~,p16] = ttest(corr_jrd_segments,[],'tail','right');
[~,p17] = ttest(corr_jrd_quadrants,[],'tail','right');
[~,p18] = ttest(corr_jrd_orth_segments,[],'tail','right');
[~,p19] = ttest(corr_jrd_y_axis,[],'tail','right');
[~,p20] = ttest(corr_jrd_x_axis,[],'tail','right');
[~,p21] = ttest(corr_jrd_y_axis_within_segment,[],'tail','right');
[~,p22] = ttest(corr_jrd_x_axis_within_segment,[],'tail','right');
[~,p23] = ttest(corr_jrd_y_axis_between_segments,[],'tail','right');
[~,p24] = ttest(corr_jrd_x_axis_between_segments,[],'tail','right');
[~,p25] = ttest(corr_jrd_dist_within_segment,[],'tail','right');
[~,p26] = ttest(corr_jrd_dist_between_segments,[],'tail','right');
[~,p27] = ttest(corr_jrd_schema_btwn_seg_adj_quad,[],'tail','right');
[~,p28] = ttest(corr_jrd_flipping_btwn_seg_adj_quad,[],'tail','right');
[~,p29] = ttest(corr_jrd_overlay_btwn_seg_adj_quad,[],'tail','right');
pp1_jrd=[p1;p2;p3;p4;p5;nan(1,n);p6;p7;nan(1,n);nan(1,n);p8;p9;p9b;p10;p11;p12;p13;p14;p15;nan(1,n);nan(1,n);p16;p17;p18;nan(1,n);p19;p20;p21;p22;p23;p24;nan(1,n);p25;p26;nan(1,n);p27;p28;p29];
fdr_pp1_jrd = nan(size(pp1_jrd));       % Calculating FDR-corrected p-values
for i=1:size(fdr_pp1_jrd,1)
    fdr_pp1_jrd(i,[1:4,8]) = mafdr(pp1_jrd(i,[1:4,8]),'BHFDR','true');
end

means2_objview_d2minusd1 = [nanmean(corr_objview_d2minusd1_distances);...
    nanmean(corr_objview_d2minusd1_dist_within_quadrant);nanmean(corr_objview_d2minusd1_dist_within_segment_adj_quad);nanmean(corr_objview_d2minusd1_dist_between_segments_adj_quad);nanmean(corr_objview_d2minusd1_dist_diagonal_quadrants);...
    nan(1,n);nanmean(corr_objview_d2minusd1_quadrants_dist);nanmean(corr_objview_d2minusd1_estimated_dist);...
    nan(1,n);nan(1,n);nanmean(corr_objview_d2minusd1_schema);nanmean(corr_objview_d2minusd1_schema_btwn);nanmean(corr_objview_d2minusd1_schema_y_axis_only);...
    nanmean(corr_objview_d2minusd1_flipping);nanmean(corr_objview_d2minusd1_flipping_btwn);nanmean(corr_objview_d2minusd1_flipping_x_axis_only);...
    nanmean(corr_objview_d2minusd1_overlay);nanmean(corr_objview_d2minusd1_overlay_btwn);nanmean(corr_objview_d2minusd1_overlay_x_axis_only);...
    nan(1,n);nan(1,n);nanmean(corr_objview_d2minusd1_segments);nanmean(corr_objview_d2minusd1_quadrants);nanmean(corr_objview_d2minusd1_orth_segments);...
    nan(1,n);nanmean(corr_objview_d2minusd1_y_axis);nanmean(corr_objview_d2minusd1_x_axis);
    nanmean(corr_objview_d2minusd1_y_axis_within_segment);nanmean(corr_objview_d2minusd1_x_axis_within_segment);nanmean(corr_objview_d2minusd1_y_axis_between_segments);nanmean(corr_objview_d2minusd1_x_axis_between_segments);...
    nan(1,n);nanmean(corr_objview_d2minusd1_dist_within_segment);nanmean(corr_objview_d2minusd1_dist_between_segments);
    nan(1,n);nanmean(corr_objview_d2minusd1_schema_btwn_seg_adj_quad);nanmean(corr_objview_d2minusd1_flipping_btwn_seg_adj_quad);nanmean(corr_objview_d2minusd1_overlay_btwn_seg_adj_quad)];
[~,p1] = ttest(corr_objview_d2minusd1_distances,[],'tail','right');
[~,p2] = ttest(corr_objview_d2minusd1_dist_within_quadrant,[],'tail','right');
[~,p3] = ttest(corr_objview_d2minusd1_dist_within_segment_adj_quad,[],'tail','right');
[~,p4] = ttest(corr_objview_d2minusd1_dist_between_segments_adj_quad,[],'tail','right');
[~,p5] = ttest(corr_objview_d2minusd1_dist_diagonal_quadrants,[],'tail','right');
[~,p6] = ttest(corr_objview_d2minusd1_quadrants_dist,[],'tail','right');
[~,p7] = ttest(corr_objview_d2minusd1_estimated_dist,[],'tail','right');
[~,p8] = ttest(corr_objview_d2minusd1_schema,[],'tail','right');
[~,p9] = ttest(corr_objview_d2minusd1_schema_btwn,[],'tail','right');
[~,p9b] = ttest(corr_objview_d2minusd1_schema_y_axis_only,[],'tail','right');
[~,p10] = ttest(corr_objview_d2minusd1_flipping,[],'tail','right');
[~,p11] = ttest(corr_objview_d2minusd1_flipping_btwn,[],'tail','right');
[~,p12] = ttest(corr_objview_d2minusd1_flipping_x_axis_only,[],'tail','right');
[~,p13] = ttest(corr_objview_d2minusd1_overlay,[],'tail','right');
[~,p14] = ttest(corr_objview_d2minusd1_overlay_btwn,[],'tail','right');
[~,p15] = ttest(corr_objview_d2minusd1_overlay_x_axis_only,[],'tail','right');
[~,p16] = ttest(corr_objview_d2minusd1_segments,[],'tail','right');
[~,p17] = ttest(corr_objview_d2minusd1_quadrants,[],'tail','right');
[~,p18] = ttest(corr_objview_d2minusd1_orth_segments,[],'tail','right');
[~,p19] = ttest(corr_objview_d2minusd1_y_axis,[],'tail','right');
[~,p20] = ttest(corr_objview_d2minusd1_x_axis,[],'tail','right');
[~,p21] = ttest(corr_objview_d2minusd1_y_axis_within_segment,[],'tail','right');
[~,p22] = ttest(corr_objview_d2minusd1_x_axis_within_segment,[],'tail','right');
[~,p23] = ttest(corr_objview_d2minusd1_y_axis_between_segments,[],'tail','right');
[~,p24] = ttest(corr_objview_d2minusd1_x_axis_between_segments,[],'tail','right');
[~,p25] = ttest(corr_objview_d2minusd1_dist_within_segment,[],'tail','right');
[~,p26] = ttest(corr_objview_d2minusd1_dist_between_segments,[],'tail','right');
[~,p27] = ttest(corr_objview_d2minusd1_schema_btwn_seg_adj_quad,[],'tail','right');
[~,p28] = ttest(corr_objview_d2minusd1_flipping_btwn_seg_adj_quad,[],'tail','right');
[~,p29] = ttest(corr_objview_d2minusd1_overlay_btwn_seg_adj_quad,[],'tail','right');
pp1_objview_d2minusd1=[p1;p2;p3;p4;p5;nan(1,n);p6;p7;nan(1,n);nan(1,n);p8;p9;p9b;p10;p11;p12;p13;p14;p15;nan(1,n);nan(1,n);p16;p17;p18;nan(1,n);p19;p20;p21;p22;p23;p24;nan(1,n);p25;p26;nan(1,n);p27;p28;p29];
fdr_pp1_objview_d2minusd1 = nan(size(pp1_objview_d2minusd1));       % Calculating FDR-corrected p-values
for i=1:size(fdr_pp1_objview_d2minusd1,1)
    fdr_pp1_objview_d2minusd1(i,[4:8]) = mafdr(pp1_objview_d2minusd1(i,[4:8]),'BHFDR','true');
end

means3_objview_day2 = [nanmean(corr_objview_day2_distances);...
    nanmean(corr_objview_day2_dist_within_quadrant);nanmean(corr_objview_day2_dist_within_segment_adj_quad);nanmean(corr_objview_day2_dist_between_segments_adj_quad);nanmean(corr_objview_day2_dist_diagonal_quadrants);...
    nan(1,n);nanmean(corr_objview_day2_quadrants_dist);nanmean(corr_objview_day2_estimated_dist);...
    nan(1,n);nan(1,n);nanmean(corr_objview_day2_schema);nanmean(corr_objview_day2_schema_btwn);nanmean(corr_objview_day2_schema_y_axis_only);...
    nanmean(corr_objview_day2_flipping);nanmean(corr_objview_day2_flipping_btwn);nanmean(corr_objview_day2_flipping_x_axis_only);...
    nanmean(corr_objview_day2_overlay);nanmean(corr_objview_day2_overlay_btwn);nanmean(corr_objview_day2_overlay_x_axis_only);...
    nan(1,n);nan(1,n);nanmean(corr_objview_day2_segments);nanmean(corr_objview_day2_quadrants);nanmean(corr_objview_day2_orth_segments);...
    nan(1,n);nanmean(corr_objview_day2_y_axis);nanmean(corr_objview_day2_x_axis);
    nanmean(corr_objview_day2_y_axis_within_segment);nanmean(corr_objview_day2_x_axis_within_segment);nanmean(corr_objview_day2_y_axis_between_segments);nanmean(corr_objview_day2_x_axis_between_segments);...
    nan(1,n);nanmean(corr_objview_day2_dist_within_segment);nanmean(corr_objview_day2_dist_between_segments);
    nan(1,n);nanmean(corr_objview_day2_schema_btwn_seg_adj_quad);nanmean(corr_objview_day2_flipping_btwn_seg_adj_quad);nanmean(corr_objview_day2_overlay_btwn_seg_adj_quad)];
[~,p1] = ttest(corr_objview_day2_distances,[],'tail','right');
[~,p2] = ttest(corr_objview_day2_dist_within_quadrant,[],'tail','right');
[~,p3] = ttest(corr_objview_day2_dist_within_segment_adj_quad,[],'tail','right');
[~,p4] = ttest(corr_objview_day2_dist_between_segments_adj_quad,[],'tail','right');
[~,p5] = ttest(corr_objview_day2_dist_diagonal_quadrants,[],'tail','right');
[~,p6] = ttest(corr_objview_day2_quadrants_dist,[],'tail','right');
[~,p7] = ttest(corr_objview_day2_estimated_dist,[],'tail','right');
[~,p8] = ttest(corr_objview_day2_schema,[],'tail','right');
[~,p9] = ttest(corr_objview_day2_schema_btwn,[],'tail','right');
[~,p9b] = ttest(corr_objview_day2_schema_y_axis_only,[],'tail','right');
[~,p10] = ttest(corr_objview_day2_flipping,[],'tail','right');
[~,p11] = ttest(corr_objview_day2_flipping_btwn,[],'tail','right');
[~,p12] = ttest(corr_objview_day2_flipping_x_axis_only,[],'tail','right');
[~,p13] = ttest(corr_objview_day2_overlay,[],'tail','right');
[~,p14] = ttest(corr_objview_day2_overlay_btwn,[],'tail','right');
[~,p15] = ttest(corr_objview_day2_overlay_x_axis_only,[],'tail','right');
[~,p16] = ttest(corr_objview_day2_segments,[],'tail','right');
[~,p17] = ttest(corr_objview_day2_quadrants,[],'tail','right');
[~,p18] = ttest(corr_objview_day2_orth_segments,[],'tail','right');
[~,p19] = ttest(corr_objview_day2_y_axis,[],'tail','right');
[~,p20] = ttest(corr_objview_day2_x_axis,[],'tail','right');
[~,p21] = ttest(corr_objview_day2_y_axis_within_segment,[],'tail','right');
[~,p22] = ttest(corr_objview_day2_x_axis_within_segment,[],'tail','right');
[~,p23] = ttest(corr_objview_day2_y_axis_between_segments,[],'tail','right');
[~,p24] = ttest(corr_objview_day2_x_axis_between_segments,[],'tail','right');
[~,p25] = ttest(corr_objview_day2_dist_within_segment,[],'tail','right');
[~,p26] = ttest(corr_objview_day2_dist_between_segments,[],'tail','right');
[~,p27] = ttest(corr_objview_day2_schema_btwn_seg_adj_quad,[],'tail','right');
[~,p28] = ttest(corr_objview_day2_flipping_btwn_seg_adj_quad,[],'tail','right');
[~,p29] = ttest(corr_objview_day2_overlay_btwn_seg_adj_quad,[],'tail','right');
pp1_objview_day2=[p1;p2;p3;p4;p5;nan(1,n);p6;p7;nan(1,n);nan(1,n);p8;p9;p9b;p10;p11;p12;p13;p14;p15;nan(1,n);nan(1,n);p16;p17;p18;nan(1,n);p19;p20;p21;p22;p23;p24;nan(1,n);p25;p26;nan(1,n);p27;p28;p29];
fdr_pp1_objview_day2 = nan(size(pp1_objview_day2));       % Calculating FDR-corrected p-values
for i=1:size(fdr_pp1_objview_day2,1)
    fdr_pp1_objview_day2(i,[4:8]) = mafdr(pp1_objview_day2(i,[4:8]),'BHFDR','true');
end


% All JRD distance data, in one column
all_dists_jrd_column = []; all_s = {}; counter = 1;
for subj=1:length(subj_data_dirs)
    for roi=1:num_ROIs/2
        all_dists_jrd_column = [all_dists_jrd_column; 
            corr_jrd_distances(subj,roi);
            corr_jrd_y_axis(subj,roi); corr_jrd_x_axis(subj,roi);
            corr_jrd_y_axis_within_segment(subj,roi); corr_jrd_x_axis_within_segment(subj,roi);
            corr_jrd_y_axis_between_segments(subj,roi); corr_jrd_x_axis_between_segments(subj,roi);
            corr_jrd_schema_y_axis_only(subj,roi);
            corr_jrd_flipping_x_axis_only(subj,roi); corr_jrd_overlay_x_axis_only(subj,roi)];
        for i=1:10, all_s{counter} = [subj_data_dirs(subj).name]; counter=counter+1;end
    end
end
all_s=all_s';
all_dists_group = [nanmean(corr_jrd_distances); 
    nanmean(corr_jrd_y_axis); nanmean(corr_jrd_x_axis);
    nanmean(corr_jrd_y_axis_within_segment); nanmean(corr_jrd_x_axis_within_segment);
    nanmean(corr_jrd_y_axis_between_segments); nanmean(corr_jrd_x_axis_between_segments);
    nanmean(corr_jrd_schema_y_axis_only);
    nanmean(corr_jrd_flipping_x_axis_only); nanmean(corr_jrd_overlay_x_axis_only)];
all_dists_group = all_dists_group(:,1:num_ROIs/2);
all_dists_group = all_dists_group(:);


%% Calculate noise ceiling
% The noise ceiling is the average correlation between subject-level neural
% pattern similarity matrices and the neural similariry matrix averaged
% across subjects; the assumption is that the group-average matrix is the
% best possible model that could be achieved, so the subject-level
% matricees correlation to the group matrix reflect the maximum amount of 
% precision possible, taking into account noise in the data


% Rank-transforming the original dissimilarity matrices (Nili & Kriegeskorte 2014)
all_transformed_neural_corr_ROI_jrd_nondiag = tiedrank(all_neural_corr_ROI_jrd_nondiag);
all_transformed_neural_corr_objview_day2_nondiag = tiedrank(all_neural_corr_objview_day2_nondiag);
all_transformed_neural_corr_objview_d2minusd1_nondiag = tiedrank(all_neural_corr_objview_d2minusd1_nondiag);


noise_ceiling_JRD = nan(2, num_ROIs);
noise_ceiling_objview_day2 = nan(2, num_ROIs);
noise_ceiling_objview_d2minusd1 = nan(2, num_ROIs);
for roi = 1:num_ROIs
    current_ROI_corrs_JRD = nan(length(subj_data_dirs), 2);
    current_ROI_corrs_objview_day2 = nan(length(subj_data_dirs), 2);
    current_ROI_corrs_objview_d2minusd1 = nan(length(subj_data_dirs), 2);
    
    % Calculating group mean 
    for subj = 1:length(subj_data_dirs)
        % Calculating all subjects' mean correlation matrix
        allsubjs_mean_corr_JRD = nanmean(all_transformed_neural_corr_ROI_jrd_nondiag(:,roi,:), 3);
        allsubjs_mean_corr_objview_day2 = nanmean(all_transformed_neural_corr_objview_day2_nondiag(:,roi,:), 3);
        allsubjs_mean_corr_objview_d2minusd1 = nanmean(all_transformed_neural_corr_objview_d2minusd1_nondiag(:,roi,:), 3);
        % Calculating mean correlation matrix of all subjects except for the current subject
        allsubjs_leaveoneout_mean_corr_JRD = nanmean(all_transformed_neural_corr_ROI_jrd_nondiag(:,roi,[1:subj-1, subj+1:end]), 3);
        allsubjs_leaveoneout_mean_corr_objview_day2 = nanmean(all_transformed_neural_corr_objview_day2_nondiag(:,roi,[1:subj-1, subj+1:end]), 3);
        allsubjs_leaveoneout_mean_corr_objview_d2minusd1 = nanmean(all_transformed_neural_corr_objview_d2minusd1_nondiag(:,roi,[1:subj-1, subj+1:end]), 3);
        
        % Upper bound - correlation of current matrix to mean of all subjects' matrices
        current_ROI_corrs_JRD(subj, 1) = corr(all_transformed_neural_corr_ROI_jrd_nondiag(:,roi,subj), allsubjs_mean_corr_JRD, 'Type', 'Spearman', 'rows', 'complete');
        current_ROI_corrs_objview_day2(subj, 1) = corr(all_transformed_neural_corr_objview_day2_nondiag(:,roi,subj), allsubjs_mean_corr_objview_day2, 'Type', 'Spearman', 'rows', 'complete');
        current_ROI_corrs_objview_d2minusd1(subj, 1) = corr(all_transformed_neural_corr_objview_d2minusd1_nondiag(:,roi,subj), allsubjs_mean_corr_objview_d2minusd1, 'Type', 'Spearman', 'rows', 'complete');
        % Lower bound - correlation of current matrix to mean of all subjects' matrices, excluding this subject
        current_ROI_corrs_JRD(subj, 2) = corr(all_transformed_neural_corr_ROI_jrd_nondiag(:,roi,subj), allsubjs_leaveoneout_mean_corr_JRD, 'Type', 'Spearman', 'rows', 'complete');
        current_ROI_corrs_objview_day2(subj, 2) = corr(all_transformed_neural_corr_objview_day2_nondiag(:,roi,subj), allsubjs_leaveoneout_mean_corr_objview_day2, 'Type', 'Spearman', 'rows', 'complete');
        current_ROI_corrs_objview_d2minusd1(subj, 2) = corr(all_transformed_neural_corr_objview_d2minusd1_nondiag(:,roi,subj), allsubjs_leaveoneout_mean_corr_objview_d2minusd1, 'Type', 'Spearman', 'rows', 'complete');
    end
    noise_ceiling_JRD(:,roi) = nanmean(current_ROI_corrs_JRD)';
    noise_ceiling_objview_day2(:,roi) = nanmean(current_ROI_corrs_objview_day2)';
    noise_ceiling_objview_d2minusd1(:,roi) = nanmean(current_ROI_corrs_objview_d2minusd1)';
end


