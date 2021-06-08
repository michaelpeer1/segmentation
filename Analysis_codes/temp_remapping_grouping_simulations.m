% Temporary simulations

num_conditions = 16;

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

% Schema preserving - with distance
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



% Flipping - with distance
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
mat_nondiag_distance_from_river = mat_distance_from_river(find(tril(ones(num_conditions),-1)));
mat_nondiag_quadrants_dist = mat_quadrants_dist(find(tril(ones(num_conditions),-1)));
mat_nondiag_within_seg_adj_quad = mat_within_seg_adj_quad(find(tril(ones(num_conditions),-1)));
mat_nondiag_between_seg_adj_quad = mat_between_seg_adj_quad(find(tril(ones(num_conditions),-1)));
mat_nondiag_segment_grouping = mat_segment_grouping(find(tril(ones(num_conditions),-1)));

%% Simulations

num_perms = 100;
grouping_factor = 50;
noise_factor = 0.5;

grouping_index = nan(num_perms, 3);
new_grouping_index = nan(num_perms, 3);
remapping_index = nan(num_perms, 3);

% new_grouping_mat = ones(num_conditions);
% new_grouping_mat(mat_segments == 1) = 1.5052;
% new_grouping_mat(mat_segments == 0) = 1;
% new_grouping_mat = new_grouping_mat(find(tril(ones(num_conditions),-1)));

for i = 1:num_perms
    remapping_sim = mat_nondiag_distance;
    remapping_sim(mat_nondiag_segments==0) = rand(1,sum(mat_nondiag_segments(:)==0));
    
    locs_new = locations_all_objects;
    locs_new([1:4,13:16],2) = locs_new([1:4,13:16],2) + grouping_factor;    % Making segments further away from each other
    mat_new = pdist2(locs_new, locs_new);         % Distances between all objects
    mat_new = mat_new / max(mat_new(:));    % Normalizing to range 0-1
    mat_new = 1 - mat_new; 
    grouping_sim = mat_new(find(tril(ones(num_conditions),-1)));
%     grouping_sim(mat_nondiag_segments==1) = grouping_sim(mat_nondiag_segments==1) + grouping_factor;
    
    noise_sim = mat_nondiag_distance;
    noise_sim = noise_sim + (rand(length(noise_sim),1)-0.5)*noise_factor;
    
    grouping_index(i, 1) = corr(remapping_sim, mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
    grouping_index(i, 2) = corr(grouping_sim, mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
    grouping_index(i, 3) = corr(noise_sim, mat_nondiag_segment_grouping, 'Type', 'Spearman', 'rows', 'complete');
        
%     new_grouping_index(i, 1) = corr(remapping_sim, new_grouping_mat, 'Type', 'Spearman', 'rows', 'complete');
%     new_grouping_index(i, 2) = corr(grouping_sim, new_grouping_mat, 'Type', 'Spearman', 'rows', 'complete');
%     new_grouping_index(i, 3) = corr(noise_sim, new_grouping_mat, 'Type', 'Spearman', 'rows', 'complete');         
    new_grouping_index(i, 1) = nanmean(remapping_sim(mat_nondiag_segments==1)) - 1.505*nanmean(remapping_sim(mat_nondiag_segments==0));
    new_grouping_index(i, 2) = nanmean(grouping_sim(mat_nondiag_segments==1)) - 1.505*nanmean(grouping_sim(mat_nondiag_segments==0));
    new_grouping_index(i, 3) = nanmean(noise_sim(mat_nondiag_segments==1)) - 1.505*nanmean(noise_sim(mat_nondiag_segments==0));
    
    remapping_index(i,1) = corr(remapping_sim, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete') - corr(remapping_sim, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
    remapping_index(i,2) = corr(grouping_sim, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete') - corr(grouping_sim, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
    remapping_index(i,3) = corr(noise_sim, mat_nondiag_dist_within_segment, 'Type', 'Spearman', 'rows', 'complete') - corr(noise_sim, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
    
    % Permutation test for a new remapping measure
    n = 100; % Number of permutations
    real_between_corr1 = corr(remapping_sim, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
    real_between_corr2 = corr(grouping_sim, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
    real_between_corr3 = corr(noise_sim, mat_nondiag_dist_between_segments, 'Type', 'Spearman', 'rows', 'complete');
    perm_corrs1 = [];    perm_corrs2 = [];    perm_corrs3 = [];
    for p = 1:n
        m = mat_nondiag_segments(randperm(length(mat_nondiag_segments)));   % Randomizing the within-between segment assignment
        mm = mat_nondiag_distance; mm(m==0) = nan;      % Removing the randomly assigned between-segment locations
        perm_corrs1(p) = corr(remapping_sim, mm, 'Type', 'Spearman', 'rows', 'complete');
        perm_corrs2(p) = corr(grouping_sim, mm, 'Type', 'Spearman', 'rows', 'complete');
        perm_corrs3(p) = corr(noise_sim, mm, 'Type', 'Spearman', 'rows', 'complete');
    end
    new_remapping_index(i,1) = sum(real_between_corr1 <= perm_corrs1) / n;
    new_remapping_index(i,2) = sum(real_between_corr2 <= perm_corrs2) / n;
    new_remapping_index(i,3) = sum(real_between_corr3 <= perm_corrs3) / n;
end


% Measuring how noise affects each model, and which models are easier to detect
p1=[];for i=1:10000,p1(i)=corr(mat_nondiag_distance, mat_nondiag_distance+rand(size(mat_nondiag_distance))*5,'Type', 'Kendall', 'rows', 'complete');end
p2=[];for i=1:10000,p2(i)=corr(mat_nondiag_overlay, mat_nondiag_overlay+rand(size(mat_nondiag_distance))*5,'Type', 'Kendall', 'rows', 'complete');end
p3=[];for i=1:10000,p3(i)=corr(mat_nondiag_flipping, mat_nondiag_flipping+rand(size(mat_nondiag_distance))*5,'Type', 'Kendall', 'rows', 'complete');end
p4=[];for i=1:10000,p4(i)=corr(mat_nondiag_schema, mat_nondiag_schema+rand(size(mat_nondiag_distance))*5,'Type', 'Kendall', 'rows', 'complete');end
p5=[];for i=1:10000,p5(i)=corr(mat_nondiag_segment_grouping(~isnan(mat_nondiag_segment_grouping)), mat_nondiag_segment_grouping(~isnan(mat_nondiag_segment_grouping))+rand(size(mat_nondiag_distance(~isnan(mat_nondiag_segment_grouping))))*5,'Type', 'Kendall', 'rows', 'complete');end
mat_combined_distance_overlay = mean([mat_nondiag_distance,mat_nondiag_overlay],2);
pc1=corr(mat_combined_distance_overlay, mat_nondiag_distance,'Type', 'Spearman', 'rows', 'complete');
pc2=corr(mat_combined_distance_overlay, mat_nondiag_overlay,'Type', 'Spearman', 'rows', 'complete');
% std(mat_nondiag_distance)
% std(mat_nondiag_overlay)
% nanstd(mat_nondiag_segment_grouping)
