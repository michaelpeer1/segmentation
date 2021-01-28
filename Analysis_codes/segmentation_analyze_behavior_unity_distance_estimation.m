% Distance estimation analysis

% Pilot data
% filenames = {};
% filenames{1} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1329.csv';  % MP
% filenames{2} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1629.csv';  % RM
% filenames{3} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1754.csv';  % ZL
% filenames{4} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1304_corrected.csv';  % AT - partial

% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v1_river\Subjects_file.xlsx');
% num_subjects = max(find(~cellfun(@isnan, subjects_file(4:end,1))));
[~, ~, subjects_file] = xlsread('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('/Users/mpeer/Desktop/Segmentation_data/fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v2_wall\Subjects_file_v2_wall.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_behavioral_pilot\Subjects_file_fMRI_behav_pilot.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_experiment\Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\5_behavioral_wall_new\Subjects_file_v5_wall_new.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\6_separate_learning\Subjects_file_v6_separate_learning.xlsx');
num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
filenames = {};
for i=1:num_subjects
    %     if ~isnan(subjects_file{3+i, 13})       % Behavioral version of the task
    %         filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 13}) '.csv'];
    if ~isnan(subjects_file{3+i, 13})       % fMRI version of the task
        %         filenames{end+1} = [fullfile(subjects_file{3+i, 8}, subjects_file{3+i, 13}) '.csv'];
        filenames{end+1} = [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 13}) '.csv'];
    end
end

answers_all_subjects = cell(1, length(filenames));

all_RTs_within_segment = zeros(length(filenames), 1);
all_RTs_between_segments = zeros(length(filenames), 1);
all_scaled_est_dists_within_segment = zeros(length(filenames), 1);
all_scaled_est_dists_between_segments = zeros(length(filenames), 1);
all_scaled_real_dists_within_segment = zeros(length(filenames), 1);
all_scaled_real_dists_between_segments = zeros(length(filenames), 1);
all_abs_error_within_segment = zeros(length(filenames), 1);
all_abs_error_between_segments = zeros(length(filenames), 1);
all_abs_error_overall = zeros(length(filenames), 1);
all_error_within_segment = zeros(length(filenames), 1);
all_error_between_segments = zeros(length(filenames), 1);
all_corr_real_est_overall = zeros(length(filenames), 1);
all_corr_real_est_within = zeros(length(filenames), 1);
all_corr_real_est_between = zeros(length(filenames), 1);
all_est_dists_within_segment = zeros(length(filenames), 1);
all_est_dists_between_segments = zeros(length(filenames), 1);
all_p_within_between_est_dists = zeros(length(filenames), 1);
all_p_within_between_average_error = zeros(length(filenames), 1);
all_p_within_between_average_abs_error = zeros(length(filenames), 1);
all_std_est_dist_within = zeros(length(filenames), 1);
all_std_est_dist_between = zeros(length(filenames), 1);
all_p_RTs_within_between = zeros(length(filenames), 1);
all_RTs_same_within_segment_answer = zeros(length(filenames), 1);
all_RTs_different_within_segment_answer = zeros(length(filenames), 1);
all_RTs_same_within_orth_segments_answer = zeros(length(filenames), 1);
all_RTs_different_within_orth_segments_answer = zeros(length(filenames), 1);
all_task_overall_times = zeros(length(filenames), 1);
all_num_answers = zeros(length(filenames), 1);
all_abs_path_error_within_segment = zeros(length(filenames), 1);
all_abs_path_error_between_segments = zeros(length(filenames), 1);
all_abs_path_error_overall = zeros(length(filenames), 1);
all_corr_path_est_overall = zeros(length(filenames), 1);
all_corr_path_est_within = zeros(length(filenames), 1);
all_corr_path_est_between = zeros(length(filenames), 1);
all_p_within_between_average_path_abs_error = zeros(length(filenames), 1);
all_corr_overlay_est_overall = zeros(length(filenames), 1);
all_corr_overlay_est_within = zeros(length(filenames), 1);
all_corr_overlay_est_between = zeros(length(filenames), 1);

for f = 1:length(filenames)
    %% Reading the data
    disp(strcat("File ", num2str(f)));
    %     [~,~,results] = xlsread(filenames{f});
    fid = fopen(filenames{f});
    data = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:length(results(:))        % Changing strings to numbers
        if ~isnan(str2double(results{i}))
            results{i} = str2double(results{i});
        end
    end
    
    locations_trials = [];
    for i=1:size(results, 1)
        if strcmp(results{i,1}, 'Current question is')
            locations_trials(end+1) = i;
        end
    end
    num_trials = length(locations_trials);
    locations_trials(end+1) = size(results,1);
    
    
    all_answers = nan(num_trials, 18);   % 1.quadrant1, 2.quadrant2, 3.object1, 4.object2, 5.time input started, 6.time input entered, 7.input, 8.between_within, 9.real_distance, 10.scaled_estimated_distance, 11.scaled_real_distance, 12.absolute_scaled_error, 13.actual_scaled_error, 14.path_distance, 15.scaled_path_distance, 16.absolute_scaled_path_error, 17.overlay_distance,18.scaled_overlay_distance
    
    for i=1:num_trials
        % Entering question parameters into the list
        all_answers(i,1:4) = cell2mat(results(locations_trials(i), 6:9));
        
        % Getting the other parameters from the next lines, if existing
        for j = locations_trials(i):locations_trials(i+1)
            if strcmp(results{j,1}, 'New input started')
                all_answers(i, 5) = results{j, 12} - results{locations_trials(i), 12};
            elseif strcmp(results{j,1}, 'Enter key pressed; Player input is')
                all_answers(i, 6) = results{j, 12} - results{locations_trials(i), 12};
                all_answers(i, 7) = results{j, 5};
            end
        end
        
        % Figuring out if comparison was between segments or within
        % 0 1 and 1 0 and 2 3 and 3 2 are between segments, 1 2 and 2 1 and 3 0 and 0 3 are within
        q1 = all_answers(i,1); q2 = all_answers(i,2);
        if ((q1 == 1 && q2 == 2) || (q1 == 2 && q2 == 1) || (q1 == 3 && q2 == 0) || (q1 == 0 && q2 == 3))
            all_answers(i,8) = 1;   % within segment
        elseif ((q1 == 0 && q2 == 1) || (q1 == 1 && q2 == 0) || (q1 == 2 && q2 == 3) || (q1 == 3 && q2 == 2))
            all_answers(i,8) = -1;   % between segments
        else
            all_answers(i,8) = 0;   % Other (within quadrant, across diagonal)
        end
        
        % Removing trials with response 0
        if all_answers(i,7) == 0
            for q=5:11, all_answers(i,q) = nan; end
        end
    end
    
    
    locs_within = find(all_answers(:,8) == 1);
    locs_between = find(all_answers(:,8) == -1);
    locs_within = locs_within(~isnan(all_answers(locs_within,7)) & all_answers(locs_within,7)~=0);
    locs_between = locs_between(~isnan(all_answers(locs_between,7)) & all_answers(locs_between,7)~=0);
    all_relevant_locations = [locs_within; locs_between];
    
    
    %% Response times
    % Removing outliers (3 SDs)
    all_current_RTs = all_answers(all_relevant_locations, 6);
    outlier_thresh_high = mean(all_current_RTs) + 3*std(all_current_RTs); outlier_thresh_low = mean(all_current_RTs) - 3*std(all_current_RTs);
    current_RTs_within = all_answers(locs_within, 6); current_RTs_between = all_answers(locs_between, 6);
    current_RTs_within = current_RTs_within(current_RTs_within < outlier_thresh_high & current_RTs_within > outlier_thresh_low);
    current_RTs_between = current_RTs_between(current_RTs_between < outlier_thresh_high & current_RTs_between > outlier_thresh_low);
    
    % Difference between within / between segment comparisons in RT to rotation start / final answer
    disp("RT within between")
    [~, p] = ttest2(current_RTs_within, current_RTs_between);
    disp([nanmean(current_RTs_within) nanmean(current_RTs_between) p])
    
    all_RTs_within_segment(f) = nanmean(current_RTs_within);
    all_RTs_between_segments(f) = nanmean(current_RTs_between);
    all_p_RTs_within_between(f) = p;
    
    
    %% Distances
    
    % Calculating a distance matrix based on real Euclidean distances
    %     locations_all_objects = [16 35; 32 63; 65 45; 41 17; ...
    %         35 -16; 63 -32; 45 -65; 17 -41; ...
    %         -16 -35; -32 -63; -65 -45; -41 -17; ...
    %         -35 16; -63 32; -45 65; -17 41];
    locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
        42 -16; 58 -50; 27 -57; 17 -25; ...
        -16 -42; -50 -58; -57 -27; -25 -17; ...
        -42 16; -58 50; -27 57; -17 25];
    all_real_distances = pdist2(locations_all_objects, locations_all_objects);
    
    % Calculating a distance matrix based on path distances
    segment=[1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1];
    bridge1_coords = [50 10; 50 -10];
    bridge2_coords = [-50 10; -50 -10];
    length_bridges = 20;
    all_path_distances = all_real_distances;
    for i=1:size(locations_all_objects, 1)
        for j=1:size(locations_all_objects, 1)
            if segment(i) ~= segment(j)
                dist_bridge1 = sqrt(sum((locations_all_objects(i,:) - bridge1_coords(segment(i),:)).^2)) + ...
                    length_bridges + ...
                    sqrt(sum((locations_all_objects(j,:) - bridge1_coords(segment(j),:)).^2));
                dist_bridge2 = sqrt(sum((locations_all_objects(i,:) - bridge2_coords(segment(i),:)).^2)) + ...
                    length_bridges + ...
                    sqrt(sum((locations_all_objects(j,:) - bridge2_coords(segment(j),:)).^2));
                all_path_distances(i,j) = min(dist_bridge1, dist_bridge2);
            end
        end
    end
    
    % Calculating a distance matrix based on an overlay model
    locations_new_overlay = locations_all_objects;
    locations_new_overlay(5:12,2) = locations_new_overlay(5:12,2) + 75;
    all_overlay_distances = pdist2(locations_new_overlay, locations_new_overlay);         % Distances between all objects
    
    
    for i=1:num_trials
        if ~isnan(all_answers(i,7))
            all_answers(i,9) = all_real_distances(all_answers(i,3)+1, all_answers(i,4)+1);
            all_answers(i,14) = all_path_distances(all_answers(i,3)+1, all_answers(i,4)+1);
            all_answers(i,17) = all_overlay_distances(all_answers(i,3)+1, all_answers(i,4)+1);
        end
    end
    
    % scaling all estimated distances
    %     all_answers(:,10) = all_answers(:,7) - min(all_answers(all_relevant_locations, 7));
    %     all_answers(:,10) = all_answers(:,10) / max(all_answers(all_relevant_locations, 10)) * 100;     % Scaling to the range 0-100
    all_answers(all_relevant_locations, 10) = zscore(all_answers(all_relevant_locations, 7));     % Scaling using zscore
    
    % scaling all real distances
    %     all_answers(:,11) = all_answers(:,9) - min(all_answers(all_relevant_locations, 9));
    %     all_answers(:,11) = all_answers(:,11) / max(all_answers(all_relevant_locations, 11)) * 100;     % Scaling to the range 0-100
    all_answers(all_relevant_locations, 11) = zscore(all_answers(all_relevant_locations, 9));     % Scaling using zscore
    all_answers(all_relevant_locations, 15) = zscore(all_answers(all_relevant_locations, 14));     % Scaling using zscore
    all_answers(all_relevant_locations, 18) = zscore(all_answers(all_relevant_locations, 17));     % Scaling using zscore
    
    disp("Average estimated distance - within between")
    % Difference between within / between segment comparisons in estimated distance
    [~, p] = ttest2(all_answers(locs_within, 10), all_answers(locs_between, 10));
    disp([nanmean(all_answers(locs_within, 10)) nanmean(all_answers(locs_between, 10)) p])
    
    
    % Correlation between real and estimated distances
    actual_corr = corr(all_answers(all_relevant_locations,10), all_answers(all_relevant_locations,11));
    actual_path_corr = corr(all_answers(all_relevant_locations,10), all_answers(all_relevant_locations,15));
    actual_overlay_corr = corr(all_answers(all_relevant_locations,10), all_answers(all_relevant_locations,18));
    % Permutation test to check significance
    num_iters = 1000; rand_c = zeros(1, num_iters); rand_path_c = zeros(1, num_iters);
    a = all_answers(all_relevant_locations,10); b = all_answers(all_relevant_locations, 11); c = all_answers(all_relevant_locations, 15);
    for i=1:num_iters
        rand_c(i) = corr(a(randperm(length(a))), b);
        rand_path_c(i) = corr(a(randperm(length(a))), c);
    end
    p_rand = sum(actual_corr < rand_c) / num_iters;
    p_rand_path = sum(actual_path_corr < rand_path_c) / num_iters;
    
    disp("Corr real estimated distances - all relevant locations")
    %     figure; scatter(all_answers(all_relevant_locations,10), all_answers(all_relevant_locations,11))
    disp([actual_corr p_rand]);
    disp("Corr real estimated distances - within between")
    disp([corr(all_answers(locs_within,10),all_answers(locs_within,11)) corr(all_answers(locs_between,10),all_answers(locs_between,11))]);
    
    
    % Errors
    all_answers(:,12) = abs(all_answers(:,10) - all_answers(:,11));     % Absolute error
    all_answers(:,13) = all_answers(:,10) - all_answers(:,11);          % Actual error
    all_answers(:,16) = abs(all_answers(:,10) - all_answers(:,15));     % Absolute error - path distance
    disp("Absolute scaled errors within between")
    [~, p] = ttest2(all_answers(locs_within, 12), all_answers(locs_between, 12));
    disp([nanmean(all_answers(locs_within, 12)) nanmean(all_answers(locs_between, 12)) p]);
    % Errors - path distance
    
    
    all_scaled_est_dists_within_segment(f) = nanmean(all_answers(locs_within, 10));
    all_scaled_est_dists_between_segments(f) = nanmean(all_answers(locs_between, 10));
    all_scaled_real_dists_within_segment(f) = nanmean(all_answers(locs_within, 11));
    all_scaled_real_dists_between_segments(f) = nanmean(all_answers(locs_between, 11));
    all_abs_error_within_segment(f) = nanmean(all_answers(locs_within, 12));
    all_abs_error_between_segments(f) = nanmean(all_answers(locs_between, 12));
    all_abs_error_overall(f) = nanmean(all_answers(all_relevant_locations, 12));
    all_error_within_segment(f) = nanmean(all_answers(locs_within, 13));
    all_error_between_segments(f) = nanmean(all_answers(locs_between, 13));
    all_abs_path_error_within_segment(f) = nanmean(all_answers(locs_within, 16));
    all_abs_path_error_between_segments(f) = nanmean(all_answers(locs_between, 16));
    all_abs_path_error_overall(f) = nanmean(all_answers(all_relevant_locations, 16));
    all_corr_real_est_overall(f) = actual_corr;
    all_corr_real_est_within(f) = corr(all_answers(locs_within,10),all_answers(locs_within,11));
    all_corr_real_est_between(f) = corr(all_answers(locs_between,10),all_answers(locs_between,11));
    all_corr_path_est_overall(f) = actual_path_corr;
    all_corr_path_est_within(f) = corr(all_answers(locs_within,10),all_answers(locs_within,15));
    all_corr_path_est_between(f) = corr(all_answers(locs_between,10),all_answers(locs_between,15));
    all_corr_overlay_est_overall(f) = actual_overlay_corr;
    all_corr_overlay_est_within(f) = corr(all_answers(locs_within,10),all_answers(locs_within,18));
    all_corr_overlay_est_between(f) = corr(all_answers(locs_between,10),all_answers(locs_between,18));
    all_est_dists_within_segment(f) = nanmean(all_answers(locs_within, 7));
    all_est_dists_between_segments(f) = nanmean(all_answers(locs_between, 7));
    
    [~, p] = ttest2(all_answers(locs_within, 10), all_answers(locs_between, 10));
    all_p_within_between_est_dists(f) = p;
    
    [~, p] = ttest2(all_answers(locs_within, 13), all_answers(locs_between, 13));
    all_p_within_between_average_error(f) = p;
    
    [~, p] = ttest2(all_answers(locs_within, 12), all_answers(locs_between, 12));
    all_p_within_between_average_abs_error(f) = p;
    
    [~, p] = ttest2(all_answers(locs_within, 16), all_answers(locs_between, 16));
    all_p_within_between_average_path_abs_error(f) = p;
    
    all_std_est_dist_within(f) = std(all_answers(locs_within, 10));
    all_std_est_dist_between(f) = std(all_answers(locs_between, 10));
    
    %     %% Creating a dissimilarity matrix and doing MDS - using non-scaled distances
    %     est_dists_mat = nan(16);
    %     for i=1:size(est_dists_mat, 1), est_dists_mat(i,i) = 0; end
    %     for i = 1:num_trials
    %         current_object1 = all_answers(i, 3) + 1;
    %         current_object2 = all_answers(i, 4) + 1;
    %         est_dists_mat(current_object1, current_object2) = all_answers(i, 7);
    %         est_dists_mat(current_object2, current_object1) = all_answers(i, 7);
    %     end
    %
    %     opts = statset('Display','final', 'MaxIter', 400);
    %     Y = mdscale(est_dists_mat, 2, 'Start', 'random', 'Replicates',20, 'Options', opts);
    %     all_mds{f} = Y;
    
    %
    %     Y = all_mds{1}; figure;scatter(Y(:,1), Y(:,2)); text(Y(:,1),Y(:,2),labels, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    
    %% RT priming - same segment comparison
    RTs_same_within_segment_answer = []; RTs_different_within_segment_answer = [];
    RTs_same_within_orth_segments_answer = []; RTs_different_within_orth_segments_answer = [];
    for i = 2:num_trials
        if all_answers(i,8) == 1 && all_answers(i-1,8) == 1    % if within segment, and the previous comparison was also within segment
            if all_answers(i,1) == all_answers(i-1,1) || all_answers(i,1) == all_answers(i-1,2)     % if the current and previous segments are the same
                RTs_same_within_segment_answer(end+1) = all_answers(i,6);
            else    % if the previous segment was different than the current one
                RTs_different_within_segment_answer(end+1) = all_answers(i,6);
            end
        elseif all_answers(i,8) == -1 && all_answers(i-1,8) == -1    % if between segments, and the previous comparison was also between segments
            if all_answers(i,1) == all_answers(i-1,1) || all_answers(i,1) == all_answers(i-1,2)     % if the current and previous comparisons are the same
                RTs_same_within_orth_segments_answer(end+1) = all_answers(i,6);
            else    % if the previous segment was different than the current one
                RTs_different_within_orth_segments_answer(end+1) = all_answers(i,6);
            end
        end
    end
    
    disp("RT priming - same different within segment comparison")
    if ~isempty(RTs_same_within_segment_answer) && ~isempty(RTs_different_within_segment_answer)
        [~, p] = ttest2(RTs_same_within_segment_answer, RTs_different_within_segment_answer);
        disp([nanmean(RTs_same_within_segment_answer) nanmean(RTs_different_within_segment_answer) p]);
    else
        disp([nanmean(RTs_same_within_segment_answer) nanmean(RTs_different_within_segment_answer)]);
    end
    disp("RT priming - same different within orthogonal segment comparison")
    if ~isempty(RTs_same_within_orth_segments_answer) && ~isempty(RTs_different_within_orth_segments_answer)
        [~, p] = ttest2(RTs_same_within_orth_segments_answer, RTs_different_within_orth_segments_answer);
        disp([nanmean(RTs_same_within_orth_segments_answer) nanmean(RTs_different_within_orth_segments_answer) p]);
    else
        disp([nanmean(RTs_same_within_orth_segments_answer) nanmean(RTs_different_within_orth_segments_answer)]);
    end
    
    all_RTs_same_within_segment_answer(f) = nanmean(RTs_same_within_segment_answer);
    all_RTs_different_within_segment_answer(f) = nanmean(RTs_different_within_segment_answer);
    all_RTs_same_within_orth_segments_answer(f) = nanmean(RTs_same_within_orth_segments_answer);
    all_RTs_different_within_orth_segments_answer(f) = nanmean(RTs_different_within_orth_segments_answer);
    
    % Calculating task overall time
    task_start_time = find(strcmp(results(:,1), 'L key pressed')); task_start_time = results{task_start_time(1), 12};
    task_end_time = find(strcmp(results(:,1), 'Finish message displayed'));
    if ~isempty(task_end_time)
        task_end_time = results{task_end_time(1), 12};
        all_task_overall_times(f) = (task_end_time - task_start_time) / 60;
    else
        all_task_overall_times(f) = nan;
    end
    
    answers_all_subjects{f} = all_answers;
    all_num_answers(f) = sum(~isnan(all_answers(:,7)));
    
end


% disp("All RTs within between")
% a = cell2mat(all_RTs_within_segment'); b = cell2mat(all_RTs_between_segments');
% [~, p] = ttest2(a, b);
% disp([nanmean(a) nanmean(b) p]);
%


%% Create subject-specific estimated distance matrices
% Answers - 3 object1, 4 object2, 7 input
for f = 1:length(filenames)
    dist_matrix = nan(16);
    for i = 1:16, dist_matrix(i,i) = 0; end     % Zeroing the diagonal
    for i = 1:size(answers_all_subjects{f}, 1)
        curr_object1 = answers_all_subjects{f}(i,3) + 1;
        curr_object2 = answers_all_subjects{f}(i,4) + 1;
        dist_matrix(curr_object1, curr_object2) = answers_all_subjects{f}(i,7);
        dist_matrix(curr_object2, curr_object1) = answers_all_subjects{f}(i,7);
    end
    
    % Saving the matrix
    current_subject_JRD_analysis_dir = fullfile(fullfile('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data', subjects_file{3+f, 29}), '/Analysis/Analysis_JRD');
    save(fullfile(current_subject_JRD_analysis_dir, 'estimated_distances_mat.mat'), 'dist_matrix');
end

