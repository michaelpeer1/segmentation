% Analyzing the free recall stages

% Pilot data analysis - old
% filenames = {};
% filenames{1} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1142.csv';     % MP1
% filenames{2} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1328.csv';     % MP2
% filenames{3} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1546.csv';     % RM1
% filenames{4} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1628.csv';     % RM2
% filenames{5} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1734(1)_corrected.csv';     % ZL1
% filenames{6} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1751_corrected.csv';     % ZL2
% filenames{7} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1243.csv';     % AT1
% filenames{8} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1303.csv';     % AT2

% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v1_river\Subjects_file.xlsx');
% num_subjects = max(find(~cellfun(@isnan, subjects_file(4:end,1))));
[~, ~, subjects_file] = xlsread('/Users/mpeer/Desktop/Segmentation_data/fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v2_wall\Subjects_file_v2_wall.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_behavioral_pilot\Subjects_file_fMRI_behav_pilot.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_experiment\Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\5_behavioral_wall_new\Subjects_file_v5_wall_new.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\6_separate_learning\Subjects_file_v6_separate_learning.xlsx');
num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
filenames = {}; learning_stage_filenames = {};
for i=1:num_subjects
%     if ~isnan(subjects_file{3+i, 10})       % Behavioral version of the task
%         filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 10}) '.csv'];
%         filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 12}) '.csv'];
%         learning_stage_filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 9}) '.csv'];
    if ~isnan(subjects_file{3+i, 12})       % fMRI version of the task
%         filenames{end+1} = [fullfile(subjects_file{3+i, 8}, subjects_file{3+i, 12}) '.csv'];
        filenames{end+1} = [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 12}) '.csv'];
%         learning_stage_filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 9}) '.csv'];      % to obtain the objects order
%         learning_stage_filenames{end+1} = [fullfile(subjects_file{3+i, 8}, subjects_file{3+i, 10}) '.csv'];      % to obtain the objects order
        learning_stage_filenames{end+1} = [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 10}) '.csv'];      % to obtain the objects order
    end
end

answers_all_subjects = cell(1, length(filenames));
all_RTs_within_segment_start = zeros(length(filenames), 1);
all_RTs_between_segments_start = zeros(length(filenames), 1);
all_RTs_within_segment_answer = zeros(length(filenames), 1);
all_RTs_between_segments_answer = zeros(length(filenames), 1);
all_RTs_within_quadrant_start = zeros(length(filenames), 1);
all_RTs_between_quadrants_start = zeros(length(filenames), 1);
all_RTs_within_quadrant_answer = zeros(length(filenames), 1);
all_RTs_between_quadrants_answer = zeros(length(filenames), 1);
all_RTs_within_orth_segment_start = zeros(length(filenames), 1);
all_RTs_between_orth_segments_start = zeros(length(filenames), 1);
all_RTs_within_orth_segment_answer = zeros(length(filenames), 1);
all_RTs_between_orth_segments_answer = zeros(length(filenames), 1);
all_task_times = zeros(length(filenames), 1);
all_num_transitions_within_segment = zeros(length(filenames), 1);
all_num_transitions_between_segments = zeros(length(filenames), 1);
all_num_transitions_within_orth_segment = zeros(length(filenames), 1);
all_num_transitions_between_orth_segments = zeros(length(filenames), 1);
all_num_transitions_within_quadrant = zeros(length(filenames), 1);
all_num_transitions_between_quadrants = zeros(length(filenames), 1);
all_num_transitions_within_segment_adj_quadrant = zeros(length(filenames), 1);
all_num_transitions_within_orth_segment_adj_quadrant = zeros(length(filenames), 1);
all_mean_transitions_distances = zeros(length(filenames), 1);
all_mean_path_transitions_distances = zeros(length(filenames), 1);
all_schema_preserving_transitions_percentage = zeros(length(filenames), 1);
all_overlay_preserving_transitions_percentage = zeros(length(filenames), 1);
all_num_identified_words = zeros(length(filenames), 1);

for f = 1:length(filenames)
    disp(strcat("File ", num2str(f)));
    
    %% Reading the data
    %     [~,~,results] = xlsread(filenames{f});
    fid = fopen(filenames{f});
    data = textscan(fid, '%s %s %s %s %s %s %s', 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:size(results, 1), results{i,6} = str2double(results{i,6}); results{i,7} = str2double(results{i,7}); end   % Changing strings to numbers

    locations_words_entered = [];
    for i=1:size(results, 1)
        if strcmp(results{i,1}, 'Enter key pressed; new word is')
            locations_words_entered(end+1) = i;
        end
    end
    num_words = length(locations_words_entered);
    
    %     object_names = {'cone','lamp','motorcycle','guitar',...
    %         'hydrant','pyramid','book','statue',...
    %         'bell','table','chest','umbrella',...
    %         'tree','ship','snowman','chair'};
    
    % Getting the objects order (could be different from default by randomization)
%     [~, ~, learning_stage_results] = xlsread(learning_stage_filenames{ceil(f/2)});        % for behavioral version when there are two files per subject
%     temp_loc_randomization = find(strcmp(learning_stage_results(:,1), 'Object randomization seed'));
%     object_names = lower(learning_stage_results((temp_loc_randomization+1):(temp_loc_randomization+16), 8));
    fid = fopen(learning_stage_filenames{f});
    learning_stage_results = textscan(fid, '%s %s %s %s %s %s %s %s %s', 'delimiter', ',');
    fclose(fid);
    temp_loc_randomization = find(strcmp(learning_stage_results{:,1}, 'Object randomization seed'));
    object_names = lower(learning_stage_results{8}((temp_loc_randomization+1):(temp_loc_randomization+16)));
    
%     % accounting for words that were deleted by the "delete last word" button
    entered_words_indices = cell2mat(results(locations_words_entered, 6));
%     locations_temp = zeros(1, length(object_names));
%     for i = 1:length(object_names)
%         locations_temp(i) = max(find(entered_words_indices == i));  % Choosing the last location where the word was entered, if entered multiple times
%     end
%     locations_words_entered = locations_words_entered(locations_temp);
    num_words = length(locations_words_entered);
    
    all_answers = cell(num_words, 6);   % word, RT to start, RT to answer, num_object, quadrant, segment, orthogonal segment
    
    for i = 1:num_words
        all_answers{i,1} = results{locations_words_entered(i), 5};      % word
        
        if sum(entered_words_indices == i) == 1 && sum(entered_words_indices == (i - 1)) == 1    % Response times - ignoring words that were entered more than once, as it's unclear what their RT is, and the words after them
            all_answers{i,2} = results{locations_words_entered(i)-1, 7} - results{locations_words_entered(i)-2, 7};    % RT to word start
            all_answers{i,3} = results{locations_words_entered(i), 7} - results{locations_words_entered(i)-2, 7};      % RT to word enter
        else
            all_answers{i,2} = nan; all_answers{i,3} = nan;
        end
        
        all_answers{i,4} = nan; all_answers{i,5} = nan; all_answers{i,6} = nan;
        for j=1:length(object_names)
            if strcmp(lower(all_answers{i,1}), object_names{j})
                object = j - 1;
                quadrant = floor(object / 4);
                segment = (quadrant == 1 || quadrant == 2);
                orth_segment = (quadrant == 2 || quadrant == 3);    % orthogonal separation to halves from segmentation
                
                all_answers{i,4} = object;
                all_answers{i,5} = quadrant;
                all_answers{i,6} = segment;
                all_answers{i,7} = orth_segment;
            end
        end
    end
    
    disp(strcat("Number of identified words: ", num2str(sum(~isnan([all_answers{:,4}])))));
    all_num_identified_words(f) = sum(~isnan([all_answers{:,4}]));
    
    
    %% Transition probabilities
    num_transitions_within_segment = 0; num_transitions_between_segments = 0;
    num_transitions_within_quadrant = 0; num_transitions_between_quadrants = 0;
    num_transitions_within_orth_segment = 0; num_transitions_between_orth_segments = 0;
    num_transitions_within_segment_adj_quadrant = 0; num_transitions_within_orth_segment_adj_quadrant = 0; 
    for i=2:num_words
        if ~isnan(all_answers{i,4})
            if all_answers{i,6} == all_answers{i-1, 6}  % same segment
                num_transitions_within_segment = num_transitions_within_segment + 1;
                if all_answers{i,5} ~= all_answers{i-1, 5}      % Same segment adjacent quadrant
                    num_transitions_within_segment_adj_quadrant = num_transitions_within_segment_adj_quadrant + 1;
                end
            else    % different segment
                num_transitions_between_segments = num_transitions_between_segments + 1;
            end
            
            if all_answers{i,5} == all_answers{i-1, 5}  % same quadrant
                num_transitions_within_quadrant = num_transitions_within_quadrant + 1;
            else    % different quadrant
                num_transitions_between_quadrants = num_transitions_between_quadrants + 1;
            end
            
            if all_answers{i,7} == all_answers{i-1, 7}  % same orthogonal segment
                num_transitions_within_orth_segment = num_transitions_within_orth_segment + 1;
                if all_answers{i,5} ~= all_answers{i-1, 5}      % Same segment adjacent quadrant
                    num_transitions_within_orth_segment_adj_quadrant = num_transitions_within_orth_segment_adj_quadrant + 1;
                end
            else    % different orthogonal segment
                num_transitions_between_orth_segments = num_transitions_between_orth_segments + 1;
            end
        end
    end
    
    disp("Num transitions within between segments")
    disp([num_transitions_within_segment num_transitions_between_segments])
    disp("Num transitions within between orthogonal segments")
    disp([num_transitions_within_orth_segment num_transitions_between_orth_segments])
    disp("Num transitions within between quadrants")
    disp([num_transitions_within_quadrant num_transitions_between_quadrants])
    
    all_num_transitions_within_segment(f) = num_transitions_within_segment;
    all_num_transitions_between_segments(f) = num_transitions_between_segments;
    all_num_transitions_within_orth_segment(f) = num_transitions_within_orth_segment;
    all_num_transitions_between_orth_segments(f) = num_transitions_between_orth_segments;
    all_num_transitions_within_quadrant(f) = num_transitions_within_quadrant;
    all_num_transitions_between_quadrants(f) = num_transitions_between_quadrants;
    all_num_transitions_within_segment_adj_quadrant(f) = num_transitions_within_segment_adj_quadrant;
    all_num_transitions_within_orth_segment_adj_quadrant(f) = num_transitions_within_orth_segment_adj_quadrant;

    
    %% Transition probability correlation to distance
%     locations_all_objects = [16 35; 32 63; 65 45; 41 17; ...
%         35 -16; 63 -32; 45 -65; 17 -41; ...
%         -16 -35; -32 -63; -65 -45; -41 -17; ...
%         -35 16; -63 32; -45 65; -17 41];
    % Calculating Euclidean distances between objects
    locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
        42 -16; 58 -50; 27 -57; 17 -25; ...
        -16 -42; -50 -58; -57 -27; -25 -17; ...
        -42 16; -58 50; -27 57; -17 25];
    all_real_distances = pdist2(locations_all_objects, locations_all_objects);
    
    % Calculating path distances
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
    
    transitions_distances = []; transitions_path_distances = [];
    for i = 2:num_words
        o1 = all_answers{i, 4}; o2 = all_answers{i - 1, 4};
        if ~isnan(o1) && ~isnan(o2)
            transitions_distances(end+1) = all_real_distances(o1+1, o2+1);
            transitions_path_distances(end+1) = all_path_distances(o1+1, o2+1);
        end
    end
    num_transitions = length(transitions_distances);
    
%     % Creating a random distribution of transition chains
%     num_iters = 10000;
%     rand_average_distances = zeros(1, num_iters);
%     rand_average_path_distances = zeros(1, num_iters);
%     for i=1:num_iters
%         all_current_dists = zeros(1, num_transitions);
%         all_current_path_dists = zeros(1, num_transitions);
%         all_current_objects = 1:16;
%         all_current_objects = all_current_objects(randperm(length(all_current_objects)));
%         for j = 1:num_transitions
%             o1 = all_current_objects(j); o2 = all_current_objects(j+1);
%             all_current_dists(j) = all_real_distances(o1,o2);
%             all_current_path_dists(j) = all_path_distances(o1,o2);
%         end
%         rand_average_distances(i) = mean(all_current_dists);
%         rand_average_path_distances(i) = mean(all_current_path_dists);
%     end
%     p = sum(rand_average_distances <= mean(transitions_distances)) / num_iters;
    
    disp(strcat("mean distance per transition = ", num2str(mean(transitions_distances))));
%     disp(strcat("mean random distance per transition = ", num2str(mean(rand_average_distances))));
%     disp(strcat("p = ", num2str(p)));
    
    all_mean_transitions_distances(f) = mean(transitions_distances);
    all_mean_path_transitions_distances(f) = mean(transitions_path_distances);
    
    
    %% Schema preserving transitions
    schema_preserving_transitions = 0;
    for i = 2:num_words
        o1 = all_answers{i, 4}; o2 = all_answers{i - 1, 4};
        s1 = all_answers{i, 6}; s2 = all_answers{i - 1, 6};
        if ~isnan(o1) && ~isnan(o2)
            if s1 ~= s2 && mod(o1,8) == mod(o2,8)   % Schema preserving transition
                schema_preserving_transitions = schema_preserving_transitions + 1;
            end
        end
    end
    schema_preserving_transitions_percentage = schema_preserving_transitions / num_transitions_between_segments;
    
    % Overlay model preserving transitions
    overlay_transitions = [1 8; 2 5; 3 6; 4 7; 5 2; 6 3; 7 4; 8 1; 9 16; 10 13; 11 14; 12 15; 13 10; 14 11; 15 12; 16 9];
    overlay_preserving_transitions = 0;
    for i = 2:num_words
        o1 = all_answers{i, 4}; o2 = all_answers{i - 1, 4};
        s1 = all_answers{i, 6}; s2 = all_answers{i - 1, 6};
        if ~isnan(o1) && ~isnan(o2)
            if s1 ~= s2 && overlay_transitions(o1+1, 2) == o2+1   % Overlay model preserving transition
                overlay_preserving_transitions = overlay_preserving_transitions + 1;
            end
        end
    end
    overlay_preserving_transitions_percentage = overlay_preserving_transitions / num_transitions_between_segments;

    %     % Creating a random distribution of transition chains
%     num_iters = 10000;
%     rand_schema_preserving_transitions_percentage = zeros(1, num_iters);
%     for i=1:num_iters
%         current_schema_preserving_transitions = 0;
%         current_between_segments_transitions = 0;
%         all_current_objects = 1:16;
%         all_current_objects = all_current_objects(randperm(length(all_current_objects)));
%         for j = 1:num_transitions
%             o1 = all_current_objects(j) - 1; o2 = all_current_objects(j+1) - 1; % objects
%             q1 = floor(o1/4); q2 = floor(o2/4);   % quadrant
%             s1 = (q1 == 1 || q1 == 2); s2 = (q2 == 1 || q2 == 2);       % Segment
%             if s1~=s2
%                 current_between_segments_transitions = current_between_segments_transitions + 1;
%                 if mod(o1,8) == mod(o2, 8)
%                     current_schema_preserving_transitions = current_schema_preserving_transitions + 1;
%                 end
%             end
%         end
%         rand_schema_preserving_transitions_percentage(i) = current_schema_preserving_transitions / current_between_segments_transitions;
%     end
%     p = sum(schema_preserving_transitions_percentage <= rand_schema_preserving_transitions_percentage) / num_iters;
    
    disp(strcat("percent schema preserving transitions (from between-segment) = ", num2str(schema_preserving_transitions_percentage)));
%     disp(strcat("mean random percent schema preserving transitions = ", num2str(mean(rand_schema_preserving_transitions_percentage))));
%     disp(strcat("p = ", num2str(p)));
    all_schema_preserving_transitions_percentage(f) = schema_preserving_transitions_percentage;
    all_overlay_preserving_transitions_percentage(f) = overlay_preserving_transitions_percentage;
    
    
    
    %% Response times
    
    % Removing RT outliers - setting thersholds (3 SDs)
    all_current_RTs_start = cell2mat(all_answers(:, 2));
    outlier_thresh_high_start = nanmean(all_current_RTs_start) + 3*nanstd(all_current_RTs_start); outlier_thresh_low_start = nanmean(all_current_RTs_start) - 3*nanstd(all_current_RTs_start);
    all_current_RTs_answer = cell2mat(all_answers(:, 3));
    outlier_thresh_high_answer = nanmean(all_current_RTs_answer) + 3*nanstd(all_current_RTs_answer); outlier_thresh_low_answer = nanmean(all_current_RTs_answer) - 3*nanstd(all_current_RTs_answer);
    
    for i=2:num_words
        if all_answers{i,2}<outlier_thresh_low_start || all_answers{i,2}>outlier_thresh_high_start || all_answers{i,3}<outlier_thresh_low_answer || all_answers{i,3}>outlier_thresh_high_answer || isnan(all_answers{i,4})
            all_answers{i,2} = nan; all_answers{i,3} = nan;
        end
    end
    
    % Calculating RTs
    RTs_within_quadrant_start = []; RTs_between_quadrants_start = [];
    RTs_within_quadrant_answer = []; RTs_between_quadrants_answer = [];
    RTs_within_segment_start = []; RTs_between_segments_start = [];
    RTs_within_segment_answer = []; RTs_between_segments_answer = [];
    RTs_within_orth_segment_start = []; RTs_between_orth_segments_start = [];
    RTs_within_orth_segment_answer = []; RTs_between_orth_segments_answer = [];
    for i=2:num_words
        if ~isnan(all_answers{i,2})
            if all_answers{i,5} == all_answers{i-1, 5}  % same quadrant
                RTs_within_quadrant_start(end+1) = all_answers{i,2}; RTs_within_quadrant_answer(end+1) = all_answers{i,3};
            else    % different quadrant
                RTs_between_quadrants_start(end+1) = all_answers{i,2}; RTs_between_quadrants_answer(end+1) = all_answers{i,3};
            end
            
            if all_answers{i,6} == all_answers{i-1, 6}  % same segment
                RTs_within_segment_start(end+1) = all_answers{i,2}; RTs_within_segment_answer(end+1) = all_answers{i,3};
            else    % different segment
                RTs_between_segments_start(end+1) = all_answers{i,2}; RTs_between_segments_answer(end+1) = all_answers{i,3};
            end
            
            if all_answers{i,7} == all_answers{i-1, 7}  % same orthogonal segment
                RTs_within_orth_segment_start(end+1) = all_answers{i,2}; RTs_within_orth_segment_answer(end+1) = all_answers{i,3};
            else    % different orthogonal segment
                RTs_between_orth_segments_start(end+1) = all_answers{i,2}; RTs_between_orth_segments_answer(end+1) = all_answers{i,3};
            end
        end
    end
    
    disp("RT transitions within between segments")
    %     [~, p] = ttest2(RTs_within_segment_start, RTs_between_segments_start);
    %     disp([mean(RTs_within_segment_start) mean(RTs_between_segments_start) p])
    [~, p] = ttest2(RTs_within_segment_answer, RTs_between_segments_answer);
    disp([mean(RTs_within_segment_answer) mean(RTs_between_segments_answer) p])
    
    
    all_RTs_within_quadrant_start(f) = nanmean(RTs_within_quadrant_start);
    all_RTs_between_quadrants_start(f) = nanmean(RTs_between_quadrants_start);
    all_RTs_within_quadrant_answer(f) = nanmean(RTs_within_quadrant_answer);
    all_RTs_between_quadrants_answer(f) = nanmean(RTs_between_quadrants_answer);
    all_RTs_within_segment_start(f) = nanmean(RTs_within_segment_start);
    all_RTs_between_segments_start(f) = nanmean(RTs_between_segments_start);
    all_RTs_within_segment_answer(f) = nanmean(RTs_within_segment_answer);
    all_RTs_between_segments_answer(f) = nanmean(RTs_between_segments_answer);
    all_RTs_within_orth_segment_start(f) = nanmean(RTs_within_orth_segment_start);
    all_RTs_between_orth_segments_start(f) = nanmean(RTs_between_orth_segments_start);
    all_RTs_within_orth_segment_answer(f) = nanmean(RTs_within_orth_segment_answer);
    all_RTs_between_orth_segments_answer(f) = nanmean(RTs_between_orth_segments_answer);
    
    
    task_start_time = results{find(strcmp(results(:,1), 'Instructions displayed')) + 1, 7};
    task_end_time = results{find(strcmp(results(:,1), 'Finish text displayed')), 7};
    all_task_times(f) = (task_end_time - task_start_time) / 60;
    
    answers_all_subjects{f} = all_answers;
end

