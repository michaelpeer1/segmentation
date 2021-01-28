% Distance comparison analysis

% Pilot data
% filenames = {};
% filenames{1} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1349.csv';  % MP
% filenames{2} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1716.csv';  % RM
% filenames{3} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1810.csv';  % ZL


% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v1_river\Subjects_file.xlsx');
% num_subjects = max(find(~cellfun(@isnan, subjects_file(4:end,1))));
[~, ~, subjects_file] = xlsread('/Users/mpeer/Desktop/Segmentation_data/fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v2_wall\Subjects_file_v2_wall.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_behavioral_pilot\Subjects_file_fMRI_behav_pilot.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_experiment\Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\5_behavioral_wall_new\Subjects_file_v5_wall_new.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\6_separate_learning\Subjects_file_v6_separate_learning.xlsx');
num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
filenames = {};
for i=1:num_subjects
    %     if ~isnan(subjects_file{3+i, 14})       % Behavioral version of the task
    %         filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 14}) '.csv'];
    if ~isnan(subjects_file{3+i, 14})       % fMRI version of the task
        %         filenames{end+1} = [fullfile(subjects_file{3+i, 8}, subjects_file{3+i, 14}) '.csv'];
        filenames{end+1} = [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 14}) '.csv'];
    end
end



answers_all_subjects = cell(1, length(filenames));

all_RTs_within = zeros(length(filenames), 1);
all_RTs_between = zeros(length(filenames), 1);
all_iscorrect_within = zeros(length(filenames), 1);
all_iscorrect_between = zeros(length(filenames), 1);
all_RT_same_quadrant = zeros(length(filenames), 1);
all_RT_different_quadrants = zeros(length(filenames), 1);
all_RT_same_segment = zeros(length(filenames), 1);
all_RT_different_segments = zeros(length(filenames), 1);
all_p_RTs_within_between = zeros(length(filenames), 1);
all_success_overall = zeros(length(filenames), 1);
all_success_within = zeros(length(filenames), 1);
all_success_between = zeros(length(filenames), 1);
all_RT_same_orth_segment = zeros(length(filenames), 1);
all_RT_different_orth_segments = zeros(length(filenames), 1);
all_task_overall_times = zeros(length(filenames), 1);
all_path_success_overall = zeros(length(filenames), 1);
all_path_success_within = zeros(length(filenames), 1);
all_path_success_between = zeros(length(filenames), 1);


for f = 1:length(filenames)
    %% Reading the data
    disp(strcat("File ", num2str(f)));
    %     [~,~,results] = xlsread(filenames{f});
    fid = fopen(filenames{f});
    data = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:length(results(:))        % Changing strings to numbers
        if ~isnan(str2double(results{i}))
            results{i} = str2double(results{i});
        end
    end
    
    locations_trials = [];
    for i=1:size(results, 1)
        if strcmp(results{i,1}, 'Question changed')
            locations_trials(end+1) = i;
        end
    end
    num_trials = length(locations_trials);
    locations_trials(end+1) = size(results,1);
    
    
    all_answers = nan(num_trials, 9);   % 1.quadrant1, 2.quadrant2, 3.object1, 4.object2, 5.object3, 6.input time, 7.input (left/right), 8.within_between, 9.correct_answer, 10.is_correct, 11.correct_path_answer, 12.is_correct_path
    
    for i=1:num_trials
        % Entering question parameters into the list
        all_answers(i,1:5) = cell2mat(results(locations_trials(i), 6:10));
        
        % Getting the other parameters from the next lines, if existing
        for j = locations_trials(i):locations_trials(i+1)
            if strcmp(results{j,1}, 'Arrow key pressed')
                all_answers(i, 6) = results{j, 14} - results{locations_trials(i), 14};
                current_answer = results{j, 5};
                if strcmp(current_answer,'Left'), all_answers(i, 7) = 0;
                elseif strcmp(current_answer,'Right'), all_answers(i, 7) = 1;
                end
            end
        end
        
        % Figuring out if comparison was between segments or within
        % 0 1 and 1 0 and 2 3 and 3 2 are between segments, 1 2 and 2 1 and 3 0 and 0 3 are within
        q1 = all_answers(i,1); q2 = all_answers(i,2);
        s1 = (q1 == 1 || q1 == 2); s2 = (q2 == 1 || q2 == 2);       % Segment
        if s1 == s2
            all_answers(i,8) = 1;   % within segments
        else
            all_answers(i,8) = -1;   % between segments
        end
    end
    
    locs_within = find(all_answers(:,8) == 1);
    locs_between = find(all_answers(:,8) == -1);
    all_relevant_locations = [locs_within; locs_between];
    
    
    %% Response times
    % Removing outliers (3 SDs)
    all_current_RTs = all_answers(all_relevant_locations, 6);
    outlier_thresh_high = nanmean(all_current_RTs) + 3*nanstd(all_current_RTs); outlier_thresh_low = nanmean(all_current_RTs) - 3*nanstd(all_current_RTs);
    current_RTs_within = all_answers(locs_within, 6); current_RTs_between = all_answers(locs_between, 6);
    current_RTs_within = current_RTs_within(current_RTs_within < outlier_thresh_high & current_RTs_within > outlier_thresh_low);
    current_RTs_between = current_RTs_between(current_RTs_between < outlier_thresh_high & current_RTs_between > outlier_thresh_low);
    
    % Difference between within / between segment comparisons in RT
    disp("RT within between")
    [~, p] = ttest2(current_RTs_within, current_RTs_between);
    disp([nanmean(current_RTs_within) nanmean(current_RTs_between) p])
    
    all_RTs_within(f) = nanmean(current_RTs_within);
    all_RTs_between(f) = nanmean(current_RTs_between);
    all_p_RTs_within_between(f) = p;
    
    
    %% Distances and errors
    %     locations_all_objects = [16 35; 32 63; 65 45; 41 17; ...
    %         35 -16; 63 -32; 45 -65; 17 -41; ...
    %         -16 -35; -32 -63; -65 -45; -41 -17; ...
    %         -35 16; -63 32; -45 65; -17 41];
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
    
    for i=1:num_trials
        % Calculating the correct answer
        real_dist_left = all_real_distances(all_answers(i,3)+1, all_answers(i,4)+1);
        real_dist_right = all_real_distances(all_answers(i,3)+1, all_answers(i,5)+1);
        if real_dist_left < real_dist_right
            all_answers(i,9) = 0;
        else
            all_answers(i,9) = 1;
        end
        all_answers(i,10) = (all_answers(i,7) == all_answers(i,9));     % is correct
        
        % Calculating the correct answer by path distance
        path_dist_left = all_path_distances(all_answers(i,3)+1, all_answers(i,4)+1);
        path_dist_right = all_path_distances(all_answers(i,3)+1, all_answers(i,5)+1);
        if path_dist_left < path_dist_right
            all_answers(i,11) = 0;
        else
            all_answers(i,11) = 1;
        end
        all_answers(i,12) = (all_answers(i,7) == all_answers(i,11));     % is correct
        
    end
    actual_success_rate = sum(all_answers(:, 10)) / num_trials;
    actual_path_success_rate = sum(all_answers(:, 12)) / num_trials;
    
    % Chance success level
    num_iters = 1000;
    rand_s = zeros(1,num_iters);
    for i=1:num_iters
        current_random_vector = randi(2,num_trials,1);
        rand_s(i) = sum(current_random_vector == all_answers(:,9)) / num_trials;
    end
    p_actual = sum(rand_s > actual_success_rate) / num_iters;
    
    disp("Success rate overall, significance vs. chance level")
    disp([actual_success_rate p_actual])
    disp("Success rate within between")
    disp([sum(all_answers(locs_within, 10))/length(locs_within) sum(all_answers(locs_between, 10))/length(locs_between)])
    
    all_success_overall(f) = actual_success_rate;
    all_success_within(f) = sum(all_answers(locs_within, 10))/length(locs_within);
    all_success_between(f) = sum(all_answers(locs_between, 10))/length(locs_between);
    
    all_path_success_overall(f) = actual_path_success_rate;
    all_path_success_within(f) = sum(all_answers(locs_within, 12))/length(locs_within);
    all_path_success_between(f) = sum(all_answers(locs_between, 12))/length(locs_between);
    
    
    %% RT priming - quadrant
    RTs_same_quadrant_answer = []; RTs_different_quadrant_answer = [];
    for i = 2:num_trials
        if all_answers(i,1) == all_answers(i-1,1)
            RTs_same_quadrant_answer(end+1) = all_answers(i,6);
        else
            RTs_different_quadrant_answer(end+1) = all_answers(i,6);
        end
    end
    
    disp("RT priming starting quadrant within between")
    [~, p] = ttest2(RTs_same_quadrant_answer, RTs_different_quadrant_answer);
    disp([mean(RTs_same_quadrant_answer) mean(RTs_different_quadrant_answer) p])
    
    
    % RT priming - segment
    RTs_same_segment_answer = []; RTs_different_segment_answer = [];
    RTs_same_orth_segment_answer = []; RTs_different_orth_segment_answer = [];
    for i = 2:num_trials
        q1a = all_answers(i,1); q1b = all_answers(i-1,1);
        s1a = (q1a == 1 || q1a == 2); s1b = (q1b == 1 || q1b == 2);       % Segment
        % 0 1 and 1 0 and 2 3 and 3 2 are between segments, 1 2 and 2 1 and 3 0 and 0 3 are within
        if s1a == s1b
            RTs_same_segment_answer(end+1) = all_answers(i,6);
        else
            RTs_different_segment_answer(end+1) = all_answers(i,6);
        end
        
        orth_s1a = (q1a == 1 || q1a == 0); orth_s1b = (q1b == 1 || q1b == 0);       % Orthogonal segment
        if orth_s1a == orth_s1b
            RTs_same_orth_segment_answer(end+1) = all_answers(i,6);
        else
            RTs_different_orth_segment_answer(end+1) = all_answers(i,6);
        end
    end
    
    disp("RT priming starting segment within between")
    [~, p] = ttest2(RTs_same_segment_answer, RTs_different_segment_answer);
    disp([mean(RTs_same_segment_answer) mean(RTs_different_segment_answer) p])
    
    
    all_RT_same_quadrant(f) = nanmean(RTs_same_quadrant_answer);
    all_RT_different_quadrants(f) = nanmean(RTs_different_quadrant_answer);
    all_RT_same_segment(f) = nanmean(RTs_same_segment_answer);
    all_RT_different_segments(f) = nanmean(RTs_different_segment_answer);
    all_RT_same_orth_segment(f) = nanmean(RTs_same_orth_segment_answer);
    all_RT_different_orth_segments(f) = nanmean(RTs_different_orth_segment_answer);
    
    % Calculating task overall time
    task_start_time = find(strcmp(results(:,1), 'Enter key pressed')); task_start_time = results{task_start_time(1), 14};
    task_end_time = find(strcmp(results(:,1), 'Finish message displayed')); task_end_time = results{task_end_time(1), 14};
    all_task_overall_times(f) = (task_end_time - task_start_time) / 60;
    
    answers_all_subjects{f} = all_answers;
end


% disp("All RTs within between")
% a = cell2mat(all_RTs_within'); b = cell2mat(all_RTs_between');
% [~, p] = ttest2(a, b);
% disp([nanmean(a) nanmean(b) p]);
%
% disp("All success within between")
% c = cell2mat(all_iscorrect_within'); d = cell2mat(all_iscorrect_between');
% [~, p] = ttest2(c, d);
% disp([nanmean(c) nanmean(d) p]);
%
% disp("All RT priming starting quadrant same different")
% e = cell2mat(all_RT_same_quadrant); f = cell2mat(all_RT_different_quadrants);
% [~, p] = ttest2(e, f);
% disp([nanmean(e) nanmean(f) p]);
%
% disp("All RT priming starting segment same different")
% g = cell2mat(all_RT_same_segment); h = cell2mat(all_RT_different_segments);
% [~, p] = ttest2(g, h);
% disp([nanmean(g) nanmean(h) p]);
%

