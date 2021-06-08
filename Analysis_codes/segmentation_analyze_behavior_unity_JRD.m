% JRD analysis

% Pilot data
% filenames = {};
% filenames{1} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1144.csv';      % MP
% filenames{2} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1548.csv';      % RM
% filenames{3} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1737.csv';      % RM
% filenames{4} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1244.csv';      % RM

% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v1_river\Subjects_file.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v1_river\Subjects_file.xlsx');
[~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Older versions\2_Behavioral_data_v1_river\Subjects_file.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v2_wall\Subjects_file_v2_wall.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\5_behavioral_wall_new\Subjects_file_v5_wall_new.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\6_separate_learning\Subjects_file_v6_separate_learning.xlsx');
% num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
filenames = {};
for i=1:num_subjects
    if ~isnan(subjects_file{3+i, 11})
        filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 11}) '.csv'];
    end
end

answers_all_subjects = cell(1, length(filenames));

all_RTs_within_started = zeros(length(filenames), 1);
all_RTs_between_started = zeros(length(filenames), 1);
all_RTs_within_entered = zeros(length(filenames), 1);
all_RTs_between_entered = zeros(length(filenames), 1);

all_RTs_primed_same_quadrant_start = zeros(length(filenames), 1);
all_RTs_primed_different_quadrant_start = zeros(length(filenames), 1);
all_RTs_primed_same_quadrant_enter = zeros(length(filenames), 1);
all_RTs_primed_different_quadrant_enter = zeros(length(filenames), 1);

all_RTs_primed_same_segment_start = zeros(length(filenames), 1);
all_RTs_primed_different_segment_start = zeros(length(filenames), 1);
all_RTs_primed_same_segment_enter = zeros(length(filenames), 1);
all_RTs_primed_different_segment_enter = zeros(length(filenames), 1);

all_angular_errors_within = zeros(length(filenames), 1);
all_angular_errors_between = zeros(length(filenames), 1);
all_task_overall_times = zeros(length(filenames), 1);
all_angular_errors_within_orth = zeros(length(filenames), 1);
all_angular_errors_between_orth = zeros(length(filenames), 1);
all_corr_superordinate = zeros(length(filenames), 1);
all_bias_within = zeros(length(filenames), 1);
all_bias_between = zeros(length(filenames), 1);
all_segment_direction_dev = zeros(length(filenames), 1);
all_segment_direction_dev_within = zeros(length(filenames), 1);
all_segment_direction_dev_between = zeros(length(filenames), 1);
p_within_subject = zeros(length(filenames), 1);
accuracy_overall = zeros(length(filenames), 1);
all_p_angular_error_within_between = zeros(length(filenames), 1);

all_angular_errors_overall = zeros(length(filenames), 1);
all_num_answers = zeros(length(filenames), 1);



for f = 1:length(filenames)
    %% Reading the data
    disp(strcat("File ", num2str(f)));
    [~,~,results] = xlsread(filenames{f});
    
    locations_trials = [];
    for i=1:size(results, 1)
        if strcmp(results{i,1}, 'Question changed')
            locations_trials(end+1) = i;
        end
    end
    num_trials = length(locations_trials);
    locations_trials(end+1) = size(results,1);
    
    all_answers = nan(num_trials, 11);   % 1.quadrant1, 2.quadrant2, 3.object1, 4.object2, 5.object3, ...
    % 6.time to rotation start, 7.time to final answer, 8.entered direction, 9.within/between, 10.true angle, ...
    % 11.angle difference (estimated vs. actual), 12.within/between_orthogonal_segment, ...
    % 13.true angle relative to environment, 14.estimated angle relative to environment, 15.vertical deviation, 16.segment direction deviation
    
    for i=1:num_trials
        % Entering question parameters into the list
        all_answers(i,1:5) = cell2mat(results(locations_trials(i), 6:10));
        
        % Getting the other parameters from the next lines, if existing
        for j = locations_trials(i):locations_trials(i+1)
            if strcmp(results{j,1}, 'Compass rotation started')
                all_answers(i, 6) = results{j, 14} - results{locations_trials(i), 14};
            elseif strcmp(results{j,1}, 'Compass direction entered')
                all_answers(i, 7) = results{j, 14} - results{locations_trials(i), 14};
                all_answers(i, 8) = results{j, 5};
            end
        end
        
        % Figuring out if comparison was between segments or within
        % 0 1 and 1 0 and 2 3 and 3 2 are between segments, 1 2 and 2 1 and 3 0 and 0 3 are within
        q1 = all_answers(i,1); q2 = all_answers(i,2);
        s1 = (q1 == 1 || q1 == 2); s2 = (q2 == 1 || q2 == 2);       % Segment
        %         if ((q1 == 1 && q2 == 2) || (q1 == 2 && q2 == 1) || (q1 == 3 && q2 == 0) || (q1 == 0 && q2 == 3))
        if s1 == s2, all_answers(i,9) = 1;  % within segment
        else, all_answers(i,9) = -1; end    % between segments
        
        orth_s1 = (q1 == 1 || q1 == 0); orth_s2 = (q2 == 1 || q2 == 0);       % Orthogonal segment direction
        if orth_s1 == orth_s2, all_answers(i,12) = 1;  % within orthogonal segment
        else, all_answers(i,12) = -1; end    % between orthogonal segments
        
    end
    
    locs_within = find(all_answers(:,9) == 1);
    locs_between = find(all_answers(:,9) == -1);
    all_relevant_locations = [locs_within; locs_between];
    
    locs_within_orth = find(all_answers(:,12) == 1);
    locs_between_orth = find(all_answers(:,12) == -1);
    
    % Calculating task overall time
    task_start_time = find(strcmp(results(:,1), 'L key pressed')); task_start_time = results{task_start_time(1), 14};
    task_end_time = find(strcmp(results(:,1), 'Finish message displayed')); task_end_time = results{task_end_time(1), 14};
    all_task_overall_times(f) = (task_end_time - task_start_time) / 60;
    
    
    %% Response times
    % Removing outliers (3 SDs)
    all_current_RTs_start = all_answers(all_relevant_locations, 6);
    outlier_thresh_high_start = nanmean(all_current_RTs_start) + 3*nanstd(all_current_RTs_start); outlier_thresh_low_start = nanmean(all_current_RTs_start) - 3*nanstd(all_current_RTs_start);
    %     current_RTs_within_start = all_answers(locs_within, 6); current_RTs_between_start = all_answers(locs_between, 6);
    %     current_RTs_within_start = current_RTs_within_start(current_RTs_within_start < outlier_thresh_high_start & current_RTs_within_start > outlier_thresh_low_start);
    %     current_RTs_between_start = current_RTs_between_start(current_RTs_between_start < outlier_thresh_high_start & current_RTs_between_start > outlier_thresh_low_start);
    
    all_current_RTs_answer = all_answers(all_relevant_locations, 7);
    outlier_thresh_high_answer = nanmean(all_current_RTs_answer) + 3*nanstd(all_current_RTs_answer); outlier_thresh_low_answer = nanmean(all_current_RTs_answer) - 3*nanstd(all_current_RTs_answer);
    %     current_RTs_within_answer = all_answers(locs_within, 7); current_RTs_between_answer = all_answers(locs_between, 7);
    %     current_RTs_within_answer = current_RTs_within_answer(current_RTs_within_answer < outlier_thresh_high_answer & current_RTs_within_answer > outlier_thresh_low_answer);
    %     current_RTs_between_answer = current_RTs_between_answer(current_RTs_between_answer < outlier_thresh_high_answer & current_RTs_between_answer > outlier_thresh_low_answer);
    
    for i = 1:num_trials
        if isnan(all_answers(i,6)) || all_answers(i,6)<outlier_thresh_low_start || all_answers(i,6)>outlier_thresh_high_start || all_answers(i,7)<outlier_thresh_low_answer || all_answers(i,7)>outlier_thresh_high_answer
            %             for q = 6:size(all_answers,2), all_answers(i,q) = nan; end
            for q = 6:7, all_answers(i,q) = nan; end
        end
    end
    current_RTs_within_start = all_answers(locs_within,6);
    current_RTs_between_start = all_answers(locs_between,6);
    current_RTs_within_answer = all_answers(locs_within,7);
    current_RTs_between_answer = all_answers(locs_between,7);
    % Difference between within / between segment comparisons in RT to rotation start / final answer
    disp("RT within between")
    %     [~, p] = ttest2(current_RTs_within_start, current_RTs_between_start);
    %     disp([nanmean(current_RTs_within_start) nanmean(current_RTs_between_start) p])
    [~, p] = ttest2(current_RTs_within_answer, current_RTs_between_answer);
    disp([nanmean(current_RTs_within_answer) nanmean(current_RTs_between_answer) p])
    
    
    all_RTs_within_started(f) = nanmean(current_RTs_within_start);
    all_RTs_between_started(f) = nanmean(current_RTs_between_start);
    all_RTs_within_entered(f) = nanmean(current_RTs_within_answer);
    all_RTs_between_entered(f) = nanmean(current_RTs_between_answer);
    
    
    %% RT priming - quadrant
    RTs_same_quadrant_start = []; RTs_different_quadrant_start = [];
    RTs_same_quadrant_answer = []; RTs_different_quadrant_answer = [];
    for i = 2:num_trials
        if all_answers(i,1) == all_answers(i-1,1)
            %             if (all_answers(i,6) > outlier_thresh_low_start && all_answers(i,6) < outlier_thresh_high_start)
            RTs_same_quadrant_start(end+1) = all_answers(i,6);
            %             end
            %             if (all_answers(i,7) > outlier_thresh_low_answer && all_answers(i,7) < outlier_thresh_high_answer)
            RTs_same_quadrant_answer(end+1) = all_answers(i,7);
            %             end
        else
            %             if (all_answers(i,6) > outlier_thresh_low_start && all_answers(i,6) < outlier_thresh_high_start)
            RTs_different_quadrant_start(end+1) = all_answers(i,6);
            %             end
            %             if (all_answers(i,7) > outlier_thresh_low_answer && all_answers(i,7) < outlier_thresh_high_answer)
            RTs_different_quadrant_answer(end+1) = all_answers(i,7);
            %             end
        end
    end
    
    disp("RT priming quadrant same different")
    %     [~, p] = ttest2(RTs_same_quadrant_start, RTs_different_quadrant_start);
    %     disp([nanmean(RTs_same_quadrant_start) nanmean(RTs_different_quadrant_start) p])
    [~, p] = ttest2(RTs_same_quadrant_answer, RTs_different_quadrant_answer);
    disp([nanmean(RTs_same_quadrant_answer) nanmean(RTs_different_quadrant_answer) p])
    
    
    % RT priming - segment
    RTs_same_segment_start = []; RTs_different_segment_start = [];
    RTs_same_segment_answer = []; RTs_different_segment_answer = [];
    for i = 2:num_trials
        q1a = all_answers(i,1); q1b = all_answers(i-1,1);       % quadrants
        s1a = (q1a == 1 || q1a == 2); s1b = (q1b == 1 || q1b == 2);       % Segments
        
        %         % 0 1 and 1 0 and 2 3 and 3 2 are between segments, 1 2 and 2 1 and 3 0 and 0 3 are within
        %         if ((q1a == 1 && q1b == 2) || (q1a == 2 && q1b == 1) || (q1a == 3 && q1b == 0) || (q1a == 0 && q1b == 3))
        if s1a == s1b
            %             if (all_answers(i,6) > outlier_thresh_low_start && all_answers(i,6) < outlier_thresh_high_start)
            RTs_same_segment_start(end+1) = all_answers(i,6);
            %             end
            %             if (all_answers(i,7) > outlier_thresh_low_answer && all_answers(i,7) < outlier_thresh_high_answer)
            RTs_same_segment_answer(end+1) = all_answers(i,7);
            %             end
        else
            %             if (all_answers(i,6) > outlier_thresh_low_start && all_answers(i,6) < outlier_thresh_high_start)
            RTs_different_segment_start(end+1) = all_answers(i,6);
            %             end
            %             if (all_answers(i,7) > outlier_thresh_low_answer && all_answers(i,7) < outlier_thresh_high_answer)
            RTs_different_segment_answer(end+1) = all_answers(i,7);
            %             end
        end
    end
    
    disp("RT priming segment same different")
    %     [~, p] = ttest2(RTs_same_segment_start, RTs_different_segment_start);
    %     disp([nanmean(RTs_same_segment_start) nanmean(RTs_different_segment_start) p])
    [~, p] = ttest2(RTs_same_segment_answer, RTs_different_segment_answer);
    disp([nanmean(RTs_same_segment_answer) nanmean(RTs_different_segment_answer) p])
    
    
    
    all_RTs_primed_same_quadrant_start(f) = nanmean(RTs_same_quadrant_start);
    all_RTs_primed_different_quadrant_start(f) = nanmean(RTs_different_quadrant_start);
    all_RTs_primed_same_quadrant_enter(f) = nanmean(RTs_same_quadrant_answer);
    all_RTs_primed_different_quadrant_enter(f) = nanmean(RTs_different_quadrant_answer);
    
    all_RTs_primed_same_segment_start(f) = nanmean(RTs_same_segment_start);
    all_RTs_primed_different_segment_start(f) = nanmean(RTs_different_segment_start);
    all_RTs_primed_same_segment_enter(f) = nanmean(RTs_same_segment_answer);
    all_RTs_primed_different_segment_enter(f) = nanmean(RTs_different_segment_answer);
    
    
    
    
    %% Angular errors
    %     locations_all_objects = [16 35; 32 63; 65 45; 41 17; ...
    %         35 -16; 63 -32; 45 -65; 17 -41; ...
    %         -16 -35; -32 -63; -65 -45; -41 -17; ...
    %         -35 16; -63 32; -45 65; -17 41];
    locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
        42 -16; 58 -50; 27 -57; 17 -25; ...
        -16 -42; -50 -58; -57 -27; -25 -17; ...
        -42 16; -58 50; -27 57; -17 25];
    
    % figure;scatter(locations_all_objects(:,1),locations_all_objects(:,2))
    
    for i=1:num_trials
        %         if ~isnan(all_answers(i,6))     % Ignoring trials with no rotation
        % computing the real angle between the vectors
        o1 = all_answers(i,3) + 1;
        o2 = all_answers(i,4) + 1;
        o3 = all_answers(i,5) + 1;
        
        v1 = locations_all_objects(o2, :) - locations_all_objects(o1, :);   % direction between object 1 and 2
        v2 = locations_all_objects(o3, :) - locations_all_objects(o1, :);   % direction between object 1 and 3
        
        %         angle = acos(v1*v2' / (norm(v1)*norm(v2)));
        %         angle = angle * 180 / pi;       % converting to radians
        %         angle = -atan2d(v1(1)*v2(2)-v2(1)*v1(2),v1(1)*v1(2)+v2(1)*v2(2));
        %         if angle < 0, angle = 360 + angle; end
        angle = atan2d(v1(2),v1(1)) - atan2d(v2(2),v2(1));
        if abs(angle) > 180, angle = angle - 360*sign(angle); end
        if angle < 0, angle = 360 + angle; end
        
        all_answers(i, 10) = angle;
        
        % Computing the difference between real and estimated angles
        angle_diff = all_answers(i, 10) - all_answers(i,8);
        angle_diff = abs(mod(angle_diff + 180, 360) - 180);     % Converting to the smallest angle of the possible two between the two angles
        all_answers(i, 11) = angle_diff;
        %         end
        
        % Calculating true and estimated angle relative to the environment
        new_angle_o13 = atan2d(v2(2),v2(1));     % The angle of the vector from object 1 to object 3
        if abs(new_angle_o13) > 180, new_angle_o13 = new_angle_o13 - 360*sign(new_angle_o13); end
        if new_angle_o13 < 0, new_angle_o13 = 360 + new_angle_o13; end
        % Calculating the estimated angle relative to the environment (estimated direction when taking into account viewing direction)
        % The answer that is saved by the unity task is the clockwise rotation - opposite than real angle
        new_angle_o12 = atan2d(v1(2),v1(1));
        if abs(new_angle_o12) > 180, new_angle_o12 = new_angle_o12 - 360*sign(new_angle_o12); end
        if new_angle_o12 < 0, new_angle_o12 = 360 + new_angle_o12; end
        new_est_angle = new_angle_o12 - all_answers(i,8);
        if abs(new_est_angle) > 180, new_est_angle = new_est_angle - 360*sign(new_est_angle); end
        if new_est_angle < 0, new_est_angle = 360 + new_est_angle; end
        all_answers(i, 13) = new_angle_o13;
        all_answers(i, 14) = new_est_angle;
        
        % Calculating vertical bias and horizontal bias of answers - overall (regardless of direction of segment)
        diff_angles = new_est_angle - new_angle_o13;
        if abs(diff_angles) > 180, diff_angles = diff_angles - 360*sign(diff_angles); end
        % Vertical deviations (between-segments direction) are scored as positive, horizontal scored as negative
        if new_angle_o13<=90
            deviation = diff_angles;
        elseif new_angle_o13>90 && new_angle_o13<=180
            deviation = -1 * diff_angles;
        elseif new_angle_o13>180 && new_angle_o13<=270
            deviation = diff_angles;
        elseif new_angle_o13>270 && new_angle_o13<=360
            deviation = -1 * diff_angles;
        end
        all_answers(i, 15) = deviation;
        
        % Calculating relative to the actual segment direction
        if (all_answers(i,1) == 3 && all_answers(i,2) == 0) || (all_answers(i,1) == 2 && all_answers(i,2) == 1)
            segment_direction = 0;
        elseif (all_answers(i,1) == 1 && all_answers(i,2) == 0) || (all_answers(i,1) == 2 && all_answers(i,2) == 3)
            segment_direction = 90;
        elseif (all_answers(i,1) == 0 && all_answers(i,2) == 3) || (all_answers(i,1) == 1 && all_answers(i,2) == 2)
            segment_direction = 180;
        elseif (all_answers(i,1) == 3 && all_answers(i,2) == 2) || (all_answers(i,1) == 0 && all_answers(i,2) == 1)
            segment_direction = 270;
        end
        realdirection_segment_diff = abs(new_angle_o13 - segment_direction);
        estdirection_segment_diff = abs(new_est_angle - segment_direction);
        if abs(realdirection_segment_diff) > 180, realdirection_segment_diff = realdirection_segment_diff - 360*sign(realdirection_segment_diff); end
        if abs(estdirection_segment_diff) > 180, estdirection_segment_diff = estdirection_segment_diff - 360*sign(estdirection_segment_diff); end
        all_answers(i, 16) = abs(realdirection_segment_diff) - abs(estdirection_segment_diff);
    end
    
    disp("Mean angular error within between")
    % Difference between within / between segment comparisons in angular error
    [~, p] = ttest2(all_answers(locs_within, 11), all_answers(locs_between, 11));
    disp([nanmean(all_answers(locs_within, 11)) nanmean(all_answers(locs_between, 11)) p])
    all_p_angular_error_within_between(f) = p;
    
    % Overall errors
    all_angular_errors_overall(f) = nanmean(all_answers(all_relevant_locations, 11));
    all_num_answers(f) = num_trials - sum(isnan(all_answers(all_relevant_locations, 11)));
    
    % Random error distributions
    temp_all = []; num_iters = 10000;
    for i=1:num_iters, temp = randi(360); temp = temp-all_answers(randi(48),10); temp_all(i) = abs(mod(temp + 180, 360) - 180); end
    p = sum(temp_all <= mean(all_answers(:,11))) / num_iters;
    disp(strcat("average random error = ", num2str(nanmean(temp_all))));
    disp(strcat("p overall vs. random = ", num2str(p)));
    
    all_angular_errors_within(f) = nanmean(all_answers(locs_within, 11));
    all_angular_errors_between(f) = nanmean(all_answers(locs_between, 11));
    
    all_angular_errors_within_orth(f) = nanmean(all_answers(locs_within_orth, 11));
    all_angular_errors_between_orth(f) = nanmean(all_answers(locs_between_orth, 11));
    
    
    answers_all_subjects{f} = all_answers;
end



% Calculating the response accuracy in terms of only left vs. right
% judgment (to be able to compare to accuracy in fMRI behavioral results)
for f = 1:length(answers_all_subjects)
    aaa = answers_all_subjects{f}(:,[8 10]);        % Response and real angle
    correct_responses = (sum(aaa<=180,2)==2 | sum(aaa>180,2)==2);    % either both correct and given directions are to the right or to the left
    accuracy_overall(f) = sum(correct_responses) / length(correct_responses);
end


% Calculating if there is bias toward vertical axis for between or horizontal for within (bias toward superordinate segment direction)

for f = 1:length(answers_all_subjects)
    all_corr_superordinate(f) = corr(answers_all_subjects{f}(:,9), answers_all_subjects{f}(:,15));
    all_bias_within(f) = mean(answers_all_subjects{f}(find(answers_all_subjects{f}(:,9) == 1), 15));
    all_bias_between(f) = mean(answers_all_subjects{f}(find(answers_all_subjects{f}(:,9) == -1), 15));
    all_segment_direction_dev(f) = mean(answers_all_subjects{f}(:,16));
    all_segment_direction_dev_within(f) = mean(answers_all_subjects{f}(find(answers_all_subjects{f}(:,9) == 1),16));
    all_segment_direction_dev_between(f) = mean(answers_all_subjects{f}(find(answers_all_subjects{f}(:,9) == -1),16));
    [~, p_within_subject(f)] = ttest(answers_all_subjects{f}(find(answers_all_subjects{f}(:,9) == 1), 15), answers_all_subjects{f}(find(answers_all_subjects{f}(:,9) == -1), 15));
end
