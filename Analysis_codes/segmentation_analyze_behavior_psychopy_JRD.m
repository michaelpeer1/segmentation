function current_subj_results = analyze_segmentation_psychopy_JRD(logfiles)
% current_subj_results = analyze_segmentation_psychopy_JRD(logfiles)


current_subj_RTs_same_segment = [];
current_subj_RTs_different_segment = [];
current_subj_RTs_same_orth_segment = [];
current_subj_RTs_different_orth_segment = [];
current_subj_num_overall = 0;
current_subj_num_correct = 0;
current_subj_num_answered = 0;
current_subj_num_correct_within_segment = 0;
current_subj_num_correct_between_segments = 0;
current_subj_match_objects_start_quadrant = zeros(16);
current_subj_match_objects_all = zeros(16);
current_subj_angles_accuracy = zeros(360,1);
current_subj_angles_num = zeros(360,1);
current_subj_responses_per_object = zeros(16,2);
all_subject_RTs = [];
all_subject_objects_distance_to_previous = [];
all_subject_objects_distance_to_previous_same_seg_adj_quad = [];
all_subject_objects_distance_to_previous_same_orth_seg_adj_quad = [];
all_subject_RTs_same_seg_adj_quad = [];
all_subject_RTs_same_orth_seg_adj_quad = [];

for f = 1:length(logfiles)
    % Loading the data
    %     [~,~,results] = xlsread(logfiles{f});
    fid = fopen(logfiles{f});
    num_columns = 60; temp_string = ''; for i=1:num_columns, temp_string = [temp_string '%s ']; end
    data = textscan(fid, temp_string, 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:length(results(:))        % Changing strings to numbers
        if ~isnan(str2double(results{i}))
            results{i} = str2double(results{i});
        end
    end
    for i = 1:length(results(:)), if isempty(results{i}), results{i} = nan; end, end
    
    all_object1 = cell2mat(results(4:end - 1, find(strcmp(results(1,:), 'Object1_unity_location'))));
    all_object2 = cell2mat(results(4:end - 1, find(strcmp(results(1,:), 'Object2_unity_location'))));
    all_object3 = cell2mat(results(4:end - 1, find(strcmp(results(1,:), 'Object3_unity_location'))));
    %     all_responses_stimulus = results(3:end - 1, 25);      % PC task
    %     all_responses_ISI = results(3:end - 1, 40);
    %     all_RTs_stimulus = results(3:end - 1, 44);
    %     all_RTs_ISI = results(3:end - 1, 41);
    all_responses_stimulus = results(4:end - 1, find(strcmp(results(1,:), 'key_resp.keys')));            % fMRI task
    all_responses_ISI = results(4:end - 1, find(strcmp(results(1,:), 'key_resp_ISI.keys')));
    all_RTs_stimulus = results(4:end - 1, find(strcmp(results(1,:), 'key_resp.rt')));
    all_RTs_ISI = results(4:end - 1, find(strcmp(results(1,:), 'key_resp_ISI.rt')));
    num_trials = length(all_object1);
    
    % Combining responses for stimulus and ISI durations
    all_responses = nan(num_trials, 1); all_RTs = nan(num_trials, 1);
    for i = 1:num_trials
        % If there was a response during the question presentation
        %         if strcmp(all_responses_stimulus{i}, 'left')
        if strcmp(all_responses_stimulus{i}, 'b')
            all_responses(i) = -1;
            all_RTs(i) = all_RTs_stimulus{i};
            %         elseif strcmp(all_responses_stimulus{i}, 'right')
        elseif strcmp(all_responses_stimulus{i}, 'y') || strcmp(all_responses_stimulus{i}, 'g') || strcmp(all_responses_ISI{i}, 'r')
            all_responses(i) = 1;
            all_RTs(i) = all_RTs_stimulus{i};
        end
        % If there was a response during the following ISI
        %         if strcmp(all_responses_ISI{i}, 'left')
        if strcmp(all_responses_ISI{i}, 'b')
            all_responses(i) = -1;
            all_RTs(i) = all_RTs_ISI{i} + 5;
            %         elseif strcmp(all_responses_ISI{i}, 'right')
        elseif strcmp(all_responses_ISI{i}, 'y') || strcmp(all_responses_ISI{i}, 'g') || strcmp(all_responses_ISI{i}, 'r')
            all_responses(i) = 1;
            all_RTs(i) = all_RTs_ISI{i} + 5;
        end
    end
    
    
    %% Calculate correct answers
    locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
        42 -16; 58 -50; 27 -57; 17 -25; ...
        -16 -42; -50 -58; -57 -27; -25 -17; ...
        -42 16; -58 50; -27 57; -17 25];
    objects_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects

    current_subj_num_overall = current_subj_num_overall + num_trials;
    
    for i = 1:num_trials
        if ~isnan(all_object1(i))
            % Figuring out the real angle of stimulus 3 compared to the viewing direction
            o1 = all_object1(i) + 1;
            o2 = all_object2(i) + 1;
            o3 = all_object3(i) + 1;
            
            v1 = locations_all_objects(o2, :) - locations_all_objects(o1, :);   % direction between object 1 and 2
            v2 = locations_all_objects(o3, :) - locations_all_objects(o1, :);   % direction between object 1 and 3
            
            angle = atan2d(v1(2),v1(1)) - atan2d(v2(2),v2(1));
            if abs(angle) > 180, angle = angle - 360*sign(angle); end
            if angle < 0, angle = 360 + angle; end
            
            % Figuring out if there was a response
            if ~isnan(all_responses(i))
                current_subj_num_answered = current_subj_num_answered + 1;
            end
            
            % Figuring out if the response was correct
            is_correct = nan;
            if angle < 180 && all_responses(i) == 1     % Correct answer - stimulus to the right
                current_subj_num_correct = current_subj_num_correct + 1;
                is_correct = 1;
            elseif angle > 180 && all_responses(i) == -1     % Correct answer - stimulus to the left
                current_subj_num_correct = current_subj_num_correct + 1;
                is_correct = 1;
            elseif ~isnan(all_responses(i))     % Wrong answer
                is_correct = 0;
            end
            
            % Figuring out if the comparison is within or between segments
            if (o1 < 4 || o1 > 11), s1 = 1; elseif (o1 >= 4 || o1 <= 11), s1 = 2; end
            if (o3 < 4 || o3 > 11), s2 = 1; elseif (o3 >= 4 || o3 <= 11), s2 = 2; end
            if is_correct == 1 && s1 == s2
                current_subj_num_correct_within_segment = current_subj_num_correct_within_segment + 1;
            elseif is_correct == 1 && s1 ~= s2
                current_subj_num_correct_between_segments = current_subj_num_correct_between_segments + 1;
            end
            
            % Recording correct answers by angle
            if ~isnan(is_correct)
                current_subj_angles_accuracy(round(angle)) = current_subj_angles_accuracy(round(angle)) + is_correct;
            end
            current_subj_angles_num(round(angle)) = current_subj_angles_num(round(angle)) + 1;
            
            % Recording responses per object
            if all_responses(i) == 1
                current_subj_responses_per_object(o1,1) = current_subj_responses_per_object(o1,1) + 1;
            elseif all_responses(i) == -1
                current_subj_responses_per_object(o1,2) = current_subj_responses_per_object(o1,2) + 1;
            end
        end
    end
    
    
    %% Calculate RT priming - object1 location (imagined standing location)
    all_segments = nan(num_trials, 1); all_orth_segments = nan(num_trials, 1); all_quadrants = nan(num_trials, 1);
    all_segments(all_object1 < 4 | all_object1 > 11) = 1;
    all_segments(all_object1 >= 4 & all_object1 <= 11) = 2;
    all_orth_segments(all_object1 < 8) = 1;
    all_orth_segments(all_object1 >= 8) = 2;
    all_quadrants(all_object1 <= 3) = 1; all_quadrants(all_object1 >= 4 & all_object1 <= 7) = 2;
    all_quadrants(all_object1 >= 8 & all_object1 <= 11) = 3; all_quadrants(all_object1 >= 12 & all_object1 <= 15) = 4;

    % Figuring out which response had the same or different segment before it
    all_same_segment = [nan; (all_segments(1:end-1) == all_segments(2:end))];
    all_same_orth_segment = [nan; (all_orth_segments(1:end-1) == all_orth_segments(2:end))];
    all_same_quadrant = [nan; (all_quadrants(1:end-1) == all_quadrants(2:end))];

    % Getting RTs
    all_RTs_same_segment = all_RTs(all_same_segment == 1);
    all_RTs_different_segment = all_RTs(all_same_segment == 0);
    all_RTs_same_orth_segment = all_RTs(all_same_orth_segment == 1);
    all_RTs_different_orth_segment = all_RTs(all_same_orth_segment == 0);
    
    current_subj_RTs_same_segment = [current_subj_RTs_same_segment; all_RTs_same_segment];
    current_subj_RTs_different_segment = [current_subj_RTs_different_segment; all_RTs_different_segment];
    current_subj_RTs_same_orth_segment = [current_subj_RTs_same_orth_segment; all_RTs_same_orth_segment];
    current_subj_RTs_different_orth_segment = [current_subj_RTs_different_orth_segment; all_RTs_different_orth_segment];
    
    % Distance between object and previous object
    objects_distance_to_previous = nan(size(all_object1));
    for i = 2:length(all_object1)
        if ~isnan(all_object1(i) + all_object1(i-1))
            objects_distance_to_previous(i) = objects_distances(all_object1(i)+1, all_object1(i-1)+1);
        end
    end
    all_subject_RTs = [all_subject_RTs; all_RTs];
    all_subject_objects_distance_to_previous = [all_subject_objects_distance_to_previous; objects_distance_to_previous];
    all_subject_objects_distance_to_previous_same_seg_adj_quad = [all_subject_objects_distance_to_previous_same_seg_adj_quad; objects_distance_to_previous(all_same_segment == 1 & all_same_quadrant == 0)];
    all_subject_objects_distance_to_previous_same_orth_seg_adj_quad = [all_subject_objects_distance_to_previous_same_orth_seg_adj_quad; objects_distance_to_previous(all_same_orth_segment == 1 & all_same_quadrant == 0)];
    all_subject_RTs_same_seg_adj_quad = [all_subject_RTs_same_seg_adj_quad; all_RTs(all_same_segment == 1 & all_same_quadrant == 0)];
    all_subject_RTs_same_orth_seg_adj_quad = [all_subject_RTs_same_orth_seg_adj_quad; all_RTs(all_same_orth_segment == 1 & all_same_quadrant == 0)];
    
    % Checking match between objects in task
    for i=1:num_trials
        current_subj_match_objects_start_quadrant(all_object1(i)+1,all_object2(i)+1) = current_subj_match_objects_start_quadrant(all_object1(i)+1,all_object2(i)+1)+1;
        current_subj_match_objects_start_quadrant(all_object2(i)+1,all_object1(i)+1) = current_subj_match_objects_start_quadrant(all_object2(i)+1,all_object1(i)+1)+1;
        
        current_subj_match_objects_all(all_object1(i)+1,all_object2(i)+1) = current_subj_match_objects_all(all_object1(i)+1,all_object2(i)+1)+1;
        current_subj_match_objects_all(all_object2(i)+1,all_object1(i)+1) = current_subj_match_objects_all(all_object2(i)+1,all_object1(i)+1)+1;
        current_subj_match_objects_all(all_object1(i)+1,all_object3(i)+1) = current_subj_match_objects_all(all_object1(i)+1,all_object3(i)+1)+1;
        current_subj_match_objects_all(all_object3(i)+1,all_object1(i)+1) = current_subj_match_objects_all(all_object3(i)+1,all_object1(i)+1)+1;
        current_subj_match_objects_all(all_object2(i)+1,all_object3(i)+1) = current_subj_match_objects_all(all_object2(i)+1,all_object3(i)+1)+1;
        current_subj_match_objects_all(all_object3(i)+1,all_object2(i)+1) = current_subj_match_objects_all(all_object3(i)+1,all_object2(i)+1)+1;
    end
    
end

current_subj_results.num_overall = current_subj_num_overall;
current_subj_results.num_answered = current_subj_num_answered;
current_subj_results.num_correct = current_subj_num_correct;
current_subj_results.num_correct_within_segment = current_subj_num_correct_within_segment;
current_subj_results.num_correct_between_segments = current_subj_num_correct_between_segments;
current_subj_results.RTs_same_segment = nanmean(current_subj_RTs_same_segment);
current_subj_results.RTs_different_segment = nanmean(current_subj_RTs_different_segment);
current_subj_results.RTs_same_orth_segment = nanmean(current_subj_RTs_same_orth_segment);
current_subj_results.RTs_different_orth_segment = nanmean(current_subj_RTs_different_orth_segment);
current_subj_results.match_objects_start_quadrant = current_subj_match_objects_start_quadrant;
current_subj_results.match_objects_all = current_subj_match_objects_all;
current_subj_results.angles_accuracy = current_subj_angles_accuracy;
current_subj_results.angles_num = current_subj_angles_num;
current_subj_results.responses_per_object = current_subj_responses_per_object;

% Calculating the correlation between RT and distance to previous item
current_subj_results.corr_distances_RTs = corr(all_subject_objects_distance_to_previous, all_subject_RTs, 'rows', 'complete');
current_subj_results.corr_distances_RTs_within_seg_adj_quad = corr(all_subject_objects_distance_to_previous_same_seg_adj_quad, all_subject_RTs_same_seg_adj_quad, 'rows', 'complete');
current_subj_results.corr_distances_RTs_within_orth_seg_adj_quad = corr(all_subject_objects_distance_to_previous_same_orth_seg_adj_quad, all_subject_RTs_same_orth_seg_adj_quad, 'rows', 'complete');
