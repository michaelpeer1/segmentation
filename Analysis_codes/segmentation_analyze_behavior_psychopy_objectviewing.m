function current_subj_results = analyze_segmentation_psychopy_objectviewing(logfiles)
% Analyzing the object viewing stage behavioral results (psychopy experiment)
%
% current_subj_results = analyze_segmentation_objectviewing(logfiles)
% logfiles - cell array with the names of logfiles from the experiment

current_subj_num_correct = 0;
current_subj_num_overall = 0;
current_subj_RTs_same_segment = [];
current_subj_RTs_different_segment = [];
current_subj_RTs_same_orth_segment = [];
current_subj_RTs_different_orth_segment = [];
all_subject_RTs = [];
all_subject_objects_distance_to_previous = [];
all_subject_objects_distance_to_previous_same_seg_adj_quad = [];
all_subject_objects_distance_to_previous_same_orth_seg_adj_quad = [];
all_subject_RTs_same_seg_adj_quad = [];
all_subject_RTs_same_orth_seg_adj_quad = [];

% The distances between all objects
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
objects_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects


for f = 1:length(logfiles)
    %% Loading the data
    %     [~,~,results] = xlsread(logfiles{f});
    fid = fopen(logfiles{f});
    num_columns = 46; temp_string = ''; for i=1:num_columns, temp_string = [temp_string '%s ']; end
    data = textscan(fid, temp_string, 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:length(results(:))        % Changing strings to numbers
        if ~isnan(str2double(results{i}))
            results{i} = str2double(results{i});
        end
    end

    all_stimuli_locations = results(4:end-1, find(strcmp(results(1,:), 'Stimulus_unity_position')));
    all_isrotated = results(4:end-1, find(strcmp(results(1,:), 'Is_rotated')));
    %     all_responses = results(3:end, 21);       % PC responses
    %     all_RTs = results(3:end, 22);
    all_responses_stimulus = results(4:end-1, find(strcmp(results(1,:), 'key_resp.keys')));             % Response box responses
    all_RTs_stimulus = results(4:end-1, find(strcmp(results(1,:), 'key_resp.rt')));
    all_responses_ISI = results(4:end-1, find(strcmp(results(1,:), 'key_resp_ISI.keys')));
    all_RTs_ISI = results(4:end-1, find(strcmp(results(1,:), 'key_resp_ISI.rt')));
    num_trials = length(all_stimuli_locations);

    % Combining responses for stimulus and ISI durations
    all_responses = cell(num_trials, 1); all_RTs = cell(num_trials, 1);
    for i = 1:num_trials
        % If there was a response during the question presentation
        if ~strcmp(all_responses_stimulus{i}, 'None')
            all_responses{i} = all_responses_stimulus{i};
            all_RTs{i} = all_RTs_stimulus{i};
        % If there was a response during the following ISI
        elseif ~strcmp(all_responses_ISI{i}, 'None')
            all_responses{i} = all_responses_ISI{i};
            all_RTs{i} = all_RTs_ISI{i} + 1;
        else
            all_responses{i} = nan;
            all_RTs{i} = nan;
        end
    end
    
    % Changing the responses to numbers
    %     all_responses(strcmp('left', all_responses)) = {-1};
    %     all_responses(strcmp('down', all_responses)) = {0};
    %     all_responses(strcmp('right', all_responses)) = {1};
    %     all_responses(strcmp('None', all_responses)) = {nan};
    all_responses(strcmp('b', all_responses)) = {-1};
    all_responses(strcmp('y', all_responses)) = {0};
    all_responses(strcmp('g', all_responses)) = {1};
    all_responses(strcmp('r', all_responses)) = {1};
    
    % Changing blank stimuli to nan
    all_stimuli_locations(strcmp('None', all_stimuli_locations)) = {nan};
    for i = 1:length(all_stimuli_locations), if isempty(all_stimuli_locations{i}), all_stimuli_locations{i} = nan; end, end
    for i = 1:length(all_responses), if isempty(all_responses{i}), all_responses{i} = nan; end, end
    for i = 1:length(all_isrotated), if isempty(all_isrotated{i}), all_isrotated{i} = nan; end, end
    for i = 1:length(all_RTs), if isempty(all_RTs{i}), all_RTs{i} = nan; end, end

    % Changing all cell arrays to matrices
    all_stimuli_locations = cell2mat(all_stimuli_locations);
    all_isrotated = cell2mat(all_isrotated);
    all_responses = cell2mat(all_responses);
    all_RTs = cell2mat(all_RTs);
    
    
    %% Analyzing percent correct
    num_correct = nansum(all_isrotated == all_responses);
    %     num_overall = sum(~isnan(all_stimuli_locations));
    num_overall = sum(~isnan(all_responses));
    
    current_subj_num_correct = current_subj_num_correct + num_correct;
    current_subj_num_overall = current_subj_num_overall + num_overall;
    
    
    %% Analyzing RT priming by segment
    % Identifying the segment of each stimulus
    all_segments = nan(length(all_isrotated), 1); all_orth_segments = nan(length(all_isrotated), 1); all_quadrants = nan(length(all_isrotated), 1); 
    all_segments(all_stimuli_locations < 4 | all_stimuli_locations > 11) = 1;
    all_segments(all_stimuli_locations >= 4 & all_stimuli_locations <= 11) = 2;
    all_orth_segments(all_stimuli_locations < 8) = 1;
    all_orth_segments(all_stimuli_locations >= 8) = 2;
    all_quadrants(all_stimuli_locations <= 3) = 1; all_quadrants(all_stimuli_locations >= 4 & all_stimuli_locations <= 7) = 2;
    all_quadrants(all_stimuli_locations >= 8 & all_stimuli_locations <= 11) = 3; all_quadrants(all_stimuli_locations >= 12 & all_stimuli_locations <= 15) = 4;

    
    % Figuring out which response had the same or different segment before it
    all_same_segment = [nan; (all_segments(1:end-1) == all_segments(2:end))];
    all_same_orth_segment = [nan; (all_orth_segments(1:end-1) == all_orth_segments(2:end))];
    all_same_quadrant = [nan; (all_quadrants(1:end-1) == all_quadrants(2:end))];
    % Removing irrelevant locations (during / after blank, same stimulus appearing twice in a row, etc.)
    all_relevant_locations = [0; all_stimuli_locations(1:end-1) ~= all_stimuli_locations(2:end) & ~isnan(all_stimuli_locations(1:end-1) + all_stimuli_locations(2:end))];
    all_same_segment(all_relevant_locations == 0) = nan;
    all_same_orth_segment(all_relevant_locations == 0) = nan;
    all_same_quadrant(all_relevant_locations == 0) = nan;
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
    objects_distance_to_previous = nan(size(all_stimuli_locations));
    for i = 2:length(all_stimuli_locations)
        if ~isnan(all_stimuli_locations(i) + all_stimuli_locations(i-1))
            objects_distance_to_previous(i) = objects_distances(all_stimuli_locations(i)+1, all_stimuli_locations(i-1)+1);
        end
    end
    all_subject_RTs = [all_subject_RTs; all_RTs];
    all_subject_objects_distance_to_previous = [all_subject_objects_distance_to_previous; objects_distance_to_previous];
    all_subject_objects_distance_to_previous_same_seg_adj_quad = [all_subject_objects_distance_to_previous_same_seg_adj_quad; objects_distance_to_previous(all_same_segment == 1 & all_same_quadrant == 0)];
    all_subject_objects_distance_to_previous_same_orth_seg_adj_quad = [all_subject_objects_distance_to_previous_same_orth_seg_adj_quad; objects_distance_to_previous(all_same_orth_segment == 1 & all_same_quadrant == 0)];
    all_subject_RTs_same_seg_adj_quad = [all_subject_RTs_same_seg_adj_quad; all_RTs(all_same_segment == 1 & all_same_quadrant == 0)];
    all_subject_RTs_same_orth_seg_adj_quad = [all_subject_RTs_same_orth_seg_adj_quad; all_RTs(all_same_orth_segment == 1 & all_same_quadrant == 0)];
end

current_subj_results  = struct();
current_subj_results.num_correct = current_subj_num_correct;
current_subj_results.num_overall = current_subj_num_overall;
current_subj_results.RTs_same_segment = nanmean(current_subj_RTs_same_segment);
current_subj_results.RTs_different_segment = nanmean(current_subj_RTs_different_segment);
current_subj_results.RTs_same_orth_segment = nanmean(current_subj_RTs_same_orth_segment);
current_subj_results.RTs_different_orth_segment = nanmean(current_subj_RTs_different_orth_segment);


% Calculating the correlation between RT and distance to previous item
current_subj_results.corr_distances_RTs = corr(all_subject_objects_distance_to_previous, all_subject_RTs, 'rows', 'complete');
current_subj_results.corr_distances_RTs_within_seg_adj_quad = corr(all_subject_objects_distance_to_previous_same_seg_adj_quad, all_subject_RTs_same_seg_adj_quad, 'rows', 'complete');
current_subj_results.corr_distances_RTs_within_orth_seg_adj_quad = corr(all_subject_objects_distance_to_previous_same_orth_seg_adj_quad, all_subject_RTs_same_orth_seg_adj_quad, 'rows', 'complete');



% logfiles1 = {'C:\Users\michaelpeer1\Desktop\Segmentation_fMRI\Subjects_fMRI_pilot\Pilot_subj01_seed_1\Tasks\data\s1_Segmentation_2_objectviewing_run1_2019_Apr_18_2114.csv',...
%     'C:\Users\michaelpeer1\Desktop\Segmentation_fMRI\Subjects_fMRI_pilot\Pilot_subj01_seed_1\Tasks\data\s1_Segmentation_3_objectviewing_run2_2019_Apr_18_2126.csv'};
% logfiles2 = {'C:\Users\michaelpeer1\Desktop\Segmentation_fMRI\Subjects_fMRI_pilot\Pilot_subj03_seed_3\Tasks\data_corrected\s3_Segmentation_2_objectviewing_run1_2019_Apr_18_2323_corrected.csv',...
%     'C:\Users\michaelpeer1\Desktop\Segmentation_fMRI\Subjects_fMRI_pilot\Pilot_subj03_seed_3\Tasks\data_corrected\s3_Segmentation_3_objectviewing_run2_2019_Apr_18_2334_corrected.csv'};
% 
% subj1_results = analyze_segmentation_psychopy_objectviewing(logfiles1);
% subj2_results = analyze_segmentation_psychopy_objectviewing(logfiles2);
% 
% disp([nanmean(subj1_results.RTs_same_segment), nanmean(subj1_results.RTs_different_segment)])
% disp([nanmean(subj2_results.RTs_same_segment), nanmean(subj2_results.RTs_different_segment)])
% 
% [~, p] = ttest2(subj1_results.RTs_same_segment, subj1_results.RTs_different_segment); disp(p)
% [~, p] = ttest2(subj2_results.RTs_same_segment, subj2_results.RTs_different_segment); disp(p)

