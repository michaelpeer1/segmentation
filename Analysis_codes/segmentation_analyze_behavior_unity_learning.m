% Analyzing the learning stage

% Pilot data analysis - old
% filenames = {};
% filenames{1} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1106.csv';     % MP
% filenames{2} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190225_1509.csv';     % RM
% filenames{3} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190226_1608.csv';     % ZL_first_attempt
% filenames{4} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1630.csv';     % ZL1
% filenames{5} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1642.csv';     % ZL2
% filenames{6} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1656.csv';     % ZL3
% filenames{7} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190301_1703.csv';     % ZL4
% filenames{8} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1144.csv';     % AT1
% filenames{9} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1154.csv';     % AT2
% filenames{10} = 'C:\Unity projects\Segmentation_paradigm\Assets\recorded_events_20190304_1209.csv';     % AT3

% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v1_river\Subjects_file.xlsx');
% num_subjects = max(find(~cellfun(@isnan, subjects_file(4:end,1))));
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation of space\Behavioral_data_v2_wall\Subjects_file_v2_wall.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_behavioral_pilot\Subjects_file_fMRI_behav_pilot.xlsx');
[~, ~, subjects_file] = xlsread('/Users/mpeer/Desktop/Segmentation_data/fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx');
[~, ~, subjects_file] = xlsread('C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\Segmentation project\Segmentation_data\fMRI_experiment_logfiles\Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\fMRI_experiment\Subjects_file_fMRI.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\5_behavioral_wall_new\Subjects_file_v5_wall_new.xlsx');
% [~, ~, subjects_file] = xlsread('E:\Epstein_lab\Projects\1 - Segmentation\Subjects_data\Other versions\6_separate_learning\Subjects_file_v6_separate_learning.xlsx');
num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
filenames = {};
for i=1:num_subjects
%     if ~isnan(subjects_file{3+i, 9})       % Behavioral version of the task
%         filenames{end+1} = [fullfile(subjects_file{3+i, 7}, subjects_file{3+i, 9}) '.csv'];
    if ~isnan(subjects_file{3+i, 11})       % fMRI version of the task
        filenames{end+1} = [fullfile(subjects_file{3+i, 8}, subjects_file{3+i, 10}) '.csv'];        % PC version
        filenames{end+1} = [fullfile(subjects_file{3+i, 8}, subjects_file{3+i, 11}) '.csv'];
%         filenames{end+1} = [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 10}) '.csv'];        % Mac version
%         filenames{end+1} = [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 11}) '.csv'];
    end
end


answers_all_subjects = cell(1, length(filenames));
all_mistakes_num = zeros(length(filenames), 1);
all_mistakes_segment_schema_preserving = zeros(length(filenames), 1);
all_mistakes_segment_overlay_preserving = zeros(length(filenames), 1);
all_mistakes_quadrant_schema_preserving = zeros(length(filenames), 1);
all_mistakes_same_segment = zeros(length(filenames), 1);
all_mistakes_same_orth_segment = zeros(length(filenames), 1);
all_mistakes_same_quadrant = zeros(length(filenames), 1);
all_mistakes_percent_segment_schema_preserving = zeros(length(filenames), 1);
all_mistakes_percent_segment_overlay_preserving = zeros(length(filenames), 1);
all_mistakes_percent_quadrant_schema_preserving = zeros(length(filenames), 1);
all_mistakes_percent_same_segment = zeros(length(filenames), 1);
all_mistakes_percent_same_quadrant = zeros(length(filenames), 1);
all_mistakes_same_segment_different_quadrant = zeros(length(filenames), 1);
all_mistakes_same_orth_segment_different_quadrant = zeros(length(filenames), 1);
all_mistakes_diagonal_quadrant = zeros(length(filenames), 1);
num_levels = 5;
all_mistakes_per_level = nan(length(filenames), num_levels);
all_times_per_level = nan(length(filenames), num_levels);
all_objects_per_level = nan(length(filenames), num_levels);

for f = 1:length(filenames)
    %% Reading the data
    %     [~,~,results] = xlsread(filenames{f});
    fid = fopen(filenames{f});
    data = textscan(fid,'%s %s %s %s %s %s %s %s %s', 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:size(results, 1), results{i,5} = str2double(results{i,5}); results{i,9} = str2double(results{i,9}); end   % Changing strings to numbers

    locations_object_encountered = [];
    for i = 1:size(results, 1)
        if strcmp(results{i,1}, 'Correct object encountered') || strcmp(results{i,1}, 'Wrong object encountered')
            locations_object_encountered(end+1) = i;
        elseif strcmp(results{i,1}, 'Object randomization seed')
            temp_loc_randomization = i;
        end
    end
    num_object_encounters = length(locations_object_encountered);
    all_num_object_encounters(f) = num_object_encounters;
    
    % Calculating the max level reached - for fMRI task
    temp_levels = cell2mat(results(25:end-2, 5)); num_levels = max(temp_levels);
    
    % Getting the objects order (could be different from default by randomization)
    temp_loc_randomization = find(strcmp(results(:,1), 'Object randomization seed'));
    objects_order = results((temp_loc_randomization+1):(temp_loc_randomization+16), 7);
    
    all_answers = cell(num_object_encounters, 6);   % Object, object to search, correct/wrong, object_quadrant, object_to_search_quadrant, time to find, training_level
    
    % Recording the data
    for i = 1:num_object_encounters
        all_answers{i,1} = results{locations_object_encountered(i), 7};     % Encountered object
%         all_answers{i,1} = str2num(all_answers{i,1}(7:8)); % Converting to number
        all_answers{i,1} = find(cellfun(@(x) strcmp(all_answers{i,1}, x), objects_order)) - 1;  % finding where the object is according to objects_order
        
        % Finding the current object to search
        temp_index = locations_object_encountered(i);
        while ~strcmp(results{temp_index, 1}, 'Current object to search:')
            temp_index = temp_index - 1;
            if temp_index == 0
                break;
            end
        end
        if temp_index~=0
            all_answers{i,2} = results{temp_index, 7};                          % Object to search
%             all_answers{i,2} = str2num(all_answers{i,2}(7:8)); % Converting to number
            all_answers{i,2} = find(cellfun(@(x) strcmp(all_answers{i,2}, x), objects_order)) - 1;  % finding where the object is according to objects_order
            
            all_answers{i,3} = (all_answers{i,1} == all_answers{i,2});          % Correct or wrong answer
            
            all_answers{i,4} = floor(all_answers{i,1} / 4);                     % Object quadrant
            all_answers{i,5} = floor(all_answers{i,2} / 4);                     % Object to search quadrant
            
            all_answers{i,6} = (results{locations_object_encountered(i), 9}) - (results{temp_index, 9});    % time to find
            
            all_answers{i,7} = results{temp_index, 5};
        else
            for q = 1:size(all_answers, 2), all_answers{i,q} = nan; end
        end
    end
    
    
    %% Analyzing types of mistakes - excluding level 1 (wrong object bumped into by mistake, since everything is visible)
    mistakes_segment_schema_preserving = 0;
    mistakes_segment_overlay_preserving = 0;
    mistakes_quadrant_schema_preserving = 0;
    mistakes_same_quadrant = 0;
    mistakes_same_segment = 0;
    mistakes_same_orth_segment = 0;
    mistakes_same_segment_different_quadrant = 0;
    mistakes_same_orth_segment_different_quadrant = 0;
    mistakes_diagonal_quadrant = 0;
    
    mistakes_per_level = zeros(1, num_levels);
    
    locations_wrong_answers = find([all_answers{:,3}] == 0);
    num_mistakes = length(locations_wrong_answers);
    for i = 1:num_mistakes
        o1 = all_answers{locations_wrong_answers(i),1}; o2 = all_answers{locations_wrong_answers(i),2};   % Object
        q1 = all_answers{locations_wrong_answers(i),4}; q2 = all_answers{locations_wrong_answers(i),5};   % Quadrant
        s1 = (q1 == 1 || q1 == 2); s2 = (q2 == 1 || q2 == 2);       % Segment
        orth_s1 = (q1 == 0 || q1 == 1); orth_s2 = (q2 == 0 || q2 == 1);     % Orthogonal direction to segmentation
        current_level = all_answers{locations_wrong_answers(i),7};
        
        % Schema preserving mistakes, segment level - not same segment but same location
        if (s1 ~= s2 && mod(o1,8) == mod(o2,8) && current_level ~= 1)
            mistakes_segment_schema_preserving = mistakes_segment_schema_preserving + 1;
        end
        
        % Overlay preserving mistakes, segment level - not same segment but same location (assuming overlay and not rotation)
        overlay_transitions = [1 8; 2 5; 3 6; 4 7; 5 2; 6 3; 7 4; 8 1; 9 16; 10 13; 11 14; 12 15; 13 10; 14 11; 15 12; 16 9];
        if (s1 ~= s2 && overlay_transitions(o1+1, 2) == o2+1 && current_level ~= 1)
            mistakes_segment_overlay_preserving = mistakes_segment_overlay_preserving + 1;
        end
        
        % Schema preserving mistakes, quadrant level - not same quadrant but same location
        if (mod(o1,4) == mod(o2,4)  && current_level ~= 1)
            mistakes_quadrant_schema_preserving = mistakes_quadrant_schema_preserving + 1;
        end
        
        % Mistakes by quadrant
        if (q1 == q2 && current_level ~= 1)         % Same quadrant mistakes
            mistakes_same_quadrant = mistakes_same_quadrant + 1;
        elseif (s1 == s2 && current_level ~= 1)     % Same segment, different quadrant
            mistakes_same_segment_different_quadrant = mistakes_same_segment_different_quadrant + 1;
        elseif (orth_s1 == orth_s2 && current_level ~= 1)   % Same orthogonal segment, different quadrant
            mistakes_same_orth_segment_different_quadrant = mistakes_same_orth_segment_different_quadrant + 1;
        else        % Across-diagonal quadrant
            mistakes_diagonal_quadrant = mistakes_diagonal_quadrant + 1;
        end
        
        % Same segment mistakes
        if (s1 == s2  && current_level ~= 1)
            mistakes_same_segment = mistakes_same_segment + 1;
        end
        
        % Same orthogonal segment mistakes
        if (orth_s1 == orth_s2  && current_level ~= 1)
            mistakes_same_orth_segment = mistakes_same_orth_segment + 1;
        end
        
        mistakes_per_level(current_level) = mistakes_per_level(current_level) + 1;
    end
    
    time_per_level = zeros(1, num_levels);
    for i = 1:num_levels
        start_level = min(find(cellfun(@(x) x==i, results(4:end,5)))) + 3;
        end_level = max(find(cellfun(@(x) x==i, results(4:end,5)))) + 3 - 2;
        time_per_level(i) = (results{end_level, 9} - results{start_level, 9}) / 60;
    end
    
    objects_per_level = zeros(1, num_levels);
    for i = 1:num_levels
        objects_per_level(i) = length(find(cell2mat(all_answers(:,7))==i));
    end    
    
    disp(strcat("Num mistakes: ", num2str(num_mistakes)));
    disp(strcat("Segment schema preserving mistakes percentage: ", num2str(mistakes_segment_schema_preserving / num_mistakes)));
    disp(strcat("Quadrant schema preserving mistakes percentage: ", num2str(mistakes_quadrant_schema_preserving / num_mistakes)));
    disp(strcat("Same segment mistakes percentage: ", num2str(mistakes_same_segment / num_mistakes)));
    disp(strcat("Same quadrant mistakes percentage: ", num2str(mistakes_same_quadrant / num_mistakes)));
    
    all_mistakes_num(f) = num_mistakes;
    all_mistakes_segment_schema_preserving(f) = mistakes_segment_schema_preserving;
    all_mistakes_segment_overlay_preserving(f) = mistakes_segment_overlay_preserving;
    all_mistakes_quadrant_schema_preserving(f) = mistakes_quadrant_schema_preserving;
    all_mistakes_same_segment(f) = mistakes_same_segment;
    all_mistakes_same_orth_segment(f) = mistakes_same_orth_segment;
    all_mistakes_same_quadrant(f) = mistakes_same_quadrant;
    all_mistakes_percent_segment_schema_preserving(f) = mistakes_segment_schema_preserving / num_mistakes;
    all_mistakes_percent_segment_overlay_preserving(f) = mistakes_segment_overlay_preserving / num_mistakes;
    all_mistakes_percent_quadrant_schema_preserving(f) = mistakes_quadrant_schema_preserving / num_mistakes;
    all_mistakes_percent_same_segment(f) = mistakes_same_segment / num_mistakes;
    all_mistakes_percent_same_quadrant(f) = mistakes_same_quadrant / num_mistakes;
    all_mistakes_same_segment_different_quadrant(f) = mistakes_same_segment_different_quadrant;
    all_mistakes_same_orth_segment_different_quadrant(f) = mistakes_same_orth_segment_different_quadrant;
    all_mistakes_diagonal_quadrant(f) = mistakes_diagonal_quadrant;
    all_mistakes_per_level(f, 1:length(mistakes_per_level)) = mistakes_per_level;
    all_times_per_level(f, 1:length(time_per_level)) = time_per_level;
    all_objects_per_level(f, 1:length(objects_per_level)) = objects_per_level;
    
    %%
    answers_all_subjects{f} = all_answers;
end

% all_times_per_level = cell2mat(all_times_per_level');
% all_mistakes_per_level = cell2mat(all_mistakes_per_level');
