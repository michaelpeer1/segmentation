function segmentation_create_SPM_design_JRD(psychopy_log_filename, output_mat_filename)
% segmentation_create_SPM_design_JRD(psychopy_log_filename, output_mat_filename)
%
% This script analyzes the psychopy logfile for the JRD
% paradigm in the segmentation project, and extracts the relevant stimulus 
% appearance times for SPM analysis.
% Each object appears for 5s with varying ISI.
% - The time for onsets and durations is saved in seconds.
% - Object names are by their position in unity in this subject, not their
% original object index.

num_objects = 16;
duration_object_appearance = 5;
duration_start_rest = 10;

% Initializing arrays
names = cell(1, num_objects);
onsets = cell(1, num_objects);
durations = cell(1, num_objects);

for i = 1:num_objects
    names{i} = ['object' num2str(i - 1)];
end

% Reading the log file
temp_file = fopen(psychopy_log_filename, 'r');
data = textscan(temp_file, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'Delimiter',',','HeaderLines',3);     % Read the first 12 columns of the log file 
fclose(temp_file);

counter_current_time = duration_start_rest;
for i = 1:length(data{1})     % Going over each line of the log file
    if ~isempty(data{1}{i})
        current_stimulus_num = str2double(data{9}{i}) + 1;     % The unity index of the stimulus
        onsets{current_stimulus_num} = [onsets{current_stimulus_num}, counter_current_time];
        counter_current_time = counter_current_time + duration_object_appearance + str2double(data{12}{i});
        durations{current_stimulus_num} = [durations{current_stimulus_num}, duration_object_appearance];
    end
end

save(output_mat_filename, 'names', 'onsets', 'durations');




% Saving all of the information
all_info = struct();
all_info.object1_original = str2double(data{4});
all_info.object2_original = str2double(data{5});
all_info.object3_original = str2double(data{6});
all_info.quadrant_start = str2double(data{7});
all_info.quadrant_target = str2double(data{8});
all_info.object1_unity_location = str2double(data{9});
all_info.object2_unity_location = str2double(data{10});
all_info.object3_unity_location = str2double(data{11});
% Getting the response
responses = {data{29}, data{44}}; all_responses={}; 
all_responses(~strcmp(responses{1},'None')) = responses{1}(~strcmp(responses{1},'None')); 
all_responses(~strcmp(responses{2},'None')) = responses{2}(~strcmp(responses{2},'None')); 
all_responses_new = nan(1,length(all_responses)); all_responses_new(strcmp(all_responses, 'b')) = -1; all_responses_new(strcmp(all_responses, 'y') | strcmp(all_responses, 'g')  | strcmp(all_responses, 'r')) = 1;
all_info.response = all_responses_new';
% Getting the angles and whether the response is correct
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
for i = 1:length(all_info.response)
    if ~isnan(all_info.object1_original(i))
        o1 = all_info.object1_unity_location(i) + 1; o2 = all_info.object2_unity_location(i) + 1; o3 = all_info.object3_unity_location(i) + 1;
        % Figuring out the real angle of stimulus 3 compared to the viewing direction
        v1 = locations_all_objects(o2, :) - locations_all_objects(o1, :);   % direction between object 1 and 2
        v2 = locations_all_objects(o3, :) - locations_all_objects(o1, :);   % direction between object 1 and 3
        all_info.viewing_direction(i) = atan2d(v1(2),v1(1)); all_info.target_direction(i) = atan2d(v2(2),v2(1));
        % Figuring out if the response was correct
        all_info.is_correct(i) = nan;
        angle = atan2d(v1(2),v1(1)) - atan2d(v2(2),v2(1));  % Angle between viewing and target directions
        if abs(angle) > 180, angle = angle - 360*sign(angle); end
        if angle < 0, angle = 360 + angle; end
        if angle < 180 && all_info.response(i) == 1, all_info.is_correct(i) = 1;     % Correct answer - stimulus to the right
        elseif angle > 180 && all_info.response(i) == -1, all_info.is_correct(i) = 1;     % Correct answer - stimulus to the left
        elseif ~isnan(all_info.response(i)), all_info.is_correct(i) = 0; % Wrong answer
        end
    end
end
% Saving
[a, b] = fileparts(output_mat_filename);
save(fullfile(a, ['all_info_JRD_run_', b(end-3:end), '.mat']), 'all_info');
