function segmentation_create_SPM_design_object_viewing(psychopy_log_filename, output_mat_filename)
% segmentation_create_SPM_design_object_viewing(psychopy_log_filename, output_mat_filename)
%
% This script analyzes the psychopy logfile for the object viewing
% paradigm in the segmentation project, , and extracts the relevant 
% stimulus appearance times for SPM analysis.
% Each object appears for 1s with 1s rest.
% - The time for onsets and durations is saved in seconds.
% - Object names are by their position in unity in this subject, not their
% original object index.

num_objects = 16;
duration_object_appearance = 1;
duration_object_rest = 1;
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
data = textscan(temp_file, '%s %s %s %s %s %f %*[^\n]', 'Delimiter',',','HeaderLines',3);     % Read the first four columns of the log file (stimulus filename, stimulus index, stimulus name, stimulus unity position)
fclose(temp_file);

for i = 1:length(data{1})     % Going over each line of the log file
    if ~isempty(data{2}{i})
        if str2double(data{2}{i}) ~= 17      % 17 is the sign for blank
            current_stimulus_num = str2double(data{4}{i}) + 1;     % The unity index of the stimulus
            onsets{current_stimulus_num} = [onsets{current_stimulus_num}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
            durations{current_stimulus_num} = [durations{current_stimulus_num}, (duration_object_appearance + duration_object_rest)];
        end
    end
end

save(output_mat_filename, 'names', 'onsets', 'durations');
