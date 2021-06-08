function segmentation_create_SPM_design_functional_localizer(psychopy_log_filename, output_mat_filename)
% segmentation_create_SPM_design_functional_localizer(psychopy_log_filename, output_mat_filename)
%
% This script analyzes the psychopy logfile for the functional localizer
% runs in the segmentation paradigm, and extracts the relevant stimulus 
% appearance times for SPM analysis.
% Each category appears for 16s.
% - The time for onsets and durations is saved in seconds.
% - Object names are by their position in unity in this subject, not their
% original object index.

num_categories = 4;
duration_object_appearance = 16;

% Initializing arrays
names = {'Scene', 'Face', 'Object', 'Scrambled object'};
onsets = cell(1, num_categories);
durations = cell(1, num_categories);


% Reading the log file
temp_file = fopen(psychopy_log_filename, 'r');
data = textscan(temp_file, '%s %s %*[^\n]', 'Delimiter',',','HeaderLines',1);     % Read the first two columns of the log file (stimulus filename, category name)
fclose(temp_file);

counter_current_time = 0;
for i = 2:length(data{1})     % Going over each line of the log file
    if ~isempty(data{2}{i})
        if ~strcmp(data{2}{i}, data{2}{i - 1})      % If the current stimulus is different than the last one
            if strcmp(data{2}{i}, 'Fixation')
                counter_current_time = counter_current_time + duration_object_appearance;
            else
                for c = 1:num_categories
                    if strcmp(data{2}{i}, names{c})
                        onsets{c} = [onsets{c}, counter_current_time];
                        durations{c} = [durations{c}, duration_object_appearance];
                        counter_current_time = counter_current_time + duration_object_appearance;
                    end
                end
            end
        end
    end
end

save(output_mat_filename, 'names', 'onsets', 'durations');

