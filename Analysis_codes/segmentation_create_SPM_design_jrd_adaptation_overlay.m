function segmentation_create_SPM_design_jrd_adaptation_overlay(psychopy_log_filename, output_mat_filename)
% Create segmentation adaptation design
% Regressors:
% 1. Repeat of the same stimulus
% 2. All stimuli appearances that are not repeats
% Parametric modulation of regressor 2 - Euclidean distance between stimulus and previous one, UNDER THE SCHEMATIZATION (OVERLAY) MODEL

duration_object_appearance = 5;
duration_start_rest = 10;

% Reading the log file
temp_file = fopen(psychopy_log_filename, 'r');
data = textscan(temp_file, '%s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'Delimiter',',','HeaderLines',3);     % Read the first 12 columns of the log file 
fclose(temp_file);

% Initializing arrays
names = {'repeating_stimuli', 'all_stim_non_repeat'};
onsets = cell(1, length(names));
durations = cell(1, length(names));

% Defining the first stimulus as a repeating one, to model it and make sure there is something in the repeated stimuli regressor
onsets{1} = [onsets{1}, duration_start_rest];
durations{1} = [durations{1}, duration_object_appearance];


counter_current_time = duration_start_rest + duration_object_appearance + str2double(data{12}{1});
for i = 2:length(data{1})     % Going over each line of the log file
    if ~isempty(data{1}{i})
        current_stimulus_num = str2double(data{9}{i});     % The unity index of the stimulus
        previous_stimulus_num = str2double(data{9}{i-1});     % The unity index of the previous stimulus
        if current_stimulus_num == previous_stimulus_num
            % Define first regressor - repeating stimuli
            onsets{1} = [onsets{1}, counter_current_time];
            durations{1} = [durations{1}, duration_object_appearance];
        else
            % Define second regressor - stimuli that are not repeats of previous ones (for distance parametric modulation)
            onsets{2} = [onsets{2}, counter_current_time];
            durations{2} = [durations{2}, duration_object_appearance];
        end
        counter_current_time = counter_current_time + duration_object_appearance + str2double(data{12}{i});
    end
end

% Distances between objects, for parametric modulation
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
locations_all_objects(5:12,2) = locations_all_objects(5:12,2) + 75;                     % Overlay model
locations_distances_dsm = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
locations_distances_dsm = locations_distances_dsm / max(locations_distances_dsm(:));    % Normalizing to range 0-1

all_distances = [];
for i = 2:length(data{1})
    if (~isempty(data{1}{i}) && str2double(data{9}{i}) ~= str2double(data{9}{i-1}))
        % Calculating the distance between the previous and current object
        all_distances = [all_distances, locations_distances_dsm(str2double(data{9}{i}) + 1, str2double(data{9}{i-1}) + 1)];
    end
end
all_distances = all_distances - nanmean(all_distances);        % Demeaning the parameter

pmod = struct('name', {''}, 'param', {}, 'poly', {});
pmod(2).name{1} = 'distance';
pmod(2).poly{1} = 1;
pmod(2).param{1} = all_distances;


save(output_mat_filename, 'names', 'onsets', 'durations', 'pmod');
