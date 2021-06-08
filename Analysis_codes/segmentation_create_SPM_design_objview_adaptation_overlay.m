function segmentation_create_SPM_design_objview_adaptation_overlay(psychopy_log_filename, output_mat_filename)
% Create segmentation adaptation design
% Regressors:
% 1. Repeat of the same stimulus
% 2. Stimulus after null trial
% 3. All stimuli appearances that are not null, repeats, or after a null
% Parametric modulation of regressor 3 - Euclidean distance between stimulus and previous one, UNDER THE SCHEMATIZATION (OVERLAY) MODEL

duration_object_appearance = 1;
duration_object_rest = 1;
duration_start_rest = 10;

% Reading the log file
temp_file = fopen(psychopy_log_filename, 'r');
data = textscan(temp_file, '%s %s %s %s %s %f %*[^\n]', 'Delimiter',',','HeaderLines',3);     % Read the first four columns of the log file (stimulus filename, stimulus index, stimulus name, stimulus unity position)
fclose(temp_file);

% Initializing arrays
names = {'repeat_stimulus', 'stim_after_null', 'stims_non_repeat_null'};
onsets = cell(1, length(names));
durations = cell(1, length(names));

for i = 2:length(data{1})     % Going over each line of the log file
    if ~isempty(data{2}{i})
        if str2double(data{2}{i}) ~= 17      % 17 is the sign for blank
            current_stimulus_num = str2double(data{4}{i});     % The unity index of the stimulus
            if str2double(data{2}{i-1}) ~= 17   % If the previous stimulus is not null
                previous_stimulus_num = str2double(data{4}{i-1});       % The unity index of the previous stimulus
                if current_stimulus_num == previous_stimulus_num        % Define first regressor - is it the same stimulus as the previous one
                    onsets{1} = [onsets{1}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
                    durations{1} = [durations{1}, (duration_object_appearance + duration_object_rest)];
                end
            else    % Defining the second regressor - if the stimulus comes after a null trial    
                onsets{2} = [onsets{2}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
                durations{2} = [durations{2}, (duration_object_appearance + duration_object_rest)];
            end
        end
    end
end

% Distances between objects, for parametric modulation
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
locations_all_objects(5:12,2) = locations_all_objects(5:12,2) + 75;             % Overlay model
locations_distances_dsm = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
locations_distances_dsm = locations_distances_dsm / max(locations_distances_dsm(:));    % Normalizing to range 0-1

all_distances = [];
for i = 2:length(data{1})
    if (~isempty(data{2}{i}) && str2double(data{2}{i}) ~= 17 && str2double(data{2}{i-1}) ~= 17 && str2double(data{4}{i}) ~= str2double(data{4}{i-1}))
        % Defining the third regressors - appearance of stimuli that are not after null, null, or repeats (for parametric modulation)
        onsets{3} = [onsets{3}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
        durations{3} = [durations{3}, (duration_object_appearance + duration_object_rest)];
        % Calculating the distance between the previous and current object
        all_distances = [all_distances, locations_distances_dsm(str2double(data{4}{i}) + 1, str2double(data{4}{i-1}) + 1)];
    end
end
all_distances = all_distances - nanmean(all_distances);        % Demeaning the parameter

pmod = struct('name', {''}, 'param', {}, 'poly', {});
pmod(3).name{1} = 'distance';
pmod(3).poly{1} = 1;
pmod(3).param{1} = all_distances;


save(output_mat_filename, 'names', 'onsets', 'durations', 'pmod');
