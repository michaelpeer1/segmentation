function segmentation_create_SPM_design_jrd_adaptation_grmp_overlay(psychopy_log_filename, output_mat_filename)
% Create segmentation adaptation design
% Regressors:
% 1. Different segment than before (1 - different, 0 - similar)
% 2. Same segment as before (1 - similar, 0 - different)
% 3. Repeat of the same stimulus
% Parametric modulation of regressors 1,2 - Euclidean distance between stimulus and previous on, UNDER THE SCHEMATIZATION MODEL, SEPARATELY FOR WITHIN-SEGMENT AND BETWEEN-SEGMENT PAIRS

num_objects = 16;
duration_object_appearance = 5;
duration_start_rest = 10;

% Reading the log file
temp_file = fopen(psychopy_log_filename, 'r');
data = textscan(temp_file, '%s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'Delimiter',',','HeaderLines',3);     % Read the first 12 columns of the log file 
fclose(temp_file);

% Initializing arrays
names = {'diff_segment', 'same_segment', 'repeating_stimuli'};
onsets = cell(1, length(names));
durations = cell(1, length(names));

% Defining the first stimulus as a repeating one, to model it and make sure
% there is something in the repeated stimuli regressor
onsets{3} = [onsets{3}, duration_start_rest];
durations{3} = [durations{3}, duration_object_appearance];


% Distances between objects, for parametric modulation
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
locations_all_objects(5:12,2) = locations_all_objects(5:12,2) + 75;                     % Overlay model
locations_distances_dsm = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
locations_distances_dsm = locations_distances_dsm / max(locations_distances_dsm(:));    % Normalizing to range 0-1

all_distances_within = [];
all_distances_between = [];
counter_current_time = duration_start_rest + duration_object_appearance + str2double(data{12}{1});
for i = 2:length(data{1})     % Going over each line of the log file
    if ~isempty(data{1}{i})
        current_stimulus_num = str2double(data{9}{i});     % The unity index of the stimulus
        previous_stimulus_num = str2double(data{9}{i-1});     % The unity index of the previous stimulus
        if current_stimulus_num ~= previous_stimulus_num
            segment_current_stim = (current_stimulus_num>=4 & current_stimulus_num<=11);
            segment_previous_stim = (previous_stimulus_num>=4 & previous_stimulus_num<=11);
            if segment_current_stim ~= segment_previous_stim    % Define first regressor - is the segment different than the previous one
                onsets{1} = [onsets{1}, counter_current_time];
                durations{1} = [durations{1}, duration_object_appearance];
                all_distances_between = [all_distances_between, locations_distances_dsm(str2double(data{9}{i}) + 1, str2double(data{9}{i-1}) + 1)];
            else        % Define second regressor - same segment as the previous one
                onsets{2} = [onsets{2}, counter_current_time];
                durations{2} = [durations{2}, duration_object_appearance];
                all_distances_within = [all_distances_within, locations_distances_dsm(str2double(data{9}{i}) + 1, str2double(data{9}{i-1}) + 1)];
            end
        else
            % Define third regressor - repeating stimuli
            onsets{3} = [onsets{3}, counter_current_time];
            durations{3} = [durations{3}, duration_object_appearance];
        end
        counter_current_time = counter_current_time + duration_object_appearance + str2double(data{12}{i});
    end

end

% Parametric modulation by distance
all_distances_between = all_distances_between - nanmean(all_distances_between);        % Demeaning the parameter
all_distances_within = all_distances_within - nanmean(all_distances_within);        % Demeaning the parameter
pmod = struct('name', {''}, 'param', {}, 'poly', {});
pmod(1).name{1} = 'distance_between';
pmod(1).poly{1} = 1;
pmod(1).param{1} = all_distances_between;
pmod(2).name{1} = 'distance_within';
pmod(2).poly{1} = 1;
pmod(2).param{1} = all_distances_within;

save(output_mat_filename, 'names', 'onsets', 'durations', 'pmod');
