function segmentation_create_SPM_design_objview_adaptation_group_remap(psychopy_log_filename, output_mat_filename)
% Create segmentation adaptation design
% Regressors:
% 1. Different segment than before (1 - different, 0 - similar)
% 2. Same segment as before (1 - similar, 0 - different)
% 3. Repeat of the same stimulus
% 4. Stimulus after null trial
% Parametric modulation of regressors 1,2 - Euclidean distance between stimulus and previous on, SEPARATELY FOR WITHIN-SEGMENT AND BETWEEN-SEGMENT PAIRS

duration_object_appearance = 1;
duration_object_rest = 1;
duration_start_rest = 10;

% Reading the log file
temp_file = fopen(psychopy_log_filename, 'r');
data = textscan(temp_file, '%s %s %s %s %s %f %*[^\n]', 'Delimiter',',','HeaderLines',3);     % Read the first four columns of the log file (stimulus filename, stimulus index, stimulus name, stimulus unity position)
fclose(temp_file);

% Initializing arrays
names = {'diff_segment', 'same_segment', 'repeat_stimulus', 'stim_after_null'};
onsets = cell(1, length(names));
durations = cell(1, length(names));



% Distances between objects, for parametric modulation
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
locations_distances_dsm = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
locations_distances_dsm = locations_distances_dsm / max(locations_distances_dsm(:));    % Normalizing to range 0-1

all_distances_within = [];
all_distances_between = [];
for i = 2:length(data{1})     % Going over each line of the log file
    if ~isempty(data{2}{i})
        if str2double(data{2}{i}) ~= 17      % 17 is the sign for blank
            current_stimulus_num = str2double(data{4}{i});     % The unity index of the stimulus            
            if str2double(data{2}{i-1}) ~= 17   % If the previous stimulus is not null
                previous_stimulus_num = str2double(data{4}{i-1});       % The unity index of the previous stimulus
                if current_stimulus_num ~= previous_stimulus_num
                    segment_current_stim = (current_stimulus_num>=4 & current_stimulus_num<=11);
                    segment_previous_stim = (previous_stimulus_num>=4 & previous_stimulus_num<=11);
                    if segment_current_stim ~= segment_previous_stim        % Define first regressor - is the segment different than the previous one
                        onsets{1} = [onsets{1}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
                        durations{1} = [durations{1}, (duration_object_appearance + duration_object_rest)];
                        all_distances_between = [all_distances_between, locations_distances_dsm(str2double(data{4}{i}) + 1, str2double(data{4}{i-1}) + 1)];        % Calculating the distance between the previous and current object
                    else            % Define second regressor - is the segment different than the previous one
                        onsets{2} = [onsets{2}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
                        durations{2} = [durations{2}, (duration_object_appearance + duration_object_rest)];
                        all_distances_within = [all_distances_within, locations_distances_dsm(str2double(data{4}{i}) + 1, str2double(data{4}{i-1}) + 1)];        % Calculating the distance between the previous and current object
                    end
                else    % Define third regressor - is it the same stimulus as the previous one
                    onsets{3} = [onsets{3}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
                    durations{3} = [durations{3}, (duration_object_appearance + duration_object_rest)];
                end
            else        % Defining the fourth regressor - if the stimulus comes after a null trial
                onsets{4} = [onsets{4}, duration_start_rest + (sum(data{6}(1:i)) - 1) * 2];
                durations{4} = [durations{4}, (duration_object_appearance + duration_object_rest)];
            end
        end
    end
end

% Parametric modulation by distance from previous object
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
