function current_subj_results = analyze_segmentation_psychopy_localizer(logfiles)
% current_subj_results = analyze_segmentation_psychopy_localizer(logfiles)
%
% Analyzing the results from the functional localizer logfiles (success
% rate in the one-back task)

current_subj_num_correct = 0;
current_subj_num_overall = 0;
current_subj_false_positives = 0;

for f = 1:length(logfiles)
    % Loading the data
    %     [~,~,results] = xlsread(logfiles{f});
    fid = fopen(logfiles{f});
    num_columns = 26; temp_string = ''; for i=1:num_columns, temp_string = [temp_string '%s ']; end
    data = textscan(fid, temp_string, 'delimiter', ',');
    fclose(fid);
    results = {}; for i=1:length(data), results(:,i) = data{i}; end
    for i = 1:length(results(:))        % Changing strings to numbers
        if ~isnan(str2double(results{i}))
            results{i} = str2double(results{i});
        end
    end

    all_stimuli_names = results(3:end, 1);
    all_responses = results(3:end, find(strcmp(results(1,:), 'response.keys')));
    
    % Calculating success
    for i = 2:length(all_stimuli_names)
        if (strcmp(all_stimuli_names{i}, all_stimuli_names{i-1}) && ~strcmp(all_stimuli_names{i}, 'None'))   % If the same stimulus appears twice and is not fixation (subject should press something in this case)
            current_subj_num_overall = current_subj_num_overall + 1;
            if ~strcmp(all_responses{i}, 'None')        % If the subject pressed any key when this happened
                current_subj_num_correct = current_subj_num_correct + 1;
            end
        elseif ~strcmp(all_responses{i}, 'None')        % If the subject pressed any key, but shouldn't have (the stimulus didn't repeat the previous one)
            current_subj_false_positives = current_subj_false_positives + 1;
        end
    end
end

current_subj_results.num_overall = current_subj_num_overall;
current_subj_results.num_correct = current_subj_num_correct;
current_subj_results.false_positives = current_subj_false_positives;
