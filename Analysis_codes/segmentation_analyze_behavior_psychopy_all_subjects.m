[~, ~, subjects_file] = xlsread('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/fMRI_experiment_logfiles/Subjects_file_fMRI.xlsx');
num_subjects = sum(cellfun(@length, subjects_file(4:end,1)) > 1);
filenames_objectviewing_day1 = {};
filenames_objectviewing_day2 = {};
filenames_localizer = {};
filenames_JRD = {};

for i=1:num_subjects
    if ~isnan(subjects_file{3+i, 14})       % fMRI version of the task
        filenames_localizer{i} = {[fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 15}) '.csv'], [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 16}) '.csv']};
        filenames_objectviewing_day1{i} = {[fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 17}) '.csv'], [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 18}) '.csv']};
        filenames_objectviewing_day2{i} = {[fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 19}) '.csv'], [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 20}) '.csv']};
        if ~isnan(subjects_file{3+i, 23})
            filenames_JRD{i} = {[fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 21}) '.csv'], [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 22}) '.csv'], [fullfile(subjects_file{3+i, 27}, subjects_file{3+i, 23}) '.csv']};
        end
    end
end


all_subjects_results_localizer = struct(); 
all_subjects_results_localizer.num_overall = [];
all_subjects_results_localizer.num_correct = [];
all_subjects_results_localizer.false_positives = [];

all_subjects_results_objectviewing_day2 = struct();
all_subjects_results_objectviewing_day2.num_correct = [];
all_subjects_results_objectviewing_day2.num_overall = [];
all_subjects_results_objectviewing_day2.RTs_same_segment = [];
all_subjects_results_objectviewing_day2.RTs_different_segment = [];
all_subjects_results_objectviewing_day2.RTs_same_orth_segment = [];
all_subjects_results_objectviewing_day2.RTs_different_orth_segment = [];
all_subjects_results_objectviewing_day2.corr_distances_RTs = [];
all_subjects_results_objectviewing_day2.corr_distances_RTs_within_seg_adj_quad = [];
all_subjects_results_objectviewing_day2.corr_distances_RTs_within_orth_seg_adj_quad = [];

all_subjects_results_objectviewing_day1 = struct();
all_subjects_results_objectviewing_day1.num_correct = [];
all_subjects_results_objectviewing_day1.num_overall = [];
all_subjects_results_objectviewing_day1.RTs_same_segment = [];
all_subjects_results_objectviewing_day1.RTs_different_segment = [];
all_subjects_results_objectviewing_day1.RTs_same_orth_segment = [];
all_subjects_results_objectviewing_day1.RTs_different_orth_segment = [];
all_subjects_results_objectviewing_day1.corr_distances_RTs = [];
all_subjects_results_objectviewing_day1.corr_distances_RTs_within_seg_adj_quad = [];
all_subjects_results_objectviewing_day1.corr_distances_RTs_within_orth_seg_adj_quad = [];

all_subjects_results_JRD = struct();
all_subjects_results_JRD.num_overall = [];
all_subjects_results_JRD.num_answered = [];
all_subjects_results_JRD.num_correct = [];
all_subjects_results_JRD.num_correct_within_segment = [];
all_subjects_results_JRD.num_correct_between_segments = [];
all_subjects_results_JRD.RTs_same_segment = [];
all_subjects_results_JRD.RTs_different_segment = [];
all_subjects_results_JRD.RTs_same_orth_segment = [];
all_subjects_results_JRD.RTs_different_orth_segment = [];
all_subjects_results_JRD.match_objects_start_quadrant = [];
all_subjects_results_JRD.match_objects_all = [];
all_subjects_results_JRD.angles_accuracy = [];
all_subjects_results_JRD.angles_num = [];
all_subjects_results_JRD.responses_per_object = [];
all_subjects_results_JRD.corr_distances_RTs = [];
all_subjects_results_JRD.corr_distances_RTs_within_seg_adj_quad = [];
all_subjects_results_JRD.corr_distances_RTs_within_orth_seg_adj_quad = [];

for subj = 1:num_subjects
    disp(subj)
    all_subjects_results_localizer(subj) = analyze_segmentation_psychopy_localizer(filenames_localizer{subj});
    all_subjects_results_objectviewing_day1(subj) = analyze_segmentation_psychopy_objectviewing(filenames_objectviewing_day1{subj});
    all_subjects_results_objectviewing_day2(subj) = analyze_segmentation_psychopy_objectviewing(filenames_objectviewing_day2{subj});
    if ~isempty(filenames_JRD{subj})
        all_subjects_results_JRD(subj) = analyze_segmentation_psychopy_JRD(filenames_JRD{subj});
    end
end


% Saving the object co-viewing in JRD matrix
for i=1:num_subjects
    if ~isnan(subjects_file{3+i, 14})       % fMRI version of the task
        if ~isempty(filenames_JRD{i})
            current_subject_JRD_analysis_dir = fullfile(fullfile('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data', subjects_file{3+i, 29}), '/Analysis/Analysis_JRD');
            current_subject_coviewing_mat = all_subjects_results_JRD(i).match_objects_all;
            save(fullfile(current_subject_JRD_analysis_dir, 'object_coviewing_JRD.mat'), 'current_subject_coviewing_mat');
        end
    end
end


% Saving the responses similarity matrix
for i=1:num_subjects
    if ~isnan(subjects_file{3+i, 14})       % fMRI version of the task
        if ~isempty(filenames_JRD{i})
            % Normalizing the responses to percent of answers for each specific direction
            current_responses = all_subjects_results_JRD(i).responses_per_object;
            current_num_responses_per_object = sum(current_responses, 2);
            current_responses = [current_responses(:,1) ./ current_num_responses_per_object, current_responses(:,2) ./ current_num_responses_per_object];
            % Looking at the difference between left and right (in percent responses)
            current_responses = current_responses(:,1) - current_responses(:,2);
            % Calculating the difference between the responses
            current_responses_diff = pdist2(current_responses, current_responses);
            % Saving
            current_subject_JRD_analysis_dir = fullfile(fullfile('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data', subjects_file{3+i, 29}), '/Analysis/Analysis_JRD');
            save(fullfile(current_subject_JRD_analysis_dir, 'responses_similarity_JRD.mat'), 'current_responses_diff');
        end
    end
end


% Analyzing the JRD accuracy by angle
a = cell2mat({all_subjects_results_JRD.angles_accuracy});
b = cell2mat({all_subjects_results_JRD.angles_num});
c = sum(a,2)./sum(b,2);
d = nanmean(a./b,2);
lab={}; for i=1:360,lab{i}=num2str(i);end

plot(c(~isnan(c)))
xticks(1:length(find(~isnan(c))))
xticklabels(lab(~isnan(c)))
xtickangle(90)

e = []; for i=1:12, e(i) = nanmean(c((i-1)*30+1:i*30));end
