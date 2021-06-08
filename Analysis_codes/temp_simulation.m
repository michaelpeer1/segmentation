%% Modeling the MVPA results

use_Euclidean_distance = 0;     % Whether to use Euclidean or correlation distance in computation of neural RDMs

locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
% Normalizing locations to range 0-1
locations_all_objects(:,1) = locations_all_objects(:,1) - min(locations_all_objects(:,1));
locations_all_objects(:,2) = locations_all_objects(:,2) - min(locations_all_objects(:,2));
locations_all_objects(:,1) = locations_all_objects(:,1) / max(locations_all_objects(:,1));
locations_all_objects(:,2) = locations_all_objects(:,2) / max(locations_all_objects(:,2));

% Creating a distance model matrix (true relation between patterns)
mat_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
mat_distances = mat_distances / max(mat_distances(:));    % Normalizing to range 0-1
mat_distances = 1 - mat_distances;  % Converting to similarity matrix


% Getting the non-diagonal elements
num_objects = size(locations_all_objects, 1);
relevant_locs = find(triu(ones(num_objects),1));     % The non-diagonal matrix elements
mat_distances_nondiag = mat_distances(relevant_locs);


num_runs = 3;
num_ROI_voxels = 200;
num_subjects = 10;
num_permutations = 20;

for noise_level = 1:10

for p = 1:num_permutations
    for subj = 1:num_subjects
        %% Simulating the neural patterns according to the model matrix
        % Simulating neural patterns for the X and Y axis (location codes) - a
        % random pattern for the maximum and minimum of each axis; all other
        % patterns are a weighted average of these patterns, based on location
        location_x_min_pattern = rand(1, num_ROI_voxels); location_x_max_pattern = rand(1, num_ROI_voxels);
        location_y_min_pattern = rand(1, num_ROI_voxels); location_y_max_pattern = rand(1, num_ROI_voxels);
        % Getting the pattern for each object
        all_object_patterns = nan(num_objects, num_ROI_voxels, num_runs);
        for i = 1:num_objects
            curr_item_x = locations_all_objects(i,1);
            curr_item_y = locations_all_objects(i,2);
            curr_x_pattern = location_x_max_pattern * curr_item_x + location_x_min_pattern * (1 - curr_item_x);
            curr_y_pattern = location_y_max_pattern * curr_item_y + location_y_min_pattern * (1 - curr_item_y);
            for run=1:num_runs
                all_object_patterns(i,:,run) = mean([curr_x_pattern; curr_y_pattern]);
            end
        end
        % Adding noise to each run
        all_object_patterns_w_noise = all_object_patterns + rand(size(all_object_patterns)) * (0.5+noise_level/2);
        
        
        %% Simulating the RSA
        
        % averaging of patterns across runs
        mean_pattern_across_runs = mean(all_object_patterns_w_noise, 3)';
        if use_Euclidean_distance == 0
            all_patterns_corrmat = fisherz(corr(mean_pattern_across_runs, mean_pattern_across_runs));
        else
            all_patterns_corrmat = nan(num_objects, num_objects);
            for c1 = 1:size(all_patterns_corrmat, 1)
                for c2 = 1:size(all_patterns_corrmat, 2)
                    all_patterns_corrmat(c1,c2) = sqrt(sum((mean_pattern_across_runs(:,c1) - mean_pattern_across_runs(:,c2)).^2));
                end
            end
        end
        all_patterns_corrmat_nondiag = all_patterns_corrmat(relevant_locs);
        RSA_result_averaging(subj) = corr(all_patterns_corrmat_nondiag, mat_distances_nondiag, 'type', 'Spearman');
        all_subjs_patterns_RDM_averaging(subj,:) = all_patterns_corrmat_nondiag;
        
        % Cross-validation across runs
        current_corrs = nan(num_objects, num_objects, num_runs);
        counter = 1;
        for r1 = 1:num_runs
            for r2 = r1+1:num_runs
                % Going over all possible run combinations
                if use_Euclidean_distance == 0
                    current_corrs(:,:,counter) = fisherz(corr(all_object_patterns_w_noise(:,:,r1)', all_object_patterns_w_noise(:,:,r2)'));
                else
                    for c1 = 1:size(all_patterns_corrmat_nondiag, 1)
                        for c2 = 1:size(all_patterns_corrmat_nondiag, 2)
                            current_corrs(c1,c2,counter) = sqrt(sum((all_object_patterns_w_noise(:,c1,r1) - all_object_patterns_w_noise(:,c2,r2)).^2));
                        end
                    end
                end
                counter = counter + 1;
            end
        end
        all_patterns_corrmat = nanmean(current_corrs, 3);
        all_patterns_corrmat_nondiag = all_patterns_corrmat(relevant_locs);
        RSA_result_crossval(subj) = corr(all_patterns_corrmat_nondiag, mat_distances_nondiag, 'type', 'Spearman');
        all_subjs_patterns_RDM_crossval(subj,:) = all_patterns_corrmat_nondiag;
    end
    
    RSA_results_avg_allsubjs(p) = mean(RSA_result_averaging);
    RSA_results_cv_allsubjs(p) = mean(RSA_result_crossval);
    [~,RSA_pval_avg(p)] = ttest(RSA_result_averaging,[],'tail','right');
    [~,RSA_pval_cv(p)] = ttest(RSA_result_crossval,[],'tail','right');
    
    
    % Reliability analysis - correspondence between subjects' matrices
    reliability_averaging_current = corr(all_subjs_patterns_RDM_averaging,'type','Spearman');
    reliability_crossval_current = corr(all_subjs_patterns_RDM_crossval,'type','Spearman');
    reliability_averaging(p) = mean(reliability_averaging_current(find(triu(ones(num_subjects),1))));
    reliability_crossval(p) = mean(reliability_crossval_current(find(triu(ones(num_subjects),1))));
end

reliability_averaging_all(noise_level) = mean(reliability_averaging);
reliability_crossval_all(noise_level) = mean(reliability_crossval);
RSA_corr_avg_all(noise_level) = mean(RSA_results_avg_allsubjs);
RSA_corr_cv_all(noise_level) = mean(RSA_results_cv_allsubjs);
RSA_pval_avg_all(noise_level) = mean(RSA_pval_avg);
RSA_pval_cv_all(noise_level) = mean(RSA_pval_cv);

end
