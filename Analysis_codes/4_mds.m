% MDS and neural correlations

% Defining general variables
betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data/sub-seg01/Analysis/Analysis_objview/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/Segmentation project/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview'};
%     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC','rHC','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview'};
% all_ROI_names = {'bilat_OPA', 'rOPA', 'bilat_RSC', 'rRSC', 'bilat_OPA_activation_jrd', 'rOPA_activation_jrd', 'bilat_RSC_activation_jrd', 'rRSC_activation_jrd', 'bilat_HC_activation_objview', 'lHC_activation_objview'};
% all_ROI_names = {'bilat_HC_activation_jrd', 'bilat_HC_activation_objview', 'bilat_RSC', 'bilat_PPA', 'bilat_OPA', 'bilat_LOC'};
% % all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_PPA_activation_jrd', 'bilat_OPA_activation_jrd', 'bilat_HC_activation_jrd','bilat_RSC_activation_objview', 'bilat_PPA_activation_objview', 'bilat_OPA_activation_objview', 'bilat_HC_activation_objview'};
% %     all_ROI_names = {'lRSC_activation_jrd', 'rRSC_activation_jrd', 'lPPA_activation_jrd', 'rPPA_activation_jrd', 'lOPA_activation_jrd', 'rOPA_activation_jrd', 'lHC_activation_jrd','rHC_activation_jrd','lRSC_activation_objview', 'rRSC_activation_objview', 'lPPA_activation_objview', 'rPPA_activation_objview', 'lOPA_activation_objview', 'rOPA_activation_objview', 'lHC_activation_objview','rHC_activation_objview'};
num_ROIs = length(all_ROI_names);

% Define locations of objects and create the distance matrix
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
mat_distances = pdist2(locations_all_objects, locations_all_objects);         % Distances between all objects
mat_distances = mat_distances / max(mat_distances(:));    % Normalizing to range 0-1
% Create the segment identity similarity matrix
mat_segments = ones(num_conditions);
mat_segments(1:4,1:4) = 0; mat_segments(13:16,13:16) = 0; mat_segments(5:12,5:12) = 0; mat_segments(13:16,1:4) = 0; mat_segments(1:4,13:16) = 0;
% Create the quadrant identity similarity matrix
mat_quadrants = ones(num_conditions);
mat_quadrants(1:4,1:4) = 0; mat_quadrants(5:8,5:8) = 0; mat_quadrants(9:12,9:12) = 0; mat_quadrants(13:16,13:16) = 0;
% Create the orthogonal segment identity similarity matrix
mat_orth_segments = ones(num_conditions);
mat_orth_segments(1:8,1:8) = 0; mat_orth_segments(9:16,9:16) = 0;
% Converting all matrices from dissimilarity to similarity matrices
mat_distances = 1 - mat_distances; mat_segments = 1 - mat_segments; mat_quadrants = 1 - mat_quadrants; mat_orth_segments = 1 - mat_orth_segments;
% Distance within segment, adjacent quadrants only
mat_dist_within_segment_adj_quadrants = mat_distances;
mat_dist_within_segment_adj_quadrants(mat_segments==0) = nan;
mat_dist_within_segment_adj_quadrants(mat_quadrants==1) = nan;
% Distance between segments, adjacent quadrants only
mat_diagonal_quadrants = zeros(16); mat_diagonal_quadrants(1:4,9:12) = 1; mat_diagonal_quadrants(9:12, 1:4) = 1; mat_diagonal_quadrants(5:8,13:16) = 1; mat_diagonal_quadrants(13:16,5:8) = 1;
mat_dist_between_segments_adj_quadrants = mat_distances;
mat_dist_between_segments_adj_quadrants(mat_segments==1) = nan;
mat_dist_between_segments_adj_quadrants(mat_diagonal_quadrants==1) = nan;
% Distance between objects in diagonal quadrants only
mat_dist_diagonal_quadrants = nan(16);
mat_dist_diagonal_quadrants(mat_diagonal_quadrants==1) = mat_distances(mat_diagonal_quadrants==1);
% Distances within quadrant only
mat_dist_within_quadrant = mat_distances;
mat_dist_within_quadrant(mat_quadrants==0) = nan;


allsubjs_pattern_dist_objectviewing_day2 = cell(size(all_ROI_names));
allsubjs_pattern_dist_objectviewing_day2minus1 = cell(size(all_ROI_names));
allsubjs_pattern_dist_jrd = cell(size(all_ROI_names));
for i=1:length(all_ROI_names)
    allsubjs_pattern_dist_objectviewing_day2{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
    allsubjs_pattern_dist_objectviewing_day2minus1{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
    allsubjs_pattern_dist_jrd{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
end

all_neural_corr_ROI_jrd = nan(num_conditions, num_conditions, num_ROIs, length(subj_data_dirs));
all_neural_corr_objview_day2 = nan(num_conditions, num_conditions, num_ROIs, length(subj_data_dirs));
all_neural_corr_objview_d2minusd1 = nan(num_conditions, num_conditions, num_ROIs, length(subj_data_dirs));

% Calculation RSA in each subject and ROI
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    % Object viewing analysis dirs
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_objectviewing_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objectviewing_analysis_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run2'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Run4')};
    num_runs_objectviewing = length(curr_subj_objectviewing_analysis_dirs);
    curr_subj_objectviewing_combined_dirs = {fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day1'),...
        fullfile(curr_subj_objectviewing_analysis_dir, 'Combined_runs_day2')};
    % JRD analysis dirs
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    num_runs_jrd = length(dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat')));
    curr_subj_jrd_combined_dir = fullfile(curr_subj_jrd_analysis_dir, 'Combined_runs');
    % Functional localizer analysis dirs
    curr_subj_localizer_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_localizer');
    
    
    %% Loading the subject specific ROIs
    all_ROIs_individual = {};
    for i = 1:length(all_ROI_names)
        all_ROIs_individual{i} = spm_read_vols(spm_vol(fullfile(curr_subj_localizer_analysis_dir, [all_ROI_names{i}, '_individual.nii'])));
    end
    num_ROIs_individual = length(all_ROIs_individual);
        
    
    %% MVPA analysis in each ROI
    for roi = 1:num_ROIs_individual        % Going over ROIs
        curr_ROI = find(all_ROIs_individual{roi}(:)==1);
        if ~isempty(curr_ROI)
            % Reading object viewing data in the ROI
            all_tvalues_ROI_objview = nan(length(curr_ROI), num_conditions, length(curr_subj_objectviewing_combined_dirs));
            for r = 1:length(curr_subj_objectviewing_combined_dirs)
                t_values_files_current = dir(fullfile(curr_subj_objectviewing_combined_dirs{r}, 'spmT*.nii'));
                for b = 1:num_conditions
                    curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_objectviewing_combined_dirs{r}, t_values_files_current(b).name)));
                    all_tvalues_ROI_objview(:,b,r) = curr_tvalues_image(curr_ROI);
                end
            end
            % Removing the cocktail mean from each day (mean activity at each voxel across patterns)
            for r = 1:length(curr_subj_objectviewing_combined_dirs)
                mean_day_pattern = nanmean(all_tvalues_ROI_objview(:,:,r), 2);
                for b = 1:num_conditions
                    all_tvalues_ROI_objview(:,b,r) = all_tvalues_ROI_objview(:,b,r) - mean_day_pattern;
                end
            end
                        
            % Calculating correlation between patterns
            corr_ROI_objview_day1 = fisherz(corr(all_tvalues_ROI_objview(:,:,1), all_tvalues_ROI_objview(:,:,1), 'rows', 'complete'));                    % Correlation between day 1 patterns, averaged across runs
            corr_ROI_objview_day2 = fisherz(corr(all_tvalues_ROI_objview(:,:,2), all_tvalues_ROI_objview(:,:,2),'rows', 'complete'));                    % Correlation between day 2 patterns, averaged across runs
            
            % Neural distances
            all_neural_corr_objview_day2(:,:,roi,subj) = corr_ROI_objview_day2;
            all_neural_corr_objview_d2minusd1(:,:,roi,subj) = corr_ROI_objview_day2 - corr_ROI_objview_day1;
            
            % Pattern distances
            pattern_dist_objectviewing_day1 = 1 - corr_ROI_objview_day1;
            pattern_dist_objectviewing_day2 = 1 - corr_ROI_objview_day2;            
            allsubjs_pattern_dist_objectviewing_day2{roi}(:,:,subj) = pattern_dist_objectviewing_day2;
            allsubjs_pattern_dist_objectviewing_day2minus1{roi}(:,:,subj) = pattern_dist_objectviewing_day2 - pattern_dist_objectviewing_day1;
            
            
            % JRD
            if num_runs_jrd > 0
                % Extracting the t-values at each ROI - JRD
                all_tvalues_ROI_jrd = nan(length(curr_ROI), num_conditions);
                t_values_files_current = dir(fullfile(curr_subj_jrd_combined_dir, 'spmT*.nii'));
                for b = 1:num_conditions
                    curr_tvalues_image = spm_read_vols(spm_vol(fullfile(curr_subj_jrd_combined_dir, t_values_files_current(b).name)));
                    all_tvalues_ROI_jrd(:,b) = curr_tvalues_image(curr_ROI);
                end
                % Removing the cocktail mean (mean activity at each voxel across patterns)
                mean_jrd_pattern = nanmean(all_tvalues_ROI_jrd, 2);
                for b = 1:num_conditions
                    all_tvalues_ROI_jrd(:,b) = all_tvalues_ROI_jrd(:,b) - mean_jrd_pattern;
                end
                
                % Calculating correlation between patterns
                corr_ROI_jrd = fisherz(corr(all_tvalues_ROI_jrd, all_tvalues_ROI_jrd, 'rows', 'complete'));
                if ~isempty(corr_ROI_jrd)
                    % Neural corr
                    all_neural_corr_ROI_jrd(:,:,roi,subj) = corr_ROI_jrd;
                    
                    % Pattern dist
                    pattern_dist_jrd = 1 - corr_ROI_jrd;
                    allsubjs_pattern_dist_jrd{roi}(:,:,subj) = pattern_dist_jrd;
                end
            end
        end
    end
end


% Summarizing pattern distances
for roi=1:length(all_ROI_names)
    allsubjs_mean_pattern_dist_jrd{roi} = nanmean(allsubjs_pattern_dist_jrd{roi},3);
    allsubjs_mean_pattern_dist_objectviewing_day2{roi} = nanmean(allsubjs_pattern_dist_objectviewing_day2{roi},3);
    allsubjs_mean_pattern_dist_objectviewing_day2minus1{roi} = nanmean(allsubjs_pattern_dist_objectviewing_day2minus1{roi},3);
end

% Summarizing neural correlations
neural_corr_jrd_within_quadrant = nan(length(subj_data_dirs), num_ROIs);
neural_corr_jrd_within_segment_adj_quad = nan(length(subj_data_dirs), num_ROIs);
neural_corr_jrd_between_segments_adj_quad = nan(length(subj_data_dirs), num_ROIs);
neural_corr_jrd_diagonal_quadrants = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewday2_within_quadrant = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewday2_within_segment_adj_quad = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewday2_between_segments_adj_quad = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewday2_diagonal_quadrants = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewd2minusd1_within_quadrant = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewd2minusd1_within_segment_adj_quad = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewd2minusd1_between_segments_adj_quad = nan(length(subj_data_dirs), num_ROIs);
neural_corr_objviewd2minusd1_diagonal_quadrants = nan(length(subj_data_dirs), num_ROIs);
for subj = 1:length(subj_data_dirs)
    for roi=1:num_ROIs
        curr_neural_corr = all_neural_corr_ROI_jrd(:,:,roi,subj); curr_neural_corr(eye(size(curr_neural_corr)) == 1) = nan;
        neural_corr_jrd_within_quadrant(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_within_quadrant)));
        neural_corr_jrd_within_segment_adj_quad(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_within_segment_adj_quadrants)));
        neural_corr_jrd_between_segments_adj_quad(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_between_segments_adj_quadrants)));
        neural_corr_jrd_diagonal_quadrants(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_diagonal_quadrants)));
        curr_neural_corr = all_neural_corr_objview_day2(:,:,roi,subj); curr_neural_corr(eye(size(curr_neural_corr)) == 1) = nan;
        neural_corr_objviewday2_within_quadrant(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_within_quadrant)));
        neural_corr_objviewday2_within_segment_adj_quad(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_within_segment_adj_quadrants)));
        neural_corr_objviewday2_between_segments_adj_quad(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_between_segments_adj_quadrants)));
        neural_corr_objviewday2_diagonal_quadrants(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_diagonal_quadrants)));
        curr_neural_corr = all_neural_corr_objview_d2minusd1(:,:,roi,subj); curr_neural_corr(eye(size(curr_neural_corr)) == 1) = nan;
        neural_corr_objviewd2minusd1_within_quadrant(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_within_quadrant)));
        neural_corr_objviewd2minusd1_within_segment_adj_quad(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_within_segment_adj_quadrants)));
        neural_corr_objviewd2minusd1_between_segments_adj_quad(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_between_segments_adj_quadrants)));
        neural_corr_objviewd2minusd1_diagonal_quadrants(subj,roi) = nanmean(curr_neural_corr(~isnan(mat_dist_diagonal_quadrants)));
    end
end
m1=[nanmean(neural_corr_jrd_within_quadrant); nanmean(neural_corr_jrd_within_segment_adj_quad);nanmean(neural_corr_jrd_between_segments_adj_quad); nanmean(neural_corr_jrd_diagonal_quadrants)];
m2=[nanmean(neural_corr_objviewd2minusd1_within_quadrant); nanmean(neural_corr_objviewd2minusd1_within_segment_adj_quad);nanmean(neural_corr_objviewd2minusd1_between_segments_adj_quad); nanmean(neural_corr_objviewday2_diagonal_quadrants)];
m3=[nanmean(neural_corr_objviewday2_within_quadrant); nanmean(neural_corr_objviewday2_within_segment_adj_quad);nanmean(neural_corr_objviewday2_between_segments_adj_quad); nanmean(neural_corr_objviewday2_diagonal_quadrants)];



%% MDS
curr_mat = allsubjs_mean_pattern_dist_jrd;
roi = 3;

curr_mat{roi}(find(eye(size(curr_mat{roi})))) = 0;      % Zeroing the diagonal
curr_mat{roi} = curr_mat{roi} - min(curr_mat{roi}(:));     % Normalizing so all values will be positive 
curr_mat{roi}(find(eye(size(curr_mat{roi})))) = 0;      % Zeroing the diagonal
curr_mat{roi} = curr_mat{roi} / max(curr_mat{roi}(:));  % Normalizing to range 0-1

Y = cmdscale(curr_mat{roi},2);      % Multidimensional scaling

% Plotting the MDS, with gradient descent to fit it to actual distance locations as much as possible
nums={}; for i=1:16, nums{i}=num2str(i);end
locations_all_objects = [16 42; 50 58; 57 27; 25 17; ...
    42 -16; 58 -50; 27 -57; 17 -25; ...
    -16 -42; -50 -58; -57 -27; -25 -17; ...
    -42 16; -58 50; -27 57; -17 25];
l1 = locations_all_objects;
l1 = l1/150;
i = -90; R = [cosd(i) -sind(i); sind(i) cosd(i)]; l1 = l1*R; % Rotating to correct environment orientation
% Implementing a numerical gradient descent
param = {[0:360], [-0.02:0.01:0.02], [-0.02:0.01:0.02], [1,1.1], [-1 1]};    % Rotation degree, x-translation, y-translation, scaling, flipping
max_iters = 100;
max_dist_improvement = 0.0001;
curr_dist_improvement = 100;
iter = 1;
while iter <= max_iters && curr_dist_improvement > max_dist_improvement
    % Calculating the current distance
    curr_dist = pdist2(l1, Y);
    curr_dist = sum(diag(curr_dist).^2);    % sum of squared distances
    % Going over parameters to find the best ones
    all_dists = [];
    for deg = 1:length(param{1})     % Degree of rotation
        for xtrans = 1:length(param{2})    % X translation
            for ytrans = 1:length(param{3})  % Y translation
                for sc = 1:length(param{4})    % Scaling
                    for rot = 1:length(param{5})    % Flipping
                        curr_deg = param{1}(deg); curr_xtrans = param{2}(xtrans); curr_ytrans = param{3}(ytrans); curr_sc = param{4}(sc); curr_rot = param{5}(rot);
                        Y_new = Y;
                        Y_new(:,1) = curr_rot * Y_new(:,1);
                        R = [cosd(curr_deg) -sind(curr_deg); sind(curr_deg) cosd(curr_deg)];
                        Y_new = Y_new*R;
                        Y_new(:,1) = Y_new(:,1) + curr_xtrans;
                        Y_new(:,2) = Y_new(:,2) + curr_ytrans;
                        Y_new = Y_new * curr_sc;
                        dist_temp = pdist2(l1, Y_new);
                        all_dists(deg,xtrans,ytrans,sc,rot) = sum(diag(dist_temp).^2);
                    end
                end
            end
        end
    end
    [D,I] = min(all_dists(:));
    [I1,I2,I3,I4,I5] = ind2sub(size(all_dists),I);
    curr_deg = param{1}(I1); curr_xtrans = param{2}(I2); curr_ytrans = param{3}(I3); curr_sc = param{4}(I4); curr_rot = param{5}(I5);
    Y_new = Y;
    Y_new(:,1) = curr_rot * Y_new(:,1);
    R = [cosd(curr_deg) -sind(curr_deg); sind(curr_deg) cosd(curr_deg)];
    Y_new = Y_new*R;
    Y_new(:,1) = Y_new(:,1) + curr_xtrans;
    Y_new(:,2) = Y_new(:,2) + curr_ytrans;
    Y_new = Y_new * curr_sc;
    Y = Y_new;
    curr_dist_improvement = curr_dist - D;
    disp([iter, curr_dist_improvement])
    iter = iter + 1;
end
figure;plot(Y(:,1),Y(:,2),'*'), text(Y(:,1),Y(:,2),nums), title(['roi ' all_ROI_names{roi}]);
hold on; scatter(l1(:,1),l1(:,2),'r'), text(l1(:,1),l1(:,2),nums,'r')
diffs = Y-l1; quiver(l1(:,1),l1(:,2),diffs(:,1),diffs(:,2),0)
% figure;hold on;plot(Y(1:4,1),Y(1:4,2),'.g',Y(5:8,1),Y(5:8,2),'.b',Y(9:12,1),Y(9:12,2),'.r',Y(13:16,1),Y(13:16,2),'.y','MarkerSize',30), text(Y(:,1),Y(:,2),nums,'FontSize',12), title(['roi ' all_ROI_names{roi}]);
figure;hold on;
scatter(Y(1:2,1),Y(1:2,2),200,'g','p','filled')
scatter(Y(3:4,1),Y(3:4,2),100,'g','s','filled')
scatter(Y([5,8],1),Y([5,8],2),200,'b','p','filled')
scatter(Y([6,7],1),Y([6,7],2),100,'b','s','filled')
scatter(Y([11,12],1),Y([11,12],2),200,'r','p','filled')
scatter(Y([9,10],1),Y([9,10],2),100,'r','s','filled')
scatter(Y([14,15],1),Y([14,15],2),200,'y','p','filled')
scatter(Y([13,16],1),Y([13,16],2),100,'y','s','filled')
text(Y(:,1)+0.01,Y(:,2),nums,'FontSize',12), title(['roi ' all_ROI_names{roi}])
