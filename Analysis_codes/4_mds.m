% MDS and neural correlations

% Defining general variables
betas_file1 = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data/sub-seg01/Analysis/Analysis_objectviewing/Run1/beta_0001.nii';       % Reading one betas file for reslicing of all images to this file's resolution
data_parent_dir = '/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/1 - Segmentation/Segmentation_data';
subj_data_dirs = dir(fullfile(data_parent_dir, 'sub-seg*'));
num_conditions = 16;

all_ROI_names = {'bilat_RSC_activation_jrd', 'bilat_OPA_activation_jrd', 'lHCant'};
num_ROIs = length(all_ROI_names);


% Defining variables
allsubjs_pattern_dist_objectviewing_day2 = cell(size(all_ROI_names));
allsubjs_pattern_dist_objectviewing_day2minus1 = cell(size(all_ROI_names));
allsubjs_pattern_dist_jrd = cell(size(all_ROI_names));
for i=1:length(all_ROI_names)
    allsubjs_pattern_dist_objectviewing_day2{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
    allsubjs_pattern_dist_objectviewing_day2minus1{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
    allsubjs_pattern_dist_jrd{i} = nan(num_conditions,num_conditions,length(subj_data_dirs));
end

% Getting the neural patterns correlation (distance) matrix in each ROI
for subj = 1:length(subj_data_dirs)
    % Defining current subject directories
    disp(subj)
    % Object viewing analysis dirs
    curr_subj_analysis_dir = fullfile(fullfile(data_parent_dir, subj_data_dirs(subj).name), 'Analysis');
    curr_subj_objview_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_objectviewing');
    curr_subj_objview_analysis_dirs = {fullfile(curr_subj_objview_analysis_dir, 'Run1'),...
        fullfile(curr_subj_objview_analysis_dir, 'Run2'),...
        fullfile(curr_subj_objview_analysis_dir, 'Run3'),...
        fullfile(curr_subj_objview_analysis_dir, 'Run4')};
    num_runs_objview = length(curr_subj_objview_analysis_dirs);
    % JRD analysis dirs
    curr_subj_jrd_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD');
    curr_subj_jrd_analysis_dirs = {fullfile(curr_subj_jrd_analysis_dir, 'Run1'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run2'),...
        fullfile(curr_subj_jrd_analysis_dir, 'Run3')};
    curr_subj_jrd_adaptation_analysis_dir = fullfile(curr_subj_analysis_dir, 'Analysis_JRD_adaptation');
    num_runs_jrd = length(dir(fullfile(curr_subj_jrd_adaptation_analysis_dir, '*confounds.mat')));
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
            % Getting the ROI xyz coordinates
            [curr_ROI_coordx,curr_ROI_coordy,curr_ROI_coordz] = ind2sub(size(all_ROIs_individual{roi}),curr_ROI);
            
            all_betas_ROI_objview = nan(length(curr_ROI), num_conditions, num_runs_objview);
            for run = 1:num_runs_objview
                % Get the already noise-normalized betas (obtained with the temp_noise_normalize_betas script)
                norm_beta_files_current = dir(fullfile(curr_subj_objview_analysis_dirs{run}, 'normalized_beta*.nii'));
                for b = 1:num_conditions
                    curr_norm_beta_image = spm_sample_vol(spm_vol(fullfile(curr_subj_objview_analysis_dirs{run}, norm_beta_files_current(b).name)), curr_ROI_coordx, curr_ROI_coordy, curr_ROI_coordz, 0);
                    all_betas_ROI_objview(:,b,run) = curr_norm_beta_image;
                end
            end
            % Removing the cocktail mean from each day (mean activity at each voxel across patterns)
            mean_run_pattern = nanmean(all_betas_ROI_objview(:,:,run), 2);
            for b = 1:num_conditions
                all_betas_ROI_objview(:,b,run) = all_betas_ROI_objview(:,b,run) - mean_run_pattern;
            end
            % Getting the average pattern across runs instead, and calculating the correlation between these average patterns
            all_betas_day1 = nanmean(all_betas_ROI_objview(:,:,1:2), 3);
            all_betas_day2 = nanmean(all_betas_ROI_objview(:,:,3:4), 3);
            corr_ROI_objview_day1 = fisherz(corr(all_betas_day1, all_betas_day1, 'rows', 'complete'));
            corr_ROI_objview_day2 = fisherz(corr(all_betas_day2, all_betas_day2, 'rows', 'complete'));
            
            % Pattern distances
            pattern_dist_objectviewing_day1 = 1 - corr_ROI_objview_day1;
            pattern_dist_objectviewing_day2 = 1 - corr_ROI_objview_day2;
            allsubjs_pattern_dist_objectviewing_day2{roi}(:,:,subj) = pattern_dist_objectviewing_day2;
            allsubjs_pattern_dist_objectviewing_day2minus1{roi}(:,:,subj) = pattern_dist_objectviewing_day2 - pattern_dist_objectviewing_day1;
            
            
            % JRD
            if num_runs_jrd > 0
                all_betas_ROI_jrd = nan(length(curr_ROI), num_conditions, num_runs_jrd);
                for run = 1:num_runs_jrd
                    
                    % Get the already noise-normalized betas (obtained with the temp_noise_normalize_betas script)
                    norm_beta_files_current = dir(fullfile(curr_subj_jrd_analysis_dirs{run}, 'normalized_beta*.nii'));
                    for b = 1:num_conditions
                        curr_norm_beta_image = spm_sample_vol(spm_vol(fullfile(curr_subj_jrd_analysis_dirs{run}, norm_beta_files_current(b).name)), curr_ROI_coordx, curr_ROI_coordy, curr_ROI_coordz, 0);
                        all_betas_ROI_jrd(:,b,run) = curr_norm_beta_image;
                    end
                    
                    % Removing the cocktail mean from each day (mean activity at each voxel across patterns)
                    mean_run_pattern = nanmean(all_betas_ROI_jrd(:,:,run), 2);
                    for b = 1:num_conditions
                        all_betas_ROI_jrd(:,b,run) = all_betas_ROI_jrd(:,b,run) - mean_run_pattern;
                    end
                end
                
                % Getting the average pattern across runs instead, and calculating the correlation between these average patterns
                all_betas_jrd = nanmean(all_betas_ROI_jrd, 3);
                corr_ROI_jrd = fisherz(corr(all_betas_jrd, all_betas_jrd, 'rows', 'complete'));
                
                if ~isempty(corr_ROI_jrd)                    
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



%% MDS
curr_mat = allsubjs_mean_pattern_dist_jrd;
roi = 1;

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
% figure;plot(Y(:,1),Y(:,2),'*'), text(Y(:,1),Y(:,2),nums), title(['roi ' all_ROI_names{roi}]);
% hold on; scatter(l1(:,1),l1(:,2),'r'), text(l1(:,1),l1(:,2),nums,'r')
% diffs = Y-l1; quiver(l1(:,1),l1(:,2),diffs(:,1),diffs(:,2),0)
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
