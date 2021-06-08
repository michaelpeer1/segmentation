rsc_files = dir('C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\1 - Segmentation\Segmentation_data\RDMs from Marchette 2014\*RSC*.mat');
opa_files = dir('C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\1 - Segmentation\Segmentation_data\RDMs from Marchette 2014\*OPA*.mat');
ppa_files = dir('C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\1 - Segmentation\Segmentation_data\RDMs from Marchette 2014\*PPA*.mat');
hc_files = dir('C:\Users\michaelpeer1\Dropbox (Epstein Lab)\Epstein Lab Team Folder\Michael_Peer\1 - Segmentation\Segmentation_data\RDMs from Marchette 2014\*Hipp*.mat');
% rsc_files = dir('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/RDMs from Marchette 2014/*RSC*.mat');
% opa_files = dir('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/RDMs from Marchette 2014/*OPA*.mat');
% ppa_files = dir('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/RDMs from Marchette 2014/*PPA*.mat');
% hc_files = dir('/Users/mpeer/Dropbox (Epstein Lab)/Epstein Lab Team Folder/Michael_Peer/RDMs from Marchette 2014/*Hipp*.mat');

corrs_rsc = []; corrs_opa = []; corrs_ppa = []; corrs_hc = [];
for i=1:length(rsc_files)
    load(fullfile(rsc_files(i).folder,rsc_files(i).name));
    corrs_rsc(:,:,i) = correlation;
end
for i=1:length(opa_files)
    load(fullfile(opa_files(i).folder,opa_files(i).name));
    corrs_opa(:,:,i) = correlation;
end
for i=1:length(ppa_files)
    load(fullfile(ppa_files(i).folder,ppa_files(i).name));
    corrs_ppa(:,:,i) = correlation;
end
for i=1:length(hc_files)
    load(fullfile(hc_files(i).folder,hc_files(i).name));
    corrs_hc(:,:,i) = correlation;
end


%% Test correlation to different models
% Full environment
coordinates_all_objects = [3,25; 7,25; 10,22; 10,15; 7,12; 3,12; 0,15; 0,22; 27,7; 27,3; 24,0; 17,0; 14,3; 14,7; 17,10; 24,10];
mat_distances_full_model = pdist2(coordinates_all_objects,coordinates_all_objects);
mat_distances_full_model = mat_distances_full_model / max(mat_distances_full_model(:));    % Normalizing to range 0-1
mat_distances_full_model = 1 - mat_distances_full_model;        % Convert to similarity matrix
% Distance between objects - overlay model
mat_distances_overlay = mat_distances_full_model(1:8,1:8);
% coordinates = [3.5,16; 9.5,16; 13,13; 13,3; 9.5,0; 3.5,0; 0,3; 0,13];
% mat_distances_overlay = pdist2(coordinates,coordinates);
% mat_distances_overlay = mat_distances_overlay / max(mat_distances_overlay(:));    % Normalizing to range 0-1
% mat_distances_overlay = 1 - mat_distances_overlay;
mat_distances_overlay = [mat_distances_overlay, mat_distances_overlay;mat_distances_overlay,mat_distances_overlay];
% Within and between segments
mat_distances_within = mat_distances_full_model; mat_distances_within(1:8,9:16) = nan; mat_distances_within(9:16,1:8) = nan;
mat_distances_overlay_between = mat_distances_overlay; mat_distances_overlay_between(1:8,1:8) = nan; mat_distances_overlay_between(9:16,9:16) = nan;
mat_distances_full_model_between = mat_distances_full_model; mat_distances_full_model_between(1:8,1:8) = nan; mat_distances_full_model_between(9:16,9:16) = nan;

global_heading_all_objects = [0,0,90,90,180,180,270,270,90,90,180,180,270,270,0,0];
local_heading_all_objects = [0,0,90,90,180,180,270,270,0,0,90,90,180,180,270,270];
mat_global_heading = zeros(16); mat_local_heading = zeros(16);
for i = 1:16
    for j = 1:16
        if global_heading_all_objects(i) == global_heading_all_objects(j), mat_global_heading(i,j) = 1; end
        if abs(global_heading_all_objects(i) - global_heading_all_objects(j))==90, mat_global_heading(i,j) = 0.5; end        
        if local_heading_all_objects(i) == local_heading_all_objects(j), mat_local_heading(i,j) = 1; end
        if abs(local_heading_all_objects(i) - local_heading_all_objects(j))==90, mat_local_heading(i,j) = 0.5; end        
    end
end


all_mats = {corrs_rsc,corrs_ppa,corrs_opa,corrs_hc};
dist_corr_within = []; dist_corr_overlay = []; dist_corr_overlay_between = []; dist_corr_full = []; dist_corr_full_between = [];
for m = 1:length(all_mats)      % Going over ROIs
    curr_mat = all_mats{m};
    curr_mat(:,16,21:24) = nan; curr_mat(16,:,21:24) = nan;     % Removing the data of four subjects without stimulus 16
    % Testing the coding of distances in different models (full models, or only within and between museums)
    for i = 1:size(curr_mat,3)
        %     curr_mean_within = (curr_mat(1:8,1:8,i) + curr_mat(9:16,9:16,i));
        %     curr_mean_within = (curr_mat(1:8,1:8,i) + curr_mat(9:16,9:16,i) + curr_mat(1:8,9:16,i) + curr_mat(9:16,1:8,i));
        %     curr_mean_within = (curr_mat(1:8,1:8,i) + curr_mat(9:16,9:16,i) + curr_mat(1:8,9:16,i) + curr_mat(9:16,1:8,i));
        %     curr_mean_within = (curr_mat(1:8,1:8,i) + curr_mat(9:16,9:16,i) + curr_mat(1:8,9:16,i) + curr_mat(9:16,1:8,i));
        %     curr_mean_within = (curr_mat(1:8,1:8,i) + curr_mat(9:16,9:16,i) + curr_mat(1:8,9:16,i) + curr_mat(9:16,1:8,i));
        curr_subject_mean_mat = curr_mat(:,:,i) + curr_mat(:,:,i)';         % Averaging the lower and upper triangle - creating a symmetric matrix
        
        % Comparing each model to the data, using the lower triangle only
        dist_corr_within(i,m) = corr(mat_distances_within(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
        dist_corr_overlay(i,m) = corr(mat_distances_overlay(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
        dist_corr_overlay_between(i,m) = corr(mat_distances_overlay_between(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
        dist_corr_full(i,m) = corr(mat_distances_full_model(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
        dist_corr_full_between(i,m) = corr(mat_distances_full_model_between(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
        heading_corr_global(i,m) = corr(mat_global_heading(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
        heading_corr_local(i,m) = corr(mat_local_heading(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))),'Type', 'Spearman','rows','complete');
    end
end
% Testing the full models
p1=[];p2=[];p6=[];p7=[];
STATS1={};STATS2={};STATS6={};STATS7={};
for i=1:size(dist_corr_full,2)
    [p1(i),~,STATS1{i}]=signrank(dist_corr_full(:,i),[],'tail','right');            % Integration model - location proximity
    [p2(i),~,STATS2{i}]=signrank(dist_corr_overlay(:,i),[],'tail','right');         % Schematization model - location proximity
    [p6(i),~,STATS6{i}]=signrank(heading_corr_local(:,i),[],'tail','right');        % Schematization model - heading similarity
    [p7(i),~,STATS7{i}]=signrank(heading_corr_global(:,i),[],'tail','right');       % Integration model - heading similarity
end
% [~,p1,~,STATS1]=ttest(dist_corr_full,[],'tail','right'); disp(['integration_model ' num2str(p1)])
% [~,p2,~,STATS2]=ttest(dist_corr_overlay,[],'tail','right'); disp(['overlay_model ' num2str(p2)])
% [~,p6,~,STATS6]=ttest(heading_corr_local,[],'tail','right'); disp(['heading_corr_local ' num2str(p6)])
% [~,p7,~,STATS7]=ttest(heading_corr_global,[],'tail','right'); disp(['heading_corr_global ' num2str(p7)])
% % % % Testing within and between museum distances under each model
% % [~,p3,~,STATS3]=ttest(dist_corr_within,[],'tail','right'); disp(['within_segment ' num2str(p3)])
% % [~,p4,~,STATS4]=ttest(dist_corr_full_between,[],'tail','right'); disp(['full_model_between_segments ' num2str(p4)])
% % [~,p5,~,STATS5]=ttest(dist_corr_overlay_between,[],'tail','right'); disp(['overlay_between_segments ' num2str(p5)])




%% Compare models using BIC
% r2_model1_distances = nan(size(curr_mat,3),length(all_mats)); r2_model2_overlay = nan(size(curr_mat,3),length(all_mats));
% aic_model1_distances = nan(size(curr_mat,3),length(all_mats)); aic_model2_overlay = nan(size(curr_mat,3),length(all_mats));
% bic_model1_distances = nan(size(curr_mat,3),length(all_mats)); bic_model2_overlay = nan(size(curr_mat,3),length(all_mats));
r2_model1_distances = nan(1, length(all_mats)); r2_model2_overlay = nan(1, length(all_mats));
aic_model1_distances = nan(1, length(all_mats)); aic_model2_overlay = nan(1, length(all_mats));
bic_model1_distances = nan(1, length(all_mats)); bic_model2_overlay = nan(1, length(all_mats));
for m = 1:length(all_mats)      % Going over ROIs
    curr_mat = all_mats{m};
    curr_mat(:,16,21:24) = nan; curr_mat(16,:,21:24) = nan;     % Removing the data of four subjects without stimulus 16
    % Testing the coding of distances in different models (full models, or only within and between museums)
    %     for i = 1:size(curr_mat,3)
    %         curr_subject_mean_mat = curr_mat(:,:,i) + curr_mat(:,:,i)';         % Averaging the lower and upper triangle - creating a symmetric matrix
    %
    %         % Comparing each model to the data, using the lower triangle only
    %         fit_model1 = fitlm(mat_distances_full_model(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))));
    %         fit_model2 = fitlm(mat_distances_overlay(find(tril(ones(16),-1))), curr_subject_mean_mat(find(tril(ones(16),-1))));
    %         r2_model1_distances(i, m) = fit_model1.Rsquared.Adjusted;
    %         r2_model2_overlay(i, m) = fit_model2.Rsquared.Adjusted;
    %         logL = fit_model1.LogLikelihood; [aic,bic] = aicbic(logL, 1, 120); aic_model1_distances(i, m) = aic; bic_model1_distances(i, m) = bic;
    %         logL = fit_model2.LogLikelihood; [aic,bic] = aicbic(logL, 1, 120); aic_model2_overlay(i, m) = aic; bic_model2_overlay(i, m) = bic;
    %     end
    
    curr_mat = nanmean(curr_mat, 3);
    curr_mat = curr_mat + curr_mat';
    
    % Comparing each model to the data, using the lower triangle only
    fit_model1 = fitlm(mat_distances_full_model(find(tril(ones(16),-1))), curr_mat(find(tril(ones(16),-1))));
    fit_model2 = fitlm(mat_distances_overlay(find(tril(ones(16),-1))), curr_mat(find(tril(ones(16),-1))));
    r2_model1_distances(i, m) = fit_model1.Rsquared.Adjusted;
    r2_model2_overlay(i, m) = fit_model2.Rsquared.Adjusted;
    logL = fit_model1.LogLikelihood; [aic,bic] = aicbic(logL, 1, 120); aic_model1_distances(m) = aic; bic_model1_distances(m) = bic;
    logL = fit_model2.LogLikelihood; [aic,bic] = aicbic(logL, 1, 120); aic_model2_overlay(m) = aic; bic_model2_overlay(m) = bic;
end





%% Testing a model of wall to the left vs. to the right
curr_mat = corrs_rsc;
curr_mat(:,16,21:24) = nan; curr_mat(16,:,21:24) = nan; % Removing four subjects without stimulus 16
dist_corr_wall = [];
mat_walls_direction = zeros(8);
for i=1:2:15
    for j=1:2:15
        mat_walls_direction(i,j) = 1;
    end
end
for i=2:2:16
    for j=2:2:16
        mat_walls_direction(i,j) = 1;
    end
end
for i = 1:size(curr_mat,3)
    curr_mean_within = curr_mat(:,:,i) + curr_mat(:,:,i)';
    dist_corr_wall(i) = corr(mat_walls_direction(find(tril(ones(16),-1))), curr_mean_within(find(tril(ones(16),-1))),'Type', 'Spearman');
end
[~,p]=ttest(dist_corr_wall,[],'tail','right')



%% MDS within and between regions
curr_mat = corrs_rsc;
curr_mat(:,16,21:24) = nan; curr_mat(16,:,21:24) = nan; % Removing four subjects without stimulus 16
nums={}; for i=1:8, nums{i}=num2str(i);end
mean_within = nanmean(curr_mat,3);
mean_within = (mean_within(1:8,1:8) + mean_within(9:16,9:16))/2;
mean_within = (mean_within + mean_within') / 2;
mean_within = mean_within - min(mean_within(:));
mean_within = 1 - mean_within;
mean_within(find(eye(8))) = 0;
Y = cmdscale(mean_within,2);
figure;plot(Y(:,1),Y(:,2),'*'), text(Y(:,1),Y(:,2),nums), title('within region')
mean_between = nanmean(curr_mat,3);
mean_between = (mean_between(1:8,9:16) + mean_between(9:16,1:8))/2;
mean_between = (mean_between + mean_between') / 2;
mean_between = mean_between - min(mean_between(:));
mean_between = 1 - mean_between;
mean_between(find(eye(8))) = 0;
Y = cmdscale(mean_between,2);
figure;plot(Y(:,1),Y(:,2),'*'), text(Y(:,1),Y(:,2),nums), title('between regions')
mean_all_overlay = nanmean(curr_mat,3);
mean_all_overlay = (mean_all_overlay(1:8,9:16) + mean_all_overlay(9:16,1:8) + mean_all_overlay(1:8,1:8) + mean_all_overlay(9:16,9:16))/2;
mean_all_overlay = (mean_all_overlay + mean_all_overlay') / 2;
mean_all_overlay = mean_all_overlay - min(mean_all_overlay(:));
mean_all_overlay = 1 - mean_all_overlay;
mean_all_overlay(find(eye(8))) = 0;
Y = cmdscale(mean_all_overlay,2);
figure;plot(Y(:,1),Y(:,2),'*'), text(Y(:,1),Y(:,2),nums), title('overlay full model')


