function compareGrps_avgRBCV_diam(fname) 
%
%written by Kira, June 2018
%function to compare average RBCV (overtime) and avg linescan diam (from 
%static image) from 2+ groups, will also save correlation coefficient
%must have previously run findAvg_RBCV_diam.m for each group separaretly
%
%INPUTS - 
%fname = this is a top directory which contains 2+ groups, e.g. V1
%and HC, or APOE3 and APOE4, each of these group folders must contain a
%'RBCV_meanstd.mat' file (from findAvg_RBCV_diam.m script) 
%
%OUTPUTS - 
%will save a mat file with all the data collated into one cell (this
%includes the results from the correlation coefficient), and will also save
%two figures - the scatterplots of data for each group, and barcharts to
%show averages 
%
%Requires the following additional functions to run:
%findFolders

%find all the folders which have the RCBV_meanstd.mat file 
find_mat_file=findFolders(fname,'RBCV_meanstd.mat');

%colour string for plots
colorstring={'b','m','k','c','y','g','r'}; 

%% plot comparative scatter plots, and run pearson's
%open figure to plot data
figure; 
%find size of screen (for figure sizing)
screenSz=get(0,'Screensize');
%make fig size of screen
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3)/2 screenSz(4)]);
hold on; 
%loop through the mat files containing averaged vel and diam data
for a = 1:size(find_mat_file,2) %loop to load data
    
    %find the individual exp dirs that contain mat files
    [expDir,matfile]=fileparts(find_mat_file{1,a}); 
    %load the data
    load([expDir,filesep,matfile]);
    
    %put data into cell for each group:
    %mean diam from static image for each session
    data{a}.meanSz_all=mean_sz;
    data{a}.meanSz_mean=nanmean(mean_sz); 
    data{a}.meanSz_std=nanstd(mean_sz); 
    [data{a}.meanSz_sem]=getSEM(mean_sz,2); 
    %mean RBCV across time for each session
    data{a}.meanVel_all=mean_vel;
    data{a}.meanVel_mean=nanmean(mean_vel); 
    data{a}.meanVel_std=nanstd(mean_vel); 
    [data{a}.meanVel_sem]=getSEM(mean_vel,2); 
    %std of RBCV across time for each session
    data{a}.stdVel_all=std_vel;
    data{a}.stdVel_mean=nanmean(std_vel); 
    data{a}.stdVel_std=nanstd(std_vel); 
    [data{a}.stdVel_sem]=getSEM(std_vel,2); 
    %include grp label so know which grp 
    data{a}.grpLabel=grpNm; 
    
    %% plot
    ax(a) = subplot(size(find_mat_file,2),round(size(find_mat_file,2)/2),a);
    scatter(ax(a),data{a}.meanSz_all,data{a}.meanVel_all,  ...
        'MarkerEdgeColor', colorstring{a}, 'MarkerFaceColor', ...
        colorstring{a}, 'Marker','*');
    %Add least-squares line to scatter plot
    data{a}.lsline=lsline(ax(a));
    title([data{a}.grpLabel]);
    xlabel('vessel diameter, um'); xlim([3 9]);
    ylabel('RBCV, mm/s'); ylim([0 10]); 
    
    %Calculate the correlation between X and Y using corrcoef for each grp 
    [data{a}.rho,data{a}.pval] = corrcoef(data{a}.meanSz_all, ...
        data{a}.meanVel_all);
      
end %end of loop to load data
%save figure
figsave = fullfile(fname,['RBCV_diam_scatter_compareGrps.png']);
saveas(gcf, figsave);
close;

%save the data
matFile = fullfile(fname, ['RBCV_diam_AllGrps.mat']);
save(matFile, 'data'); 

%% look at averages for both regions

for a = 1:size(data,2) %loop to extract data to plot
    %extract data for plots
    diam_mean(a)=data{a}.meanSz_mean;
    diam_sem(a)=data{a}.meanSz_sem;
    vel_mean(a)=data{a}.meanVel_mean;
    vel_sem(a)=data{a}.meanVel_sem;
    vel_std(a)=data{a}.stdVel_mean;
    vel_std_sem(a)=data{a}.stdVel_sem;
end %end of extracting data to plot

%plot average RBCV and average diam for two groups 
figure;
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3)/2 screenSz(4)]);
subplot(311);
bar(diam_mean,'m');
hold on;
errorbar(diam_mean, diam_sem, 'm.');
title('Average Diameter'); 
ylabel('um'); 
set(gca, 'XTickLabel',[{data{1}.grpLabel}, {data{2}.grpLabel}]);
subplot(312)
bar(vel_mean,'r');
hold on;
errorbar(vel_mean, vel_sem, 'r.');
title('Average Mean of RBCV');
ylabel('mms'); 
set(gca, 'XTickLabel',[{data{1}.grpLabel}, {data{2}.grpLabel}]);
subplot(313)
bar(vel_std,'r');
hold on;
errorbar(vel_std, vel_std_sem, 'r.');
title('Average Std of RBCV');
ylabel('mms'); 
set(gca, 'XTickLabel',[{data{1}.grpLabel}, {data{2}.grpLabel}]);
%save figure
figsave = fullfile(fname,['RBCV_diam_barchart_compareGrps.png']);
saveas(gcf, figsave);
close;


end %end of function 
