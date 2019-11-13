

%NB need to run following codes before this will work:
% - extractProbe - to extract haem data
% - cutByCMRO2Events_probe to cut around individual CMRO2 events

%this script will also call the following functions:
% - findFolders
% - sortDataByResponsiveTrials
% - tsParameters

clear all; close all; 
%'D:\Dropbox (Brain Energy Lab)\Everything\Kira\Projects\E3vsE4_HC\dataSortedByAnalysis\Oxyprobe_CMRO2peaks'
fname=cd; 

% important variables to send in for classifying the vessel as responsive
% or not:
%time to check diameter curve (to see if resp or not) - NB/ calc peak is
%aligned to 5s
preferences.respPeriod_diam = [5 10]; %seconds
%set how many standard deviations above thresh to deem as responsive
%aka the std multiplier
preferences.std_multiplier = 1; %SD %calc is 2, but the diameter peaks less high
%time in seconds to average either side of the max peak, i.e. to see if
%the peak is above the std dev of the baseline - this may need to be
%larger for slower responses (e.g. haem vs calc)
preferences.avgPeak = 0.5; %seconds

%find size of screen (for figure sizing)
screenSz=get(0,'Screensize');

%% find all the subfolders - i.e. groups - in the top dir
%only take the names which are less than 4 letters long, to exclude non
%label folders, so make sure your group names are short enough
listing=dir(fname);
for a = 1:size(listing,1)-2
    if size(listing(2+a).name,2) <= 4
        grpNames{a}=listing(2+a).name;
    else
        grpNames{a}=[];
    end
end
%remove the empty cells 
grpNames = grpNames(~cellfun('isempty',grpNames));

%% loop the groups, extract the data
for n = 1:size(grpNames,2)
    
    clear findmatfiles; 
    %inside each group subfolder, find the mat file
    findmatfiles=findFolders([fname,filesep,grpNames{1,n}], ...
        'dataCutByCMRO2peaks.mat');
    %inform user of progress
    disp(['loading data for ', grpNames{1,n}, ' grp...']); 
    
    clear a; 
    for a = 1:size(findmatfiles,2) %loop through expDirs with data
        
        %find individual exp dir
        [expDir,matfile]=fileparts(findmatfiles{1,a});
        %load data
        load(fullfile([expDir,filesep,matfile]));

        %put data into a cell - keep them separate so you can check for any
        %outliers and know where they came from
        %resize time dim so all are same size
        %calcium
        calcAll_sep{n}.cmro2{a} = imresize(calcSort.calc, ...
            [size(calcSort.calc,1) 65]);
        %loco
        calcAll_sep{n}.loco_all{a} = imresize(calcSort.loco_all, ...
            [size(calcSort.loco_all,1) 65]);
        %stim
        calcAll_sep{n}.stim_all{a} = imresize(calcSort.stim_all, ...
            [size(calcSort.stim_all,1) 65]);
        %haem real
        for i = 1:size(calcSort.haem,2) %loop all haem param
            calcAll_sep{n}.haem{a}(:,i,:) = imresize(squeeze( ...
                calcSort.haem(:,i,:)), ...
                [size(calcSort.haem,1) 65]);
        end
        %haem bootstrap
        for i = 1:size(calcSort.haemBS,2) %loop haem param
            for j = 1:size(calcSort.haemBS,3) %loop iterations
                calcAll_sep{n}.haemBS{a}(:,i,j,:) = imresize(squeeze( ...
                    calcSort.haemBS(:,i,j,:)), ...
                    [size(calcSort.haemBS,1)  65]);
            end
        end
        
        %time
        calcAll_sep{n}.time{a} = imresize(calcSort.time, [1 65]);
        %fps
        calcAll_sep{n}.fps{a} = fps;
        %grp names
        calcAll_sep{n}.grpLabel = grpNames{n}; 
        
        %binary for loco / stim traces
        calcAll_sep{n}.stim{a} = calcSort.stim; 
        calcAll_sep{n}.loco{a} = calcSort.loco; 
        
        clear counter;
        %put all trials into one cont trace
        if a == 1 %if 1st trial create var
            calcAll{n}.counter = size(calcSort.haem,1);
            %cmro2 traces
            calcAll{n}.cmro2 = imresize(calcSort.calc, ...
                [size(calcSort.calc,1) 65]);
            %loco traces
            calcAll{n}.loco_all = imresize(calcSort.loco_all, ...
                [size(calcSort.loco_all,1) 65]);
            %stim traces
            calcAll{n}.stim_all = imresize(calcSort.stim_all, ...
                [size(calcSort.stim_all,1) 65]);
            %haem traces
            for i = 1:size(calcSort.haem,2)
                calcAll{n}.haem(:,i,:) = imresize(squeeze(calcSort.haem(:,i,:)), ...
                    [size(calcSort.haem,1) 65]);
            end
            %shuffled haem traces
            for i = 1:size(calcSort.haemBS,2)
                for j = 1:size(calcSort.haemBS,3)
                    calcAll{n}.haemBS(:,i,j,:) = imresize(squeeze( ...
                        calcSort.haemBS(:,i,j,:)), ...
                        [size(calcSort.haemBS,1) 65]);
                end
            end
            
            calcAll{n}.time = calcAll_sep{1,n}.time{1,1};
            calcAll{n}.fps = repmat(fps, [1 size(calcSort.calc,1)]);
            calcAll{n}.stim = calcAll_sep{n}.stim{a};
            calcAll{n}.loco = calcAll_sep{n}.loco{a};
            
        else %post-first trial, so merge with rest of data
            
            %resize all the calcium events so they are all same across
            %events
            calcAll{n}.counter = [calcAll{n}.counter, size(calcSort.haem,1)];
            %CMRO2
            calcAll{n}.cmro2 = [calcAll{n}.cmro2; imresize(calcSort.calc, ...
                [size(calcSort.calc,1) 65])];
            %loco
            calcAll{n}.loco_all = [calcAll{n}.loco_all; imresize(calcSort.loco_all, ...
                [size(calcSort.loco_all,1) 65])];
            %stim
            calcAll{n}.stim_all = [calcAll{n}.stim_all; imresize(calcSort.stim_all, ...
                [size(calcSort.stim_all,1) 65])];
            %haem
            for i = 1:size(calcSort.haem,2)
                ttt(:,i,:) = imresize(squeeze(calcSort.haem(:,i,:)), ...
                    [size(calcSort.haem,1) 65]);
            end
            calcAll{n}.haem = [calcAll{n}.haem; ttt];
            clear ttt; 
            %haem BS
            for i = 1:size(calcSort.haemBS,2)
                for j = 1:size(calcSort.haemBS,3)
                    ttt(:,i,j,:) = imresize(squeeze(calcSort.haemBS ...
                        (:,i,j,:)),[size(calcSort.haemBS,1) 65]);
                end
            end
            calcAll{n}.haemBS = [calcAll{n}.haemBS; ttt];
            clear ttt;
            calcAll{n}.fps = [calcAll{n}.fps, ...
                repmat(fps ,[1 size(calcSort.calc,1)])];
            calcAll{n}.stim = [calcAll{n}.stim, calcAll_sep{n}.stim{a}];
            calcAll{n}.loco = [calcAll{n}.loco, calcAll_sep{n}.loco{a}];
        end
        calcAll{n}.grpLabel = grpNames{n}; 
        
        
    end %end of looping through expDirs containing needed data
    
    clearvars -except fname calcAll_sep calcAll grpNames preferences screenSz;
  
end %end of looping groups


%% ALL the data - i.e. not sorted by responsive Hbt:

%find time series parameters (e.g. max peak) for ALL the data across groups
for n = 1:size(grpNames,2) %loop grps
    for a = 1:size(calcAll{1,n}.cmro2,1) %loop cmro2 events
        
        %call function to assess the ts parameters for the calc data
        %take onset of calc until 2s post onset to assess data
        clear ttt; 
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.cmro2(a,:),...
            calcAll{1,n}.time,round(size(calcAll{1,n}.time,2)/ ...
            max(calcAll{1,n}.time)),[find(calcAll{1,n}.time>=3,1) ...
            find(calcAll{1,n}.time>=8,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_cmro2(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_cmro2(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_cmro2(a)=ttt.t2p;

        %call function to assess ts param for corresponding haem data
        %take onset of calc until 5s post onset (as haem is slow)
        for b = 1:size(calcAll{1,n}.haem,2) %loop haem types
            clear ttt;
            [ttt] = findtsParameters_calcPeak(squeeze(calcAll{1,n}.haem(a,b,:))',...
                calcAll{1,n}.time,round(size(calcAll{1,n}.time,2)/ ...
                max(calcAll{1,n}.time)),[find(calcAll{1,n}.time>=5,1) ...
                find(calcAll{1,n}.time>=8,1)]);
            %save out data of interest into struct
            calcAll{1,n}.AUC_haem(b,a)=ttt.AUC;
            calcAll{1,n}.maxPeak_haem(b,a)=ttt.maxPeak;
            calcAll{1,n}.t2p_haem(b,a)=ttt.t2p;
        end %end of loop haem types
        
    end %end of loop cmro2 events
end %end of loop grps

%get sem for the time series and the ts parameters for ALL data
for n = 1:size(grpNames,2)
    
    %cmro2 time series and param:
    [calcAll{1,n}.cmro2_SEM]=getSEM(calcAll{1,n}.cmro2,1);
    [calcAll{1,n}.cmro2Peak_SEM]=getSEM(calcAll{1,n}.maxPeak_cmro2,2);
    [calcAll{1,n}.cmro2AUC_SEM]=getSEM(calcAll{1,n}.AUC_cmro2,2);
    [calcAll{1,n}.cmro2t2p_SEM]=getSEM(calcAll{1,n}.t2p_cmro2,2);
    
    %haem time series and param
    %get SEM for each haem trace across detected cmro2 events
    for a = 1:size(calcAll{1,n}.haem,2)
        [calcAll{1,n}.haem_SEM(a,:)]=getSEM(squeeze(calcAll{1,n}.haem(:,a,:)),1);
    end
    %get SEM for each of the ts param across cmro2 events, per haem type
    [calcAll{1,n}.haemPeak_SEM]=getSEM(calcAll{1,n}.maxPeak_haem,2);
    [calcAll{1,n}.haemAUC_SEM]=getSEM(calcAll{1,n}.AUC_haem,2);
    [calcAll{1,n}.haemt2p_SEM]=getSEM(calcAll{1,n}.t2p_haem,2);
    
end

%stats on cmro2 peaks and hbt peaks between grps for ALL data:
%cmro2
[~,p] = ttest2(calcAll{1,1}.maxPeak_cmro2, ...
    calcAll{1,2}.maxPeak_cmro2);
if p < 0.05
    disp(['t-test, ALL data, cmro2 peaks: p=',num2str(p),', sig']);
else
    disp(['t-test, ALL data, cmro2 peaks: p=',num2str(p),', not sig']);
end
%flux
[~,p1] = ttest2(calcAll{1,1}.maxPeak_haem(1,:), ...
    calcAll{1,2}.maxPeak_haem(1,:));
if p1 < 0.05
    disp(['t-test, ALL data, flux peaks: p=',num2str(p1),', sig']);
else
    disp(['t-test, ALL data, flux peaks: p=',num2str(p1),', not sig']);
end
%so2
[~,p2] = ttest2(calcAll{1,1}.maxPeak_haem(3,:), ...
    calcAll{1,2}.maxPeak_haem(3,:));
if p2 < 0.05
    disp(['t-test, ALL data, so2 peaks: p=',num2str(p2),', sig']);
else
    disp(['t-test, ALL data, so2 peaks: p=',num2str(p2),', not sig']);
end
%hbt
[~,p3] = ttest2(calcAll{1,1}.maxPeak_haem(6,:), ...
    calcAll{1,2}.maxPeak_haem(6,:));
if p3 < 0.05
    disp(['t-test, ALL data, hbt peaks: p=',num2str(p3),', sig']);
else
    disp(['t-test, ALL data, hbt peaks: p=',num2str(p3),', not sig']);
end

%% plot traces for ALL data, i.e. not sorted by Hbt responsiveness yet 
%CMRO2 and Hbt
figure;
%make fig size of screen
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);
           
%%%%%%%%%%%%% CMRO2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(4,3,1);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.cmro2,1), 'g', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.cmro2,1), ...
    calcAll{1}.cmro2_SEM, 'lineProps','g');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.01 0.2],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' CMRO2']);

%main title to show it is ALL data, i.e. not sorted by resp
mtit(gcf,'ALL data','fontsize', 14,'color', [1 0 0],'xoff',1.25,'yoff',0.25);

%time series trace grp 2
ax2 = subplot(4,3,2);
plot(calcAll{1}.time, ...
    nanmean(calcAll{2}.cmro2,1), 'g', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.cmro2,1), ...
    calcAll{2}.cmro2_SEM, 'lineProps','g');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.01 0.2],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2}, ' CMRO2']); 

%link cmro2 time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for cmro2 traces
subplot(4,3,3);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_cmro2), ...
    nanmean(calcAll{1,2}.maxPeak_cmro2)],'g');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_cmro2), ...
    nanmean(calcAll{1,2}.maxPeak_cmro2)], ...
    [calcAll{1,1}.cmro2Peak_SEM, ...
    calcAll{1,2}.cmro2Peak_SEM],'g.');
xlim([0 3]);
title(['CMRO2 Peak, p=', num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}, ', n=',  ...
    num2str(size(calcAll{1,1}.cmro2,1))], ...
    [grpNames{2}, ', n=', ...
    num2str(size(calcAll{1,2}.cmro2,1))]});

%%%%%%%%%%%%% Flux: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(4,3,4);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem(:,1,:),1)), 'm', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem(:,1,:),1)), ...
    calcAll{1}.haem_SEM(1,:), 'lineProps','m');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.25],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' Flux']);

%time series trace grp 2
ax2 = subplot(4,3,5);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem(:,1,:),1)), 'm', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem(:,1,:),1)), ...
    calcAll{2}.haem_SEM(1,:), 'lineProps','m');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.25],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' Flux']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for flux traces
subplot(4,3,6);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem(1,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem(1,:))],'m');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem(1,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem(1,:))], ...
    [calcAll{1,1}.haemPeak_SEM(1), ...
    calcAll{1,2}.haemPeak_SEM(1)],'m.');
xlim([0 3]);
title(['Flux Peak, p=', num2str(p1)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%% So2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(4,3,7);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem(:,3,:),1)), 'k', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem(:,3,:),1)), ...
    calcAll{1}.haem_SEM(3,:), 'lineProps','k');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.05],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' So2']);

%time series trace grp 2
ax2 = subplot(4,3,8);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem(:,3,:),1)), 'k', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem(:,3,:),1)), ...
    calcAll{2}.haem_SEM(3,:), 'lineProps','k');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.05],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' So2']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for speed traces
subplot(4,3,9);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem(3,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem(3,:))],'k');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem(3,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem(3,:))], ...
    [calcAll{1,1}.haemPeak_SEM(3), ...
    calcAll{1,2}.haemPeak_SEM(3)],'k.');
xlim([0 3]);
title(['So2 Peak, p=',num2str(p2)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%% Hbt: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(4,3,10);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem(:,6,:),1)), 'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem(:,6,:),1)), ...
    calcAll{1}.haem_SEM(6,:), 'lineProps','r');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.008 0.008],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' Hbt']);

%time series trace grp 2
ax2 = subplot(4,3,11);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem(:,6,:),1)), 'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem(:,6,:),1)), ...
    calcAll{2}.haem_SEM(6,:), 'lineProps','r');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.008 0.008],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' Hbt']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for speed traces
subplot(4,3,12);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem(6,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem(6,:))],'r');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem(6,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem(6,:))], ...
    [calcAll{1,1}.haemPeak_SEM(6), ...
    calcAll{1,2}.haemPeak_SEM(6)],'r.');
xlim([0 3]);
title(['Hbt Peak, p=', num2str(p3)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%%%%%%%
%save plot into exp_dir as png
figsave = fullfile(fname,['CMRO2traces_HaemTraces_PeakBarchart_ALLdata.png']);
saveas(gcf, figsave);
close; 

clear p p1 p2 p3; 

%% calculate the NVCindex for ALL data

for n = 1:size(calcAll,2) %loop grps
    for b = 1:size(calcAll{n}.cmro2,1) %loop cmro2 events
        %calc NVCind
        %divide vessel peak by neuron peak 
        calcAll{1,n}.NVCind(b) = calcAll{1,n}.maxPeak_haem(6,b) ...
            / calcAll{1,n}.maxPeak_cmro2(b);
    end %end of loop cmro2 events
    %get SEM for NVCindex
    [calcAll{1,n}.NVCind_SEM] = getSEM(calcAll{1,n}.NVCind,2);
end %end of loop grps

%stats on NVCind between grps
[~,p] = ttest2(calcAll{1,1}.NVCind, calcAll{1,2}.NVCind);
if p < 0.05
    disp(['t-test, ALL data, NVCind: p=',num2str(p),', sig']);
else
    disp(['t-test, ALL data, NVCind: p=',num2str(p),', not sig']);
end

%plot NVCindex
figure; 
bar(1:2,[nanmean(calcAll{1,1}.NVCind), nanmean(calcAll{1,2}.NVCind)],'y');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.NVCind), nanmean(calcAll{1,2}.NVCind)], ...
    [calcAll{1,1}.NVCind_SEM, calcAll{1,2}.NVCind_SEM],'y.');
xlim([0 3]);
title(['NVCindex, RESP data, p=',num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});
%save plot into exp_dir as png
figsave = fullfile(fname,['NVCindex_allData.png']);
saveas(gcf, figsave);
close; 

clear p; 

%% classify the HBT responses as responsive or not to CMRO2:

%get most responsive index for haem REAL and haem BS
for b = 1:size(calcAll,2) %loop the grps
    
    disp(['calculating response rates for real and bootstrap data, grp',...
        num2str(b)]); 
    
    calcAll{b}.responsiveIndBS = nan(size(calcAll{1,b}.haemBS,1), ...
        size(calcAll{1,b}.haemBS,3));
    for c = 1:size(calcAll{b}.haem,1) %loop detected cmro2 events
        %send in the Hbt trace to decide if 'vessels' are resp or not to cmro2
        [calcAll{b}.responsiveInd(c)] = sortDataByResponsiveTrials ...
            (squeeze(calcAll{b}.haem(c,6,:))',calcAll{b}.time, ...
            preferences, 'diam');
        for d = 1:100 %loop the 100 iterations
            [calcAll{b}.responsiveIndBS(c,d)] = sortDataByResponsiveTrials ...
            (squeeze(calcAll{b}.haemBS(c,6,d,:))',calcAll{b}.time, ...
            preferences, 'diam');
        end
    end
    
    %check for any NaNs and make them zeros -- thus rendering them as
    %'unresponsive' trials - NB/ will affect the % calc slightly as more
    %deemed unresponsive 
    calcAll{b}.responsiveInd(find(isnan(calcAll{b}.responsiveInd)))=0; 
    calcAll{b}.responsiveIndBS(find(isnan(calcAll{b}.responsiveIndBS)))=0; 
    
    %find overall % of responsive haem to cmro2 events
    calcAll{b}.respPcnt = (sum(calcAll{b}.responsiveInd)/ ...
        size(calcAll{b}.haem,1))*100; 
    
    %calculate the % responsive across iterations for each calc event
    calcAll{b}.respPcntBS_all = (sum(calcAll{b}.responsiveIndBS')/ ...
        size(calcAll{b}.haem,1)) *100;
    
    %mean across iterations per calc event, to get the overall mean respond
    calcAll{b}.respPcntBS_avg = nanmean(calcAll{b}.respPcntBS_all);
      
end %end of loop grps for categorising Hbt as resp to cmro2 events or not

% %%% sanity check - check all responses not coming from same animal
% calcAll{1}.expLabel(find(calcAll{1}.responsiveInd==1))
% calcAll{2}.expLabel(find(calcAll{2}.responsiveInd==1))

% plot responsiveness rates in bar chart

%plot
y = [calcAll{1}.respPcnt, ...
    100 - calcAll{1}.respPcnt; ...
    
    calcAll{1}.respPcntBS_avg, ...
    100 - calcAll{1}.respPcntBS_avg; ...
    
    calcAll{2}.respPcnt, ...
    100 - calcAll{2}.respPcnt; ...
    
    calcAll{2}.respPcntBS_avg, ...
    100 - calcAll{2}.respPcntBS_avg]; 

datalabels = categorical({[grpNames{1},' real'],[grpNames{1}, ' BS'],...
    [grpNames{2},' real'],[grpNames{2},' BS']});
datalabels = reordercats(datalabels,{[grpNames{1},' real'],[grpNames{1}, ' BS'],...
    [grpNames{2},' real'],[grpNames{2},' BS']});

figure;
bar(datalabels, y, 'stacked'); 
legend('R','NR'); 
ylim([0 110]); ylabel('%'); 
title('Resp Rate'); 
%save plot
figsave = fullfile(fname,['responsivenessRates_Hbt_compGrps.png']);
saveas(gcf, figsave);
close; 

%% extract those traces which are haem responsive only

disp('Finding CMRO2 events which lead to responsive Hbt only...'); 

for b = 1:size(calcAll,2) %loop grps
    %extract time series data:
    %cmro2
    calcAll{b}.cmro2_resp = calcAll{b}.cmro2(find(calcAll{b}.responsiveInd),:);
    %haem
    calcAll{b}.haem_resp = calcAll{b}.haem(find(calcAll{b}.responsiveInd),:,:);
    %get errorbars:
    [calcAll{b}.cmro2_resp_SEM] = getSEM(calcAll{b}.cmro2_resp,1);
    for c = 1:size(calcAll{b}.haem_resp,2) %loop haem param
        [calcAll{b}.haem_resp_SEM(c,:)] = getSEM(squeeze( ...
            calcAll{b}.haem_resp(:,c,:)),1);
    end %end of loop haem param
end %end of loop grps

%ts parameters on responsive traces only: 
for n = 1:size(grpNames,2) %loop grps
    for a = 1:size(calcAll{1,n}.cmro2_resp,1) %loop resp events
        
        %call function to assess the ts parameters for the calc data
        %take onset of calc until 2s post onset to assess data
        clear ttt; 
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.cmro2_resp(a,:),...
            calcAll{1,n}.time,round(size(calcAll{1,n}.time,2)/ ...
            max(calcAll{1,n}.time)),[find(calcAll{1,n}.time>=3,1) ...
            find(calcAll{1,n}.time>=8,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_cmro2_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_cmro2_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_cmro2_resp(a)=ttt.t2p;

        %call function to assess ts param for corresponding haem data
        %take onset of calc until 5s post onset (as haem is slow)
        for b = 1:size(calcAll{1,n}.haem_resp,2) %loop haem param
            clear ttt;
            [ttt] = findtsParameters_calcPeak(squeeze(calcAll{1,n}.haem_resp(a,b,:))',...
                calcAll{1,n}.time,round(size(calcAll{1,n}.time,2)/ ...
                max(calcAll{1,n}.time)),[find(calcAll{1,n}.time>=5,1) ...
                find(calcAll{1,n}.time>=8,1)]);
            %save out data of interest into struct
            calcAll{1,n}.AUC_haem_resp(b,a)=ttt.AUC;
            calcAll{1,n}.maxPeak_haem_resp(b,a)=ttt.maxPeak;
            calcAll{1,n}.t2p_haem_resp(b,a)=ttt.t2p;
        end %end of loop haem param

    end %end of loop resp events
    
    %get errorbars
    %for each of the ts param across cmro2 events, per haem type
    [calcAll{1,n}.maxPeak_cmro2_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_cmro2_resp,2);
    [calcAll{1,n}.maxPeak_haem_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_haem_resp,2);
    [calcAll{1,n}.AUC_haem_resp_SEM]=getSEM(calcAll{1,n}.AUC_haem_resp,2);
    [calcAll{1,n}.t2p_haem_resp_SEM]=getSEM(calcAll{1,n}.t2p_haem_resp,2);
    
end

% plot responsive traces only: 

%stats on cmro2 peaks and hbt peaks between grps:
%cmro2
[~,p] = ttest2(calcAll{1,1}.maxPeak_cmro2_resp, ...
    calcAll{1,2}.maxPeak_cmro2_resp);
if p < 0.05
    disp(['t-test, RESP data, cmro2 peaks: p=',num2str(p),', sig']);
else
    disp(['t-test, RESP data, cmro2 peaks: p=',num2str(p),', not sig']);
end
%flux
[~,p1] = ttest2(calcAll{1,1}.maxPeak_haem_resp(1,:), ...
    calcAll{1,2}.maxPeak_haem_resp(1,:));
if p1 < 0.05
    disp(['t-test, RESP data, flux peaks: p=',num2str(p1),', sig']);
else
    disp(['t-test, RESP data, flux peaks: p=',num2str(p1),', not sig']);
end
%so2
[~,p2] = ttest2(calcAll{1,1}.maxPeak_haem_resp(3,:), ...
    calcAll{1,2}.maxPeak_haem_resp(3,:));
if p2 < 0.05
    disp(['t-test, RESP data, so2 peaks: p=',num2str(p2),', sig']);
else
    disp(['t-test, RESP data, so2 peaks: p=',num2str(p2),', not sig']);
end
%hbt
[~,p3] = ttest2(calcAll{1,1}.maxPeak_haem_resp(6,:), ...
    calcAll{1,2}.maxPeak_haem_resp(6,:));
if p3 < 0.05
    disp(['t-test, RESP data, hbt peaks: p=',num2str(p3),', sig']);
else
    disp(['t-test, RESP data, hbt peaks: p=',num2str(p3),', not sig']);
end
%hbo
[~,p4] = ttest2(calcAll{1,1}.maxPeak_haem_resp(4,:), ...
    calcAll{1,2}.maxPeak_haem_resp(4,:));
if p4 < 0.05
    disp(['t-test, RESP data, hbo peaks: p=',num2str(p4),', sig']);
else
    disp(['t-test, RESP data, hbo peaks: p=',num2str(p4),', not sig']);
end
%hbr
[~,p5] = ttest2(calcAll{1,1}.maxPeak_haem_resp(5,:), ...
    calcAll{1,2}.maxPeak_haem_resp(5,:));
if p3 < 0.05
    disp(['t-test, RESP data, hbr peaks: p=',num2str(p5),', sig']);
else
    disp(['t-test, RESP data, hbr peaks: p=',num2str(p5),', not sig']);
end

figure;
%make fig size of screen
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);
           
%%%%%%%%%%%%% CMRO2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(6,3,1);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.cmro2_resp,1), 'g', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.cmro2_resp,1), ...
    calcAll{1}.cmro2_resp_SEM, 'lineProps','g');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.01 0.2],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' CMRO2']);

%main title to show it is ALL data, i.e. not sorted by resp
mtit(gcf,'RESP0NSIVE data','fontsize', 14,'color', [1 0 0],'xoff',1.25,'yoff',0.25);

%time series trace grp 2
ax2 = subplot(6,3,2);
plot(calcAll{1}.time, ...
    nanmean(calcAll{2}.cmro2_resp,1), 'g', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.cmro2_resp,1), ...
    calcAll{2}.cmro2_resp_SEM, 'lineProps','g');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.01 0.2],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2}, ' CMRO2']); 

%link cmro2 time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for cmro2 traces
subplot(6,3,3);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_cmro2_resp), ...
    nanmean(calcAll{1,2}.maxPeak_cmro2_resp)],'g');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_cmro2_resp), ...
    nanmean(calcAll{1,2}.maxPeak_cmro2_resp)], ...
    [calcAll{1,1}.maxPeak_cmro2_resp_SEM, ...
    calcAll{1,2}.maxPeak_cmro2_resp_SEM],'g.');
xlim([0 3]);
title(['CMRO2 Peak, p=', num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}, ', n=',  ...
    num2str(size(calcAll{1,1}.cmro2_resp,1))], ...
    [grpNames{2}, ', n=', ...
    num2str(size(calcAll{1,2}.cmro2_resp,1))]});

%%%%%%%%%%%%% Flux: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(6,3,4);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem_resp(:,1,:),1)), 'm', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem_resp(:,1,:),1)), ...
    calcAll{1}.haem_resp_SEM(1,:), 'lineProps','m');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.25],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' Flux']);

%time series trace grp 2
ax2 = subplot(6,3,5);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem_resp(:,1,:),1)), 'm', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem_resp(:,1,:),1)), ...
    calcAll{2}.haem_resp_SEM(1,:), 'lineProps','m');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.25],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' Flux']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for flux traces
subplot(6,3,6);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(1,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(1,:))],'m');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(1,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(1,:))], ...
    [calcAll{1,1}.maxPeak_haem_resp_SEM(1), ...
    calcAll{1,2}.maxPeak_haem_resp_SEM(1)],'m.');
xlim([0 3]);
title(['Flux Peak, p=',num2str(p1)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%% So2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(6,3,7);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem_resp(:,3,:),1)), 'k', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem_resp(:,3,:),1)), ...
    calcAll{1}.haem_resp_SEM(3,:), 'lineProps','k');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.05],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' So2']);

%time series trace grp 2
ax2 = subplot(6,3,8);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem_resp(:,3,:),1)), 'k', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem_resp(:,3,:),1)), ...
    calcAll{2}.haem_resp_SEM(3,:), 'lineProps','k');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.05 0.05],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' So2']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for speed traces
subplot(6,3,9);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(3,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(3,:))],'k');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(3,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(3,:))], ...
    [calcAll{1,1}.maxPeak_haem_resp_SEM(3), ...
    calcAll{1,2}.maxPeak_haem_resp_SEM(3)],'k.');
xlim([0 3]);
title(['So2 Peak, p=',num2str(p2)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%% Hbt: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(6,3,10);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem_resp(:,6,:),1)), 'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem_resp(:,6,:),1)), ...
    calcAll{1}.haem_resp_SEM(6,:), 'lineProps','r');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.008 0.008],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' Hbt']);

%time series trace grp 2
ax2 = subplot(6,3,11);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem_resp(:,6,:),1)), 'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem_resp(:,6,:),1)), ...
    calcAll{2}.haem_resp_SEM(6,:), 'lineProps','r');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.008 0.008],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' Hbt']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for speed traces
subplot(6,3,12);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(6,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(6,:))],'r');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(6,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(6,:))], ...
    [calcAll{1,1}.maxPeak_haem_resp_SEM(6), ...
    calcAll{1,2}.maxPeak_haem_resp_SEM(6)],'r.');
xlim([0 3]);
title(['Hbt Peak, p=', num2str(p3)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%% HbO: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(6,3,13);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem_resp(:,4,:),1)), 'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem_resp(:,4,:),1)), ...
    calcAll{1}.haem_resp_SEM(4,:), 'lineProps','r');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.008 0.008],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' HbO']);

%time series trace grp 2
ax2 = subplot(6,3,14);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem_resp(:,4,:),1)), 'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem_resp(:,4,:),1)), ...
    calcAll{2}.haem_resp_SEM(4,:), 'lineProps','r');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.008 0.008],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' HbO']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for HbO traces
subplot(6,3,15);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(4,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(4,:))],'r');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(4,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(4,:))], ...
    [calcAll{1,1}.maxPeak_haem_resp_SEM(4), ...
    calcAll{1,2}.maxPeak_haem_resp_SEM(4)],'r.');
xlim([0 3]);
title(['HbO Peak, p=', num2str(p4)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%% Hbr: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time series trace grp 1
ax1 = subplot(6,3,16);
plot(calcAll{1}.time, ...
    squeeze(nanmean(calcAll{1}.haem_resp(:,5,:),1)), 'b', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(calcAll{1}.haem_resp(:,5,:),1)), ...
    calcAll{1}.haem_resp_SEM(5,:), 'lineProps','b');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.001 0.05],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{1},' Hbr']);

%time series trace grp 2
ax2 = subplot(6,3,17);
plot(calcAll{2}.time, ...
    squeeze(nanmean(calcAll{2}.haem_resp(:,5,:),1)), 'b', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(calcAll{2}.haem_resp(:,5,:),1)), ...
    calcAll{2}.haem_resp_SEM(5,:), 'lineProps','b');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([5 5],[-0.001 0.05],'k--','LineWidth',2);
xlim([calcAll{1}.time(1) calcAll{1}.time(end)]); xlabel('Time(s)'); 
title([grpNames{2},' Hbr']);

%link haem time series plots
linkaxes([ax1,ax2],'xy'); 
clear ax1 ax2; 

%plot the max peak values for HbO traces
subplot(6,3,18);
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(5,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(5,:))],'b');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_haem_resp(5,:)), ...
    nanmean(calcAll{1,2}.maxPeak_haem_resp(5,:))], ...
    [calcAll{1,1}.maxPeak_haem_resp_SEM(5), ...
    calcAll{1,2}.maxPeak_haem_resp_SEM(5)],'b.');
xlim([0 3]);
title(['Hbr Peak, p=', num2str(p5)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});

%%%%%%%%%%%%%%%%%%
%save plot into exp_dir as png
figsave = fullfile(fname,['CMRO2traces_HaemTraces_PeakBarchart_RESPdata.png']);
saveas(gcf, figsave);
close; 

clear p p1 p2 p3 p4 p5; 

%% calculate the NVCindex for RESP data

for n = 1:size(calcAll,2) %loop grps
    for b = 1:size(calcAll{n}.cmro2_resp,1) %loop cmro2 events
        %calc NVCind
        %divide vessel peak by neuron peak 
        calcAll{1,n}.NVCind_resp(b) = calcAll{1,n}.maxPeak_haem_resp(6,b) ...
            / calcAll{1,n}.maxPeak_cmro2_resp(b);
    end %end of loop cmro2 events
    %get SEM for NVCindex
    [calcAll{1,n}.NVCind_resp_SEM] = getSEM(calcAll{1,n}.NVCind_resp,2);
end %end of loop grps

%stats on NVCind between grps
[~,p] = ttest2(calcAll{1,1}.NVCind_resp, calcAll{1,2}.NVCind_resp);
if p < 0.05
    disp(['t-test, RESP data, NVCind: p=',num2str(p),', sig']);
else
    disp(['t-test, RESP data, NVCind: p=',num2str(p),', not sig']);
end

%plot NVCindex
figure; 
bar(1:2,[nanmean(calcAll{1,1}.NVCind_resp), nanmean(calcAll{1,2}.NVCind_resp)],'y');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.NVCind_resp), nanmean(calcAll{1,2}.NVCind_resp)], ...
    [calcAll{1,1}.NVCind_resp_SEM, calcAll{1,2}.NVCind_resp_SEM],'y.');
xlim([0 3]);
title(['NVCindex, RESP data, p=',num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}],[grpNames{2}]});
%save plot into exp_dir as png
figsave = fullfile(fname,['NVCindex_respOnly.png']);
saveas(gcf, figsave);
close; 

clear p; 

%% save out the data across both groups into top dir
matfile = fullfile([fname,filesep,'dataCutByCMRO2_allGrps']);
save(matfile, 'calcAll', 'calcAll_sep', 'grpNames', '-v7.3');

