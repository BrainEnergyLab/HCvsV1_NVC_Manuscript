
%NB need to run following codes before this will work:
% - xyFWHM to extract the vessel data
% - cellSort and/or calcAutoROI to extract calcium data
% - cutByCalcEvents_xy to cut around individual calcium events
% - findAvgDiam_manual_xy to extract the rough diameter of the vessels
% sampled

%this script will also call the following functions:
% - findFolders
% - sortDataByResponsiveTrials

clear all; close all; 
%'D:\Dropbox (Brain Energy Lab)\Everything\Kira\2P\2P_Data_sortedByAnalysis\localCalcVessel_xy'
fname=cd; 

% important variables to send in for classifying the vessel as responsive
% or not:
%time to check diameter curve (to see if resp or not) - NB/ calc peak is
%aligned to 5s
preferences.respPeriod_diam = [5 10]; %seconds
%set how many standard deviations above thresh to deem as responsive
%aka the std multiplier
preferences.std_multiplier = 1; %SD - 1 or 1.5
%calc is 2, but the diameter peaks less high
%time in seconds to average either side of the max peak, i.e. to see if
%the peak is above the std dev of the baseline - this may need to be
%larger for slower responses (e.g. haem vs calc)
preferences.layer = 0; %1-if want to extract layer info, 0-if not 

%% find all the subfolders - i.e. groups - in the top dir
%only take the names which are less than 4 letters long, to exclude non
%label folders, so make sure your group names are short enough
listing=dir(fname);
for a = 1:size(listing,1)-2
    if size(listing(2+a).name,2) <= 6
        grpNames{a}=listing(2+a).name;
    else
        grpNames{a}=[];
    end
end
%remove the empty cells 
grpNames = grpNames(~cellfun('isempty',grpNames));

%loop the groups, extract the data
for n = 1:size(grpNames,2)
    
    clear findmatfiles; 
    %inside each group subfolder, find the mat file
    findmatfiles=findFolders([fname,filesep,grpNames{1,n}], ...
        'calcSort_xyFWHM.mat');
    %inform user of progress
    disp(['loading data for ', grpNames{1,n}, ' grp...']); 
    
    clear a; 
    for a = 1:size(findmatfiles,2) %loop through expDirs with data
        
        %find individual exp dir
        [expDir,matfile]=fileparts(findmatfiles{1,a});
        
        %.ini file contains parameter info (e.g. size of frame)
        find_ini_file = findFolders(expDir, '*.ini');
        x = cell2mat(find_ini_file);
        %convert info into a struct, which contains size info
        ini_file = ini2struct(x);
        %clear old unnecessary variables from workspace
        clear x find_ini_file;
        %should be in um from surface if person set up correctly when imaging:
        zPos = str2num(ini_file.x_.setz); %check this is correct??
        calcAll_sep{n}.zPos{a} = zPos; 
        
        %load data
        load(fullfile([expDir,filesep,matfile]));
        %load rough diameter
        load(fullfile([expDir,filesep,'roughdiam.mat']));
        %need time
        cd(expDir);
        load('contData_xyFWHM.mat', 'time');
        load('contData_xyFWHM.mat', 'cont_diam');
        diameter = cont_diam; clear cont_diam; 
        cd(fname);
        
        %get average diam:
        calcAll_sep{n}.diamUM{a} = roughdiam; 

        %also load and store animal name, layer info, and neuron type:
            %animal name/dir
            calcAll_sep{n}.expLabel{a} = expDir(91:end);
        if preferences.layer
            %layer info
            [ttt] = findFolders(expDir, 'layer_*.txt');
            t = extractAfter(ttt{1,1},'layer_');
            calcAll_sep{n}.layerInfo{a} = extractBefore(t,'.txt');
            clear t ttt;
        end
        
        %put data into a cell - keep them separate so you can check for any
        %outliers and know where they came from
        %resize time dim so all are same size
        %calcium
        calcAll_sep{n}.calc{a} = imresize(calcSort.calc, ...
            [size(calcSort.calc,1) 65]);
        %loco
        calcAll_sep{n}.loco_all{a} = imresize(calcSort.loco_all, ...
            [size(calcSort.loco_all,1) 65]);
        %stim
        calcAll_sep{n}.stim_all{a} = imresize(calcSort.stim_all, ...
            [size(calcSort.stim_all,1) 65]);
        %vessel diameter
        calcAll_sep{n}.diam{a} = imresize(calcSort.diam, ...
            [size(calcSort.diam,1) 65]);
        for g = 1:size(calcSort.diamBS,2)
            calcAll_sep{n}.diamBS{a}(:,g,:) = imresize(squeeze(calcSort.diamBS(:,g,:)), ...
                [size(calcSort.diamBS,1) 65]);
        end
        
        %time
        calcAll_sep{n}.time{a} = imresize(calcSort.time, [1 65]);
        %fps
        calcAll_sep{n}.fps{a} = fps;
        %grp label
        calcAll_sep{n}.grpLabel = grpNames{n}; 
        
        clear counter;
        %put all trials into one cont trace
        if a == 1 %if 1st trial create var
            calcAll{n}.counter = size(calcSort.diam,1);
            calcAll{n}.calc = imresize(calcSort.calc, ...
                [size(calcSort.calc,1) 65]);
            calcAll{n}.loco_all = imresize(calcSort.loco_all, ...
                [size(calcSort.loco_all,1) 65]);
            calcAll{n}.stim_all = imresize(calcSort.stim_all, ...
                [size(calcSort.stim_all,1) 65]);
            calcAll{n}.diam = imresize(calcSort.diam, ...
                [size(calcSort.diam,1) 65]);
            for g = 1:size(calcSort.diamBS,2)
                calcAll{n}.diamBS(:,g,:) = imresize(squeeze(calcSort.diamBS(:,g,:)), ...
                    [size(calcSort.diamBS,1) 65]);
            end
            calcAll{n}.time = calcAll_sep{1,n}.time{1,1};
            %take average diameter reading in um
            calcAll{n}.diam_avg = repmat(roughdiam,...
                [1 size(calcSort.calc,1)]);
%             %take z depth
%             calcAll{n}.zPos = repmat(zPos,...
%                 [1 size(calcSort.calc,1)]);
            
            %get labels out
            for b = 1:size(calcSort.calc,1)
                if preferences.layer
                    calcAll{n}.layerInfo{b} = calcAll_sep{n}.layerInfo{a};
                end
                calcAll{n}.expLabel{b} = calcAll_sep{n}.expLabel{a};
            end
            calcAll{n}.diam_std = repmat(nanstd(nanmean(diameter, 1)),...
                [1 size(calcSort.calc,1)]);
            calcAll{n}.fps = repmat(fps, [1 size(calcSort.calc,1)]);
            
        else %post-first trial, so merge with rest of data
            
            %resize all the calcium events so they are all same across
            %events
            calcAll{n}.counter = [calcAll{n}.counter, size(calcSort.diam,1)];
            calcAll{n}.calc = [calcAll{n}.calc; imresize(calcSort.calc, ...
                [size(calcSort.calc,1) 65])];
            calcAll{n}.loco_all = [calcAll{n}.loco_all; imresize(calcSort.loco_all, ...
                [size(calcSort.loco_all,1) 65])];
            calcAll{n}.stim_all = [calcAll{n}.stim_all; imresize(calcSort.stim_all, ...
                [size(calcSort.stim_all,1) 65])];
            calcAll{n}.diam = [calcAll{n}.diam; imresize(calcSort.diam, ...
                [size(calcSort.diam,1) 65])];
            for g = 1:size(calcSort.diamBS,2)
                ttt(:,g,:) = imresize(squeeze(calcSort.diamBS(:,g,:)), ...
                    [size(calcSort.diamBS,1) 65]);
            end
            calcAll{n}.diamBS = [calcAll{n}.diamBS; ttt];
            clear ttt;
            %take average diameter reading in um
            calcAll{n}.diam_avg = [calcAll{n}.diam_avg, ...
                repmat(roughdiam,...
                [1 size(calcSort.calc,1)])];         
%             calcAll{n}.zPos = [calcAll{n}.zPos, ...
%                 repmat(zPos,...
%                 [1 size(calcSort.calc,1)])];
            calcAll{n}.diam_std = [calcAll{n}.diam_std, ...
                repmat(nanstd(nanmean(diameter, 1)),...
                [1 size(calcSort.calc,1)])];
            calcAll{n}.fps = [calcAll{n}.fps, ...
                repmat(fps ,[1 size(calcSort.calc,1)])];
            %get labels out
            clear counter; 
            counter = size(calcAll{n}.expLabel,2);
            for b = 1:size(calcSort.calc,1)
                if preferences.layer
                    calcAll{n}.layerInfo{b+counter} = calcAll_sep{n}.layerInfo{a};
                end
                calcAll{n}.expLabel{b+counter} = calcAll_sep{n}.expLabel{a};
            end
        end
        calcAll{n}.grpLabel = grpNames{n}; 
        
        
    end %end of looping through expDirs containing needed data
    
    clearvars -except fname calcAll_sep calcAll grpNames preferences;
  
end %end of looping groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% compare diam sizes sampled from:
%i.e. check any differences are not just because you're sampling from
%different sizes of vessel

for b = 1:size(calcAll,2) %loop grps
    [calcAll{1,b}.diam_avg_SEM]=getSEM(calcAll{1,b}.diam_avg,2);
end
%plot as bar charts
figure;
bar(1:2,[nanmean(calcAll{1}.diam_avg,2), nanmean(calcAll{2}.diam_avg,2)],'b');
hold on;
errorbar(1:2,[nanmean(calcAll{1}.diam_avg,2), nanmean(calcAll{2}.diam_avg,2)],...
    [getSEM(calcAll{1}.diam_avg,2), getSEM(calcAll{2}.diam_avg,2)], 'b.');
ylabel('diameter, um'); 
name = {grpNames{1};grpNames{2}};
set(gca,'xticklabel',name);
[~,p] = ttest2(calcAll{1}.diam_avg,calcAll{2}.diam_avg);
title(['Compare vessel sizes sampled, p=',num2str(p)]);

%% ALL the data - i.e. not sorted yet by responsive vessels: 
%find AUC and other parameters for ALL the data across groups
for n = 1:size(grpNames,2)
    for a = 1:size(calcAll{1,n}.calc,1)
        
        %call function to assess the ts parameters for the calc data
        %take onset of calc until 2s post onset to assess data
        clear ttt; 
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.calc(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=3,1) ...
            find(calcAll{1,n}.time>=8,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_calc(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_calc(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_calc(a)=ttt.t2p;

        %call function to assess ts param for corresponding diam data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.diam(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
            find(calcAll{1,n}.time>=10,1)],1);
        %save out data of interest into struct
        calcAll{1,n}.AUC_vessel(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_vessel(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_vessel(a)=ttt.t2p;
        
        %bs 
        clear ttt;
            [ttt] = findtsParameters_calcPeak(squeeze(nanmean(calcAll{1,n}.diamBS(a,:,:),2))',...
                calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
                find(calcAll{1,n}.time>=10,1)],1);
            %save out data of interest into struct
            calcAll{1,n}.AUC_vesselBS(a)=ttt.AUC;
            calcAll{1,n}.maxPeak_vesselBS(a)=ttt.maxPeak;
            calcAll{1,n}.t2p_vesselBS(a)=ttt.t2p;
        
    end
end
clear ttt; 

% plot traces for ALL (non sorted) data:
%find size of screen (for figure sizing)
screenSz=get(0,'Screensize');

%get sem
for n = 1:size(grpNames,2)
    [calcAll{1,n}.calc_SEM]=getSEM(nanmean(calcAll{1,n}.calc,2),1);
    [calcAll{1,n}.calcTS_SEM]=getSEM(calcAll{1,n}.calc,1);
    [calcAll{1,n}.calcPeak_SEM]=getSEM(calcAll{1,n}.maxPeak_calc,2);
    [calcAll{1,n}.calcAUC_SEM]=getSEM(calcAll{1,n}.AUC_calc,2);
    [calcAll{1,n}.calct2p_SEM]=getSEM(calcAll{1,n}.t2p_calc,2);
end

%avg deltaD
for n = 1:size(grpNames,2)
    [calcAll{1,n}.diam_SEM]=getSEM(nanmean(calcAll{1,n}.diam,2),1);
    [calcAll{1,n}.diamTS_SEM]=getSEM(calcAll{1,n}.diam,1);
    [calcAll{1,n}.diamPeak_SEM]=getSEM(calcAll{1,n}.maxPeak_vessel,2);
    [calcAll{1,n}.diamAUC_SEM]=getSEM(calcAll{1,n}.AUC_vessel,2);
    [calcAll{1,n}.diamt2p_SEM]=getSEM(calcAll{1,n}.t2p_vessel,2);
    %BS
    [calcAll{1,n}.diamBS_peak_SEM]=getSEM(calcAll{1,n}.maxPeak_vesselBS,2);
    [calcAll{1,n}.diamBS_AUC_SEM]=getSEM(calcAll{1,n}.AUC_vesselBS,2);
    [calcAll{1,n}.diamBS_t2p_SEM]=getSEM(calcAll{1,n}.t2p_vesselBS,2);
end

% plot all traces, i.e. resp and non resp combined - ALL data
figure;
ax(1)=subplot(221);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.calc,1), 'g', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.calc,1), ...
    calcAll{1}.calcTS_SEM, 'lineProps','g');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.1 0.4],'k','LineWidth',2);
title([grpNames{1}, ' calc']);

mtit(gcf,'ALL data','fontsize', 14,'color', [1 0 0],'xoff',1.25,'yoff',0.25);

ax(2)=subplot(222);
plot(calcAll{1}.time, ...
    nanmean(calcAll{2}.calc,1), 'g', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.calc,1), ...
    calcAll{2}.calcTS_SEM, 'lineProps','g');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.1 0.4],'k','LineWidth',2);
title([grpNames{2}, ' calc']);

ax(3)=subplot(223);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.diam,1), 'r', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.diam,1), ...
    calcAll{1}.diamTS_SEM, 'lineProps','r');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.02 0.04],'k','LineWidth',2);
title([grpNames{1}, ' diam']);

ax(4)=subplot(224);
plot(calcAll{2}.time, ...
    nanmean(calcAll{2}.diam,1), 'r', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.diam,1), ...
    calcAll{2}.diamTS_SEM, 'lineProps','r');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.02 0.04],'k','LineWidth',2);
title([grpNames{2}, ' diam']);

linkaxes([ax(1),ax(2)],'xy');
linkaxes([ax(3),ax(4)],'xy');
clear ax; 

%plot the max peak values for calc and diam traces across ALL data:
figure;
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_calc), ...
    nanmean(calcAll{1,2}.maxPeak_calc)],'g');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_calc), ...
    nanmean(calcAll{1,2}.maxPeak_calc)], ...
    [calcAll{1,1}.calcPeak_SEM, ...
    calcAll{1,2}.calcPeak_SEM],'g.');
xlim([0 3]);
[~,p] = ttest2(calcAll{1}.maxPeak_calc,calcAll{2}.maxPeak_calc);
title(['All Data, DeltaF Avg Peak, p=',num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}, ', nEvents=',  ...
    num2str(size(calcAll{1,1}.calc,1))], ...
    [grpNames{2}, ', nEvents=', ...
    num2str(size(calcAll{1,2}.calc,1))]});
ylabel('\Delta/ F/F');
%save plot into exp_dir as png
figsave = fullfile(fname,['DeltaFAvgPeak_Barchart.png']);
saveas(gcf, figsave);

figure;
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_vessel), ...
    nanmean(calcAll{1,2}.maxPeak_vessel)],'r');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_vessel), ...
    nanmean(calcAll{1,2}.maxPeak_vessel)], ...
    [calcAll{1,1}.diamPeak_SEM, ...
    calcAll{1,2}.diamPeak_SEM],'r.');
xlim([0 3]);
[~,p] = ttest2(calcAll{1}.maxPeak_vessel,calcAll{2}.maxPeak_vessel);
title(['All Data, DeltaD Avg Peak, p=',num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}, ', nEvents=',  ...
    num2str(size(calcAll{1,1}.diam,1))], ...
    [grpNames{2}, ', nEvents=', ...
    num2str(size(calcAll{1,2}.diam,1))]});
ylabel('\Delta/ D/D');
%save plot into exp_dir as png
figsave = fullfile(fname,['DeltaDAvgPeak_Barchart.png']);
saveas(gcf, figsave);

%% classify the VESSELS as responsive or not:
%get most responsive index
%call function
for b = 1:size(calcAll,2)
    
    disp(['calculating response rates for real and bootstrap data, grp', num2str(b)]); 
    calcAll{b}.responsiveIndBS=nan(100,size(calcAll{1,b}.calc,1));
    for c = 1:size(calcAll{b}.calc,1) %calc events
        [calcAll{b}.responsiveInd(c)] = sortDataByResponsiveTrials ...
            (calcAll{b}.diam(c,:),calcAll{b}.time, preferences, 'diam');
        for g = 1:100
            [calcAll{b}.responsiveIndBS(g,c)] = sortDataByResponsiveTrials ...
            (squeeze(calcAll{b}.diamBS(c,g,:))',calcAll{b}.time, preferences, 'diam');
        end
    end
    
    %check for any NaNs and make them zeros -- thus rendering them as
    %'unresponsive' trials - NB/ will affect the % calc slightly as more
    %deemed unresponsive 
    %FUTURE - FIND OUT WHY NANS - GUESSING NO DATA IN DIAM??
    calcAll{b}.responsiveInd(isnan(calcAll{b}.responsiveInd))=0; 
    calcAll{b}.responsiveIndBS(find(isnan(calcAll{b}.responsiveIndBS)))=0; 
    
    %find % of responsive diameter changes from real data
    calcAll{b}.respPcnt = (sum(calcAll{b}.responsiveInd)/ ...
        size(calcAll{b}.diam,1))*100; 
    %from bootstrap data
    %calculate the % responsive across iterations for each calc event
    calcAll{b}.respPcntBS_all = (sum(calcAll{b}.responsiveIndBS')/ ...
        size(calcAll{b}.diam,1))*100;
   %mean across iterations per calc event, to get the overall mean respond
    calcAll{b}.respPcntBS_avg = nanmean(calcAll{b}.respPcntBS_all);
      
end


%% plot the average responsive calc and vessel traces vs non-resp - get errorbars

%get errorbars for the responsive vs non resp
for b = 1:size(calcAll,2)
    [calcAll{b}.calc_respSEM]=getSEM ...
        (calcAll{b}.calc(find(calcAll{b}.responsiveInd),:),1);
    [calcAll{b}.calc_nonrespSEM]=getSEM ...
        (calcAll{b}.calc(find(calcAll{b}.responsiveInd==0),:),1);
    [calcAll{b}.diam_respSEM]=getSEM ...
        (calcAll{b}.diam(find(calcAll{b}.responsiveInd),:),1);
    [calcAll{b}.diam_nonrespSEM]=getSEM ...
        (calcAll{b}.diam(find(calcAll{b}.responsiveInd==0),:),1);
end


%% responsive vessels only: 

%plots

%plot just responsive traces
figure;
ax(1)=subplot(221);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.calc(find(calcAll{1}.responsiveInd),:),1), ...
    'g', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.calc(find(calcAll{1}. ...
    responsiveInd),:),1), calcAll{1}.calc_respSEM, ...
    'lineProps','g');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.1 0.2],'k','LineWidth',2);
title([grpNames{1},' calc']);

ax(2)=subplot(222);
plot(calcAll{1}.time, ...
    nanmean(calcAll{2}.calc(find(calcAll{2}.responsiveInd),:),1), ...
    'g', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.calc(find(calcAll{2}. ...
    responsiveInd),:),1), calcAll{2}.calc_respSEM, ...
    'lineProps','g');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.1 0.2],'k','LineWidth',2);
title([grpNames{2},' calc']);

ax(3)=subplot(223);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.diam(find(calcAll{1}.responsiveInd),:),1), ...
    'r', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.diam(find(calcAll{1}. ...
    responsiveInd),:),1), calcAll{1}.diam_respSEM, ...
    'lineProps','r');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.01 0.18],'k','LineWidth',2);
title([grpNames{1},' diam']);

ax(4)=subplot(224);
plot(calcAll{2}.time, ...
    nanmean(calcAll{2}.diam(find(calcAll{2}.responsiveInd),:),1), ...
    'r', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.diam(find(calcAll{2}. ...
    responsiveInd),:),1), calcAll{2}.diam_respSEM, ...
    'lineProps','r');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.01 0.18],'k','LineWidth',2);
title([grpNames{2},' diam']);

linkaxes([ax(1),ax(2)],'xy');
linkaxes([ax(3),ax(4)],'xy');
clear ax; 

%save variables with just resp data
for n = 1:size(calcAll,2)
    calcAll{1,n}.calc_resp = calcAll{1,n}.calc(find(calcAll{1,n}. ...
    responsiveInd),:);
    calcAll{1,n}.diam_resp = calcAll{1,n}.diam(find(calcAll{1,n}. ...
    responsiveInd),:);
end

%ts parameters on responsive traces only 
%find AUC and other parameters for the data across groups
for n = 1:size(grpNames,2)
    for a = 1:size(calcAll{1,n}.calc_resp,1)
        
        %call function to assess the ts parameters for the calc data
        %take onset of calc until 2s post onset to assess data
        clear ttt; 
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.calc_resp(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=3,1) ...
            find(calcAll{1,n}.time>=8,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_calc_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_calc_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_calc_resp(a)=ttt.t2p;

        %call function to assess ts param for corresponding diam data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.diam_resp(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
            find(calcAll{1,n}.time>=10,1)],1);
        %save out data of interest into struct
        calcAll{1,n}.AUC_vessel_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_vessel_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_vessel_resp(a)=ttt.t2p;
        
    end
    %errorbars
    [calcAll{1,n}.AUC_vessel_resp_SEM]=getSEM(calcAll{1,n}.AUC_vessel_resp,2);
    [calcAll{1,n}.maxPeak_vessel_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_vessel_resp,2);
    [calcAll{1,n}.t2p_vessel_resp_SEM]=getSEM(calcAll{1,n}.t2p_vessel_resp,2);
    [calcAll{1,n}.AUC_calc_resp_SEM]=getSEM(calcAll{1,n}.AUC_calc_resp,2);
    [calcAll{1,n}.maxPeak_calc_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_calc_resp,2);
    [calcAll{1,n}.t2p_calc_resp_SEM]=getSEM(calcAll{1,n}.t2p_calc_resp,2);
end

%ts param for responsive vessel traces
%plot max peaks for responsive calc and vessels
figure;
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_calc_resp), ...
    nanmean(calcAll{1,2}.maxPeak_calc_resp)],'g');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_calc_resp), ...
    nanmean(calcAll{1,2}.maxPeak_calc_resp)], ...
    [calcAll{1,1}.maxPeak_calc_resp_SEM, ...
    calcAll{1,2}.maxPeak_calc_resp_SEM],'g.');
xlim([0 3]);
[~,p] = ttest2(calcAll{1}.maxPeak_calc_resp,calcAll{2}.maxPeak_calc_resp);
title(['Resp Data, Calc Max Peak, p=',num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}, ', nEvents=',  ...
    num2str(size(calcAll{1,1}.maxPeak_calc_resp,2))], ...
    [grpNames{2}, ', nEvents=', ...
    num2str(size(calcAll{1,2}.maxPeak_calc_resp,2))]});
ylabel('\Delta/ F/F');

figure;
bar(1:2,[nanmean(calcAll{1,1}.maxPeak_vessel_resp), ...
    nanmean(calcAll{1,2}.maxPeak_vessel_resp)],'r');
hold on;
errorbar(1:2,[nanmean(calcAll{1,1}.maxPeak_vessel_resp), ...
    nanmean(calcAll{1,2}.maxPeak_vessel_resp)], ...
    [calcAll{1,1}.maxPeak_vessel_resp_SEM, ...
    calcAll{1,2}.maxPeak_vessel_resp_SEM],'r.');
xlim([0 3]);
[~,p] = ttest2(calcAll{1}.maxPeak_vessel_resp,calcAll{2}.maxPeak_vessel_resp);
title(['Resp Data, Diam Max Peak, p=',num2str(p)]);
set(gca, 'XTickLabel', {[grpNames{1}, ', nEvents=',  ...
    num2str(size(calcAll{1,1}.maxPeak_calc_resp,2))], ...
    [grpNames{2}, ', nEvents=', ...
    num2str(size(calcAll{1,2}.maxPeak_calc_resp,2))]});
ylabel('\Delta/ D/D');

%% plot non responsive traces: 
%i.e. sanity check looks like right trials have been excluded

figure;
ax(1)=subplot(221);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.calc(find(calcAll{1}.responsiveInd==0),:),1), ...
    'g', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.calc(find(calcAll{1}. ...
    responsiveInd==0),:),1), calcAll{1}.calc_nonrespSEM, ...
    'lineProps','g');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.1 0.2],'k','LineWidth',2);
title([grpNames{1},' calc']);

ax(2)=subplot(222);
plot(calcAll{1}.time, ...
    nanmean(calcAll{2}.calc(find(calcAll{2}.responsiveInd==0),:),1), ...
    'g', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.calc(find(calcAll{2}. ...
    responsiveInd==0),:),1), calcAll{2}.calc_nonrespSEM, ...
    'lineProps','g');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.1 0.2],'k','LineWidth',2);
title([grpNames{2},' calc']);

ax(3)=subplot(223);
plot(calcAll{1}.time, ...
    nanmean(calcAll{1}.diam(find(calcAll{1}.responsiveInd==0),:),1), ...
    'r', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.diam(find(calcAll{1}. ...
    responsiveInd==0),:),1), calcAll{1}.diam_nonrespSEM, ...
    'lineProps','r');
plot([calcAll{1}.time(1) calcAll{1}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.01 0.05],'k','LineWidth',2);
title([grpNames{1},' diam']);

ax(4)=subplot(224);
plot(calcAll{2}.time, ...
    nanmean(calcAll{2}.diam(find(calcAll{2}.responsiveInd==0),:),1), ...
    'r', 'LineWidth',2);
% ylim([-0.05 0.5]);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.diam(find(calcAll{2}. ...
    responsiveInd==0),:),1), calcAll{2}.diam_nonrespSEM, ...
    'lineProps','r');
plot([calcAll{2}.time(1) calcAll{2}.time(end)],[0 0],'k','LineWidth',2)
plot([2 2],[-0.01 0.05],'k','LineWidth',2);
title([grpNames{2},' diam']);

linkaxes([ax(1),ax(2)],'xy');
linkaxes([ax(3),ax(4)],'xy');
clear ax; 

%% label the trials by neuron type, layer info using grp numbers (i.e. binarising)
% for b = 1:2 %loop regions
% 
%     clear ttt t;
%     ttt = calcAll{b}.layerInfo;
%     for t = 1:size(ttt,2)
%         %layerInfo
%         if b == 1 %HC
%             if regexp(cell2mat(ttt(t)), 'SO')
%                 layerInfoInd{b}(t) = 0;
%             elseif regexp(cell2mat(ttt(t)), 'SP')
%                 layerInfoInd{b}(t) = 1;
%             elseif regexp(cell2mat(ttt(t)), 'SR')
%                 layerInfoInd{b}(t) = 2;
%             elseif regexp(cell2mat(ttt(t)), 'SLM')
%                 layerInfoInd{b}(t) = 3;
%             end
%         else %V1
%             if calcAll{b}.zPos(t) >= -60
%                 layerInfoInd{b}(t) = 0;
%             elseif (calcAll{b}.zPos(t) < -60) && (calcAll{b}.zPos(t) >= -280) 
%                 layerInfoInd{b}(t) = 1;
%             elseif (calcAll{b}.zPos(t) < -280) && (calcAll{b}.zPos(t) >= -410) 
%                 layerInfoInd{b}(t) = 2;
%             else
%                 layerInfoInd{b}(t) = 3;
%             end
%         end
%     end
% end

%% extract bootstrap traces:

clear BS; 
for h = 1:size(calcAll,2) %region
    for i = 1:size(calcAll{h}.diamBS,2) %bs iterations
        
        clear data indx;
        data = squeeze(calcAll{h}.diamBS(:,i,:));
        indx = calcAll{h}.responsiveIndBS(i,:);
        
        BS{h}.iteration{i}.responsive = data(indx==1,:);
        BS{h}.iteration{i}.nResp = sum(indx==1);
        BS{h}.iteration{i}.nonresponsive = data(indx==0,:);
        BS{h}.iteration{i}.nNonResp = sum(indx==0);
        
    end
end
%take out of cells into matrices
for h = 1:size(calcAll,2) %region
    for i = 1:size(calcAll{h}.diamBS,2) %bs iterations
        
        if i == 1
            %responsive
            calcAll{h}.diamBS_resp = BS{h}.iteration{i}.responsive;
            calcAll{h}.diamBS_nRespTrials = BS{h}.iteration{i}.nResp; 
            %nonresponsive
            calcAll{h}.diamBS_nonresp = BS{h}.iteration{i}.nonresponsive;
            calcAll{h}.diamBS_nNonRespTrials = BS{h}.iteration{i}.nNonResp; 
        else
            %responsive
            calcAll{h}.diamBS_resp = [calcAll{h}.diamBS_resp; ...
                BS{h}.iteration{i}.responsive]; 
            calcAll{h}.diamBS_nRespTrials = [calcAll{h}.diamBS_nRespTrials; ...
                BS{h}.iteration{i}.nResp]; 
            %non-responsive
            calcAll{h}.diamBS_nonresp = [calcAll{h}.diamBS_nonresp; ...
                BS{h}.iteration{i}.nonresponsive]; 
            calcAll{h}.diamBS_nNonRespTrials = [calcAll{h}.diamBS_nNonRespTrials; ...
                BS{h}.iteration{i}.nNonResp]; 
        end
    
    end
end

%errorbars
for i = 1:size(calcAll,2)
    for j = 1:size(calcAll{i}.diamBS,2)
        [calcAll{i}.diamBS_SEM(j,:)] = getSEM(squeeze(calcAll{i}.diamBS(:,i,:)),1);
    end
    [calcAll{i}.diamBS_resp_SEM] = getSEM(calcAll{i}.diamBS_resp);
    [calcAll{i}.diamBS_nonresp_SEM] = getSEM(calcAll{i}.diamBS_nonresp);
end

%plot ALL BS traces - i.e. all traces, not separated by responsiveness
figure;
subplot(211); 
plot(calcAll{1}.time, squeeze(nanmean(nanmean(calcAll{1}.diamBS,1),2)), ...
    'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, squeeze(nanmean(nanmean(calcAll{1}.diamBS,1),2))', ...
    nanmean(calcAll{1}.diamBS_SEM,1), 'lineProps','r');
title(['HC, BS diam all, n=', num2str(size(calcAll{1}.diamBS,1)*size(calcAll{1}.diamBS,2))]); 
xlim([0 18]); 
ylim([-0.01 0.04]); 
subplot(212); 
plot(calcAll{2}.time, squeeze(nanmean(nanmean(calcAll{2}.diamBS,1),2)), ...
    'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, squeeze(nanmean(nanmean(calcAll{2}.diamBS,1),2))', ...
    nanmean(calcAll{2}.diamBS_SEM,1), 'lineProps','r');
title(['V1, BS diam all, n=', num2str(size(calcAll{2}.diamBS,1)*size(calcAll{2}.diamBS,2))]); 
xlim([0 18]); 
ylim([-0.01 0.04]); 


%plot responsive vs non responsive BS traces: 
figure;
subplot(221);
plot(calcAll{1}.time, nanmean(calcAll{1}.diamBS_resp,1), ...
    'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.diamBS_resp,1), ...
    calcAll{1}.diamBS_resp_SEM, 'lineProps','r');
title([grpNames{1},', BS diam responsive, n=', num2str(size(calcAll{1}.diamBS_resp,1))]); 
xlim([0 18]); 
ylim([-0.01 0.09]); 

subplot(222);
plot(calcAll{1}.time, nanmean(calcAll{1}.diamBS_nonresp,1), ...
    'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{1}.time, nanmean(calcAll{1}.diamBS_nonresp,1), ...
    calcAll{1}.diamBS_nonresp_SEM, 'lineProps','r');
title([grpNames{1},', BS diam non-responsive, n=', num2str(size(calcAll{1}.diamBS_nonresp,1))]); 
xlim([0 18]); 
ylim([-0.01 0.09]); 

subplot(223);
plot(calcAll{2}.time, nanmean(calcAll{2}.diamBS_resp,1), ...
    'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.diamBS_resp,1), ...
    calcAll{2}.diamBS_resp_SEM, 'lineProps','r');
title([grpNames{2},', BS diam responsive, n=', num2str(size(calcAll{2}.diamBS_resp,1))]); 
xlim([0 18]); 
ylim([-0.01 0.09]); 

subplot(224);
plot(calcAll{2}.time, nanmean(calcAll{2}.diamBS_nonresp,1), ...
    'r', 'LineWidth',2);
hold on;
shadedErrorBar(calcAll{2}.time, nanmean(calcAll{2}.diamBS_nonresp,1), ...
    calcAll{2}.diamBS_nonresp_SEM, 'lineProps','r');
title([grpNames{2},', BS diam non-responsive, n=', num2str(size(calcAll{2}.diamBS_nonresp,1))]); 
xlim([0 18]); 
ylim([-0.01 0.09]); 

%get ts parameters for BS traces - responsive
%find AUC and other parameters for the data across groups
for n = 1:size(grpNames,2)
    for a = 1:size(calcAll{1,n}.calc_resp,1)
        
        %diam data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.diamBS_resp(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
            find(calcAll{1,n}.time>=8,1)],1);
        %save out data of interest into struct
        calcAll{1,n}.AUC_vesselBS_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_vesselBS_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_vesselBS_resp(a)=ttt.t2p;
        
    end
    %errorbars
    [calcAll{1,n}.AUC_vesselBS_resp_SEM]=getSEM(calcAll{1,n}.AUC_vesselBS_resp,2);
    [calcAll{1,n}.maxPeak_vesselBS_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_vesselBS_resp,2);
    [calcAll{1,n}.t2p_vesselBS_resp_SEM]=getSEM(calcAll{1,n}.t2p_vesselBS_resp,2);
    
end

%% frequency of responsiveness for real and shuffled data across grps

%plot
y = [(size(calcAll{1, 1}.diam_resp,1)/size(calcAll{1, 1}.calc,1))*100 ...
    100 - ((size(calcAll{1, 1}.diam_resp,1)/size(calcAll{1, 1}.calc,1))*100); ...
    
    calcAll{1}.respPcntBS_avg ...
    100-calcAll{1}.respPcntBS_avg; ...
    
    (size(calcAll{2}.diam_resp,1)/size(calcAll{2}.calc,1))*100 ...
    100 - ((size(calcAll{2}.diam_resp,1)/size(calcAll{2}.calc,1))*100);
    
    calcAll{2}.respPcntBS_avg ...
    100-calcAll{2}.respPcntBS_avg]; 

datalabels = categorical({['real',grpNames{1}],['shuffled',grpNames{1}],...
    ['real',grpNames{2}],['shuffled',grpNames{2}]});
datalabels = reordercats(datalabels,{['real',grpNames{1}],['shuffled',grpNames{1}],...
    ['real',grpNames{2}],['shuffled',grpNames{2}]});

figure;
bar(datalabels, y, 'stacked'); 
legend('R','NR'); 
ylim([0 110]); ylabel('%'); 
title('Freq of Resp'); 


%% categorise the responsive vs non-responsive data, i.e. by rest/loco, 
% stim/VR, layer, neuron type, vessel size 

%% vessel size
for b = 1:size(calcAll,2) % regions
    for c = 1:size(calcAll{b}.diam_avg,2) %calc events
        %binarise the vessel sizes - 0 = <10um, 1 = >= 10um
        if calcAll{b}.diam_avg(c) < 7 %small vessels
            calcAll{b}.vessSz(c)=0; %label 0
        elseif calcAll{b}.diam_avg(c) >=7 && calcAll{b}.diam_avg(c)<10 %middle vessels
            calcAll{b}.vessSz(c)=1; %label 1 or 2 depending on grouping
        else %largest vess
            calcAll{b}.vessSz(c)=2;
        end %end of check vess sizes
    end %end of looping calc events
end %end of looping brain regions 

%plot
y = [(sum(calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==0)/sum(calcAll{1}.vessSz==0))*100 ...
    100 - ((sum(calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==0)/sum(calcAll{1}.vessSz==0))*100); ...
    
    (sum(calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==1)/sum(calcAll{1}.vessSz==1))*100 ...
    100 - ((sum(calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==1)/sum(calcAll{1}.vessSz==1))*100);
    
    (sum(calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==2)/sum(calcAll{1}.vessSz==2))*100 ...
    100 - ((sum(calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==2)/sum(calcAll{1}.vessSz==2))*100);
    
    (sum(calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==0)/sum(calcAll{2}.vessSz==0))*100 ... 
    100 - ((sum(calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==0)/sum(calcAll{2}.vessSz==0))*100);
    
    (sum(calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==1)/sum(calcAll{2}.vessSz==1))*100 ... 
    100 - ((sum(calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==1)/sum(calcAll{2}.vessSz==1))*100);
    
    (sum(calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==2)/sum(calcAll{2}.vessSz==2))*100 ...
    100 - ((sum(calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==2)/sum(calcAll{2}.vessSz==2))*100)]; 

datalabels = categorical({['smallVess',grpNames{1}],['medVess',grpNames{1}],...
    ['largeVess',grpNames{1}],['smallVess',grpNames{2}],['medVess',grpNames{2}],...
    ['largeVess',grpNames{2}]});
datalabels = reordercats(datalabels,{['smallVess',grpNames{1}],...
    ['medVess',grpNames{1}],['largeVess',grpNames{1}], ...
    ['smallVess',grpNames{2}],['medVess',grpNames{2}],['largeVess',grpNames{2}]});

figure;
bar(datalabels, y, 'stacked'); 
legend('R','NR'); 
ylim([0 110]); ylabel('%'); 
title('Vessel Size'); 

%get time series stuff out:
for n = 1:size(grpNames,2)
    %diam
    calcAll{n}.diam_resp_small = calcAll{n}.diam_resp(find( ...
    calcAll{n}.vessSz(find(calcAll{n}.responsiveInd))==0),:);
    calcAll{n}.diam_resp_med = calcAll{n}.diam_resp(find( ...
    calcAll{n}.vessSz(find(calcAll{n}.responsiveInd))==1),:);
    calcAll{n}.diam_resp_big = calcAll{n}.diam_resp(find( ...
    calcAll{n}.vessSz(find(calcAll{n}.responsiveInd))==2),:);
    %calc
    calcAll{n}.calc_resp_small = calcAll{n}.calc_resp(find( ...
    calcAll{n}.vessSz(find(calcAll{n}.responsiveInd))==0),:);
    calcAll{n}.calc_resp_med = calcAll{n}.calc_resp(find( ...
    calcAll{n}.vessSz(find(calcAll{n}.responsiveInd))==1),:);
    calcAll{n}.calc_resp_big = calcAll{n}.calc_resp(find( ...
    calcAll{n}.vessSz(find(calcAll{n}.responsiveInd))==2),:);
end

%get ts parameters for BS traces - responsive
%find AUC and other parameters for the data across groups
for n = 1:size(grpNames,2)
    
    for a = 1:size(calcAll{1,n}.diam_resp_small,1)
        %diam data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.diam_resp_small(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
            find(calcAll{1,n}.time>=10,1)],1);
        %save out data of interest into struct
        calcAll{1,n}.AUC_diamSmall_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_diamSmall_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_diamSmall_resp(a)=ttt.t2p;
        %calc data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.calc_resp_small(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=4,1) ...
            find(calcAll{1,n}.time>=6,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_calcSmall_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_calcSmall_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_calcSmall_resp(a)=ttt.t2p;
    end
    %errorbars
    [calcAll{1,n}.AUC_diamSmall_resp_SEM]=getSEM(calcAll{1,n}.AUC_diamSmall_resp,2);
    [calcAll{1,n}.maxPeak_diamSmall_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_diamSmall_resp,2);
    [calcAll{1,n}.t2p_diamSmall_resp_SEM]=getSEM(calcAll{1,n}.t2p_diamSmall_resp,2);
    [calcAll{1,n}.AUC_calcSmall_resp_SEM]=getSEM(calcAll{1,n}.AUC_calcSmall_resp,2);
    [calcAll{1,n}.maxPeak_calcSmall_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_calcSmall_resp,2);
    [calcAll{1,n}.t2p_calcSmall_resp_SEM]=getSEM(calcAll{1,n}.t2p_calcSmall_resp,2);
    
    %medium
    for a = 1:size(calcAll{1,n}.diam_resp_med,1)
        %diam data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.diam_resp_med(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
            find(calcAll{1,n}.time>=10,1)],1);
        %save out data of interest into struct
        calcAll{1,n}.AUC_diamMed_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_diamMed_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_diamMed_resp(a)=ttt.t2p;
        %calc data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.calc_resp_med(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=4,1) ...
            find(calcAll{1,n}.time>=6,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_calcMed_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_calcMed_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_calcMed_resp(a)=ttt.t2p;
    end
    %errorbars
    [calcAll{1,n}.AUC_diamMed_resp_SEM]=getSEM(calcAll{1,n}.AUC_diamMed_resp,2);
    [calcAll{1,n}.maxPeak_diamMed_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_diamMed_resp,2);
    [calcAll{1,n}.t2p_diamMed_resp_SEM]=getSEM(calcAll{1,n}.t2p_diamMed_resp,2);
    [calcAll{1,n}.AUC_calcMed_resp_SEM]=getSEM(calcAll{1,n}.AUC_calcMed_resp,2);
    [calcAll{1,n}.maxPeak_calcMed_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_calcMed_resp,2);
    [calcAll{1,n}.t2p_calcMed_resp_SEM]=getSEM(calcAll{1,n}.t2p_calcMed_resp,2);
    
    for a = 1:size(calcAll{1,n}.diam_resp_big,1)
        %diam data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.diam_resp_big(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=5,1) ...
            find(calcAll{1,n}.time>=10,1)],1);
        %save out data of interest into struct
        calcAll{1,n}.AUC_diamBig_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_diamBig_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_diamBig_resp(a)=ttt.t2p;
        %calc data
        %take onset of calc until 5s post onset (as haem is slow)
        clear ttt;
        [ttt] = findtsParameters_calcPeak(calcAll{1,n}.calc_resp_big(a,:),...
            calcAll{1,n}.time,calcAll{1,n}.fps(a),[find(calcAll{1,n}.time>=4,1) ...
            find(calcAll{1,n}.time>=6,1)]);
        %save out data of interest into struct
        calcAll{1,n}.AUC_calcBig_resp(a)=ttt.AUC;
        calcAll{1,n}.maxPeak_calcBig_resp(a)=ttt.maxPeak;
        calcAll{1,n}.t2p_calcBig_resp(a)=ttt.t2p;
    end
    %errorbars
    [calcAll{1,n}.AUC_diamBig_resp_SEM]=getSEM(calcAll{1,n}.AUC_diamBig_resp,2);
    [calcAll{1,n}.maxPeak_diamBig_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_diamBig_resp,2);
    [calcAll{1,n}.t2p_diamBig_resp_SEM]=getSEM(calcAll{1,n}.t2p_diamBig_resp,2);
    [calcAll{1,n}.AUC_calcBig_resp_SEM]=getSEM(calcAll{1,n}.AUC_calcBig_resp,2);
    [calcAll{1,n}.maxPeak_calcBig_resp_SEM]=getSEM(calcAll{1,n}.maxPeak_calcBig_resp,2);
    [calcAll{1,n}.t2p_calcBig_resp_SEM]=getSEM(calcAll{1,n}.t2p_calcBig_resp,2);
    
end

%plot the diameter traces for responsive small and large vessels, per
%region

%calc
figure;

subplot(221); 
plot(calcAll{1}.time, nanmean(calcAll{1}.calc_resp(find( ...
    calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==0),:),1),'b');
hold on;
plot(calcAll{1}.time, nanmean(calcAll{1}.calc_resp(find( ...
    calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==1),:),1),'r');
plot(calcAll{1}.time, nanmean(calcAll{1}.calc_resp(find( ...
    calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==2),:),1),'g');
legend('smallVess','medVess','largeVess','Autoupdate','Off');
plot([calcAll{1}.time(1) calcAll{1}.time(end)], [0 0], 'k', 'LineWidth',2);
plot([5 5], [-0.01 0.03], 'k', 'LineWidth',0.5);
title([grpNames{1},', calc']); 
xlim([0 18]); 
ylim([-0.01 0.4]); 

subplot(222); 
plot(calcAll{2}.time, nanmean(calcAll{2}.calc_resp(find( ...
    calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==0),:),1),'b');
hold on;
plot(calcAll{2}.time, nanmean(calcAll{2}.calc_resp(find( ...
    calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==1),:),1),'r');
plot(calcAll{2}.time, nanmean(calcAll{2}.calc_resp(find( ...
    calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==2),:),1),'g');
legend('smallVess','medVess','largeVess','Autoupdate','Off');
plot([calcAll{2}.time(1) calcAll{2}.time(end)], [0 0], 'k', 'LineWidth',2);
plot([5 5], [-0.01 0.03], 'k', 'LineWidth',0.5);
title([grpNames{2},', calc']); 
xlim([0 18]); 
ylim([-0.01 0.4]); 

%vessels
subplot(223); 
plot(calcAll{1}.time, nanmean(calcAll{1}.diam_resp(find( ...
    calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==0),:),1),'b');
hold on;
plot(calcAll{1}.time, nanmean(calcAll{1}.diam_resp(find( ...
    calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==1),:),1),'r');
plot(calcAll{1}.time, nanmean(calcAll{1}.diam_resp(find( ...
    calcAll{1}.vessSz(find(calcAll{1}.responsiveInd))==2),:),1),'g');
legend('smallVess','medVess','largeVess','Autoupdate','Off');
plot([calcAll{1}.time(1) calcAll{1}.time(end)], [0 0], 'k', 'LineWidth',2);
plot([5 5], [-0.01 0.045], 'k', 'LineWidth',0.5);
title([grpNames{1},', diam']); 
xlim([0 18]); 
ylim([-0.05 0.2]); 

subplot(224); 
plot(calcAll{2}.time, nanmean(calcAll{2}.diam_resp(find( ...
    calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==0),:),1),'b');
hold on;
plot(calcAll{2}.time, nanmean(calcAll{2}.diam_resp(find( ...
    calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==1),:),1),'r');
plot(calcAll{2}.time, nanmean(calcAll{2}.diam_resp(find( ...
    calcAll{2}.vessSz(find(calcAll{2}.responsiveInd))==2),:),1),'g');
legend('smallVess','medVess','largeVess','Autoupdate','Off');
plot([calcAll{2}.time(1) calcAll{2}.time(end)], [0 0], 'k', 'LineWidth',2);
plot([5 5], [-0.01 0.045], 'k', 'LineWidth',0.5);
title([grpNames{2},', diam']); 
xlim([0 18]); 
ylim([-0.05 0.2]); 


%save out the data across both groups into top dir
matfile = fullfile([fname,filesep,'dataCutByCalc_xyFWHM_allGrps']);
save(matfile, 'calcAll', 'calcAll_sep', 'grpNames', '-v7.3');


