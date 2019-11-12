
function cutByCalcEvents_xy_LONGER(fname, prefs)
%
%function by Kira, May 2018, update June 2019
%written to detect all calcium events deemed as 'responsive' (>2std above
%pre-calc event baseline), and then also cut the haem during this calc
%increase - to see if local calc increase causes any vascular changes
%
%INPUTS-
%fname = top directory which contains calcAutoRoi.mat (net calc) file and
%contData_xyFWHM.mat (vessel diam from xy movie) file, NB/ the data should
%be separated by group (e.g. brain region, or APOE genotype), and you
%should be within the group directory (e.g. V1)
%prefs = 
%
%OUTPUTS-
%saves figures to show average calc peaks and corresponding diam, and saves
%mat file with data in inside individual expDirs - run later code to
%compare across groups
%
%Useful functions needed for code to run:
%findFolders, xyDiamStruct2contTrace, findLocoEvents, 
%sortDataByResponsiveTrials

if nargin < 2
    %auto set the prefs if not specified
    %# second between distinct calc events
    prefs.minDist = 0.5; %seconds, will later be converted to frames
    %plot the calc trace with onsets and offsets marked
    prefs.plotFlag = 0;
    %flag for if want to remove the flickers in the signal (not loco as
    %shorter than 1s, but may be calc as this is short) 1=remove, 0=keep 
    prefs.flickerFlag = 0;
    %for cutting around calc events - need the pre and post calc time
    %to cut
    prefs.base = [0 5]; %seconds, will later be converted to frames 0 3
    prefs.post = [5 12]; %seconds, will later be converted to frames 3 10
    %for responsive trial sorter- 
    %for classifying trials as responsive or not responsive - how many
    %standard deviations above baseline should peak of response be
    prefs.std_multiplier = 2; %std
    %for responsive trial sorter-
    prefs.respPeriod_calc = [4 6]; %seconds - 1s either side of peak 2 4
%     %check for calc response within 2s of detected onset 2 6
%     prefs.respPeriod_diam = [4 8]; %seconds
    %for responsive trial sorter- 
    %time in seconds to average either side of the max peak, i.e. to see if
    %the peak is above the std dev of the baseline - this may need to be
    %larger for slower responses (e.g. haem vs calc)
    prefs.avgPeak = 0.1; %seconds
    %if want to smooth vessel diam trace
    prefs.smoothVessFlag = 0; %1 to smooth, 0 to turn off smooth
end

%find all the data which has been run through net calc extractor
% [findmatfile]=findFolders(fname,'calcAutoRoi.mat');
[findmatfile]=findFolders(fname,'cell_sig_norm.mat');

%get screen size for plots
screenSz=get(0,'Screensize');

for n = 1:size(findmatfile,2)

    clearvars -except fname prefs findmatfile screenSz n 
    
    %find each session experimental directory containing net calc data
    [expDir,~]=fileparts(findmatfile{1,n});
    
    %inform user of progress
    disp(['processing file ',num2str(n),'/', ...
        num2str(size(findmatfile,2))]);
    disp(expDir); 
    
    %check the folder also has the xy diameter data extracted 
    if exist([expDir, filesep, 'contData_xyFWHM.mat'],'file')
        
        %load the data
        disp('loading data...'); 
        load([expDir, filesep, 'contData_xyFWHM.mat']);
        load([expDir, filesep, 'cell_sig_norm.mat']);
        load([expDir, filesep, 'stimLocoTraces.mat']);

        %check if the experiment is a stimulation experiment
        %NB stim needs to be in the name of the exp dir folder
        stimExpFlag = regexpi(expDir,'stim');
        
        %convert some of the preferences to frames, using the now loaded
        %fps
        prefs1.minDist = round(prefs.minDist * fps); 
        prefs1.base = round(prefs.base(2) * fps); 
        prefs1.post = round(prefs.post(2) * fps); 
        prefs1.plotFlag = prefs.plotFlag;
        prefs1.flickerFlag = prefs.flickerFlag;
        
        %smooth diam data - check this is okay to do
        if prefs.smoothVessFlag == 1
            for h = 1:size(diameter,1)
                diameter = smooth(cont_diam(h,:));
            end
        else
            diameter = cont_diam;
        end
        %search for any rogue zeros, and set these to NaNs, i.e. so they
        %don't effect the diam reading
        diameter(diameter==0)=NaN; 


        
        %% Bootstrap linear model
        %%%%%%%%% bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %also create multiple reshuffled traces of the vasc diam
        
        %reshuffle the diameter data - keep the calc as is, so have same calc
        %spikes detected
        nIts = 100;  % Number of repeats of randomization
        %Initialize  variable to save traces in
        shuffledDiamTrace = nan([nIts,size(diameter,2)]);
        
        for i = 1:nIts
            
            diamTraceTemp = nanmean(diameter,1);
            
            %Shuffle the rows (using randperm), and then stitch back together
            shuffledDiamTrace(i,:) = ...
                reshape(diamTraceTemp(:,randperm(size(diamTraceTemp,2))),[],1);
            
        end % Number of iterations
        
        
        %% detect calcium spikes

        %normalise calcium trace - remove drift etc.
        calcium_norm = calcium; 
        calcium_norm = calcium_norm - nanmin(calcium_norm);
        calcium_norm = calcium_norm / nanmax(calcium_norm);
        calcium_norm = calcium_norm - smooth(calcium_norm,100)';
        
        [calcEvents] = findSignalSpikesInTrace(calcium_norm, fps);
        
        %plot trace with detected spikes
        figure;
        %make fig 1/2 size of screen
        set(gcf, 'Position', [screenSz(1) screenSz(2) ...
            screenSz(3)/2 screenSz(4)]);
        plot(time, calcium, 'g');
        hold on;
%         plot(time(:,calcEvents(1,:)),calcium(:,calcEvents(1,:)), 'ko'); %start of peak
        plot(time(:,calcEvents(2,:)),calcium(:,calcEvents(2,:)), 'bo'); %peak
        title('calcium trace with peaks');
        legend('calc','peak peak'); 
        xlabel('time(s)'); %ylabel('\Delta/F/F');
        saveas(gcf, fullfile([expDir,filesep,'LONGER_peakFinder.png']));
        close;

        %check if any calcium events were detected
        if isempty(calcEvents) == 1
            
            %inform user that there are no loco events in this specific exp
            disp(['no calc events: ', expDir]);
            
        else %if there are loco events detected, proceed...
            
            %check there are no zeros in calcEvents, if so remove
            zeroIndSearch=find(calcEvents(1,:)==0);
            if ~isempty(zeroIndSearch)
                calcEvents(:,zeroIndSearch)=[];
            end
            
            figure;
            %make fig size of screen
            set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                screenSz(3) screenSz(4)]);
            ax1=subplot(311);
            plot(calcium, 'g');
            hold on;
            %plot calc events detected on calc trace
            %find the value of the calc data at the time point of the calc
            %event
            plot(calcEvents(2,:), calcium(1,calcEvents(2,:)), ...
                'ko', 'LineWidth',2);
            title('calcium, all detected events');
            xlabel('frames');
            ax2=subplot(312);
            plot(nanmean(diameter,1),'r');
            hold on;
            %plot calc events detected on diam trace
            plot(calcEvents(2,:), nanmean(diameter(:,calcEvents(2,:)), ...
                1),'ko', 'LineWidth',2);
            title('diameter');
            xlabel('frames'); ylabel('um');
            ax3=subplot(313);
            plot(locomotion);
            xlabel('frames');
            linkaxes([ax1,ax2,ax3],'x');
            saveas(gcf, fullfile([expDir,filesep, ...
                'LONGER_contPlots_calcEvents.png']));
            close;
            
            %check for calc events having enough time before and after for
            %cutting around them
            enoughTimInd=0; %preallocate the index with a zero
            for c = 1:size(calcEvents,2) %loop the calc events
                %check enough time to cut baseline before the calc peak
                if (calcEvents(2,c)-round(prefs.base(2)*fps)+1)<=0
                    %put the calc event onset frame into ind if not
                    %enough base time
                    if sum(enoughTimInd) == 0 %check if any frames in ind
                        %if no frames in ind, then put in first slot
                        enoughTimInd(1) = c;
                    else %already frames in ind
                        %put into next slot
                        enoughTimInd(size(enoughTimInd,2)+1) = c;
                    end %end of check if already frames in ind
                end %end of check if enough base time
                %check enough time to cut post-event period post-calc
                %onset:
                if (calcEvents(2,c)+round((prefs.post(2)*fps)- ...
                        (prefs.post(1)*fps))+1)>=size(calcium,2)
                    if sum(enoughTimInd) == 0 %check if any frames in ind
                        %if no frames in ind, then put in first slot
                        enoughTimInd(1) = c;
                    else %already frames in ind
                        %put into next slot
                        enoughTimInd(size(enoughTimInd,2)+1) = c;
                    end %end of check if already frames in ind
                end %end of check if enough post-stim time 
            end %end of loop calc events
            %remove the calcium events without enough time for cutting:
            if sum(enoughTimInd) > 0
                calcEvents(:,enoughTimInd) = [];
            end
            clear enoughTimInd; 
            
            %cut around calcium peak - so peaks are all aligned
            for c = 1:size(calcEvents,2)
                %calcium trace
                calcCalcCut(c,:) = calcium(:,calcEvents(2,c)-round( ...
                    prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                    prefs.post(2)*fps)- (prefs.post(1)*fps))+1 ); 
                %diameter trace
                diamCalcCut(c,:,:) = diameter(:,calcEvents(2,c)-round( ...
                    prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                    prefs.post(2)*fps)- (prefs.post(1)*fps))+1 ); 
                %bootstrap diam
                diamBSCalcCut(c,:,:) = shuffledDiamTrace(:,calcEvents(2,c)-round( ...
                    prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                    prefs.post(2)*fps)- (prefs.post(1)*fps))+1); 
                %loco trace
                locoCalcCut(c,:) = locomotion(:,calcEvents(2,c)-round( ...
                    prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                    prefs.post(2)*fps)- (prefs.post(1)*fps))+1 ); 
                %stim trace
                if stimExpFlag > 1
                    stimCalcCut(c,:) = stim(:,calcEvents(2,c)-round( ...
                    prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                    prefs.post(2)*fps)- (prefs.post(1)*fps))+1 ); 
                else
                    if exist('VR')
                        stimCalcCut(c,:) = VR(:,calcEvents(2,c)-round( ...
                            prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                            prefs.post(2)*fps)- (prefs.post(1)*fps))+1 );
                    else
                        stimCalcCut(c,:) = zeros(size(locoCalcCut(c,:)));
                    end
                end
            end
            timeCalcCut = [1:size(calcCalcCut,2)]/fps; 
            
            %normalise data for first 1s -i.e. way before peak of calc:
            %calcium
            for c = 1:size(calcCalcCut,1)
                [calcCalcCut(c,:,:)]=findDeltaF(calcCalcCut(c,:), ...
                    [round(2*fps) round(3*fps)]); 
            end
            %diameter
            for c = 1:size(diamCalcCut,1)
%                 [diamCalcCut(c,:,:)]=findDeltaF(squeeze(diamCalcCut(c,:,:)), ...
%                     [round(2*fps) round(3*fps)]); 
                  [diamCalcCut(c,:,:)]=findDeltaF(squeeze(diamCalcCut(c,:,:)), ...
                    [1 round(5*fps)]); 
            end
            %diam bootstrap
            for c = 1:size(diamBSCalcCut,1)
%                 [diamBSCalcCut(c,:,:,:)]=findDeltaF(squeeze(diamBSCalcCut(c,:,:,:)), ...
%                     [round(2*fps) round(3*fps)]); 
                [diamBSCalcCut(c,:,:,:)]=findDeltaF(squeeze(diamBSCalcCut(c,:,:,:)), ...
                    [1 round(5*fps)]); 
            end
            
            %get most responsive index - do this on the normalised trace,
            %as sloping trend could mean some calc events are deemed not
            %sig due to decreasing baseline
            %call function
            clear responsiveInd;
            for c = 1:size(calcCalcCut,1)
                [responsiveInd(c)] = sortDataByResponsiveTrials ...
                    (calcCalcCut(c,:),time,prefs,'calc');
            end

            %plot the trials which are categorised as responsive and non
            %responsive and save into exp dir
            figure;
            %make fig size of screen
            set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                screenSz(3)/2 screenSz(4)]);
            for c = 1:size(calcCalcCut,1) %loop calcium trials
                %plot
                subplot(ceil(size(calcCalcCut,1)/4),4,c)
                plot(timeCalcCut, calcCalcCut(c,:),'g');
                hold on;
                clear yLimits;
                yLimits=get(gca, 'YLim');
                plot([timeCalcCut(1) timeCalcCut(end)],[0 0],'k');
                 plot([round(prefs.post(1)) round(prefs.post(1))], ...
                    [yLimits(1) yLimits(2)],'k--'); %peaks aligned
                plot([round(prefs.respPeriod_calc(1)) ...
                    round(prefs.respPeriod_calc(1))], [yLimits(1) ...
                    yLimits(2)],'r'); %responsive period checked - on
                plot([round(prefs.respPeriod_calc(2)) ...
                    round(prefs.respPeriod_calc(2))], [yLimits(1) ...
                    yLimits(2)],'r'); %responsive period checked - off
                xlim([timeCalcCut(1) timeCalcCut(end)]);
                xlabel('time(s)'); ylabel('\DeltaF/F');
                if responsiveInd(c) == 1
                    title('responsive');
                elseif responsiveInd(c) == 0
                    title('not responsive');
                elseif isnan(responsiveInd(c))
                    title('no data');
                end
            end %end of looping calcium trials
            saveas(gcf, fullfile([expDir,filesep, ...
                'LONGER_categoriseMostResponsiveCalcEvents_', ... 
                num2str(prefs.std_multiplier),'std.png']));
            close;
            
            %check if any of the trials have been classed as responsive 
            if nansum(responsiveInd) > 0
                
                %redo the calcium and diameter traces with only the
                %responsive events marked:
                figure;
                %make fig size of screen
                set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                    screenSz(3) screenSz(4)]);
                ax1=subplot(211);
                plot(calcium, 'g');
                hold on;
                %plot calc events detected on calc trace
                %find the value of the calc data at the time point of the calc
                %event
                plot(calcEvents(2,responsiveInd==1), calcium(:, ...
                    calcEvents(2,responsiveInd==1)),'ko', 'LineWidth',2);
                title('calcium, responsive events only');
                xlabel('frames');
                ax2=subplot(212);
                plot(nanmean(diameter,1),'r');
                hold on;
                %plot calc events detected on diam trace
                plot(calcEvents(2,responsiveInd==1), nanmean(diameter(:, ...
                    calcEvents(2,responsiveInd==1)),1),'ko', 'LineWidth',2);
                title('diameter');
                xlabel('frames'); ylabel('um');
                linkaxes([ax1,ax2],'x');
                saveas(gcf, fullfile([expDir,filesep, ...
                    'LONGER_contPlots_calcEvents_responsive.png']));
                close;
                
                %remove any rows which have a NaN for responsive ind
                nancheck = isnan(responsiveInd);
                if sum(nancheck) > 0 
                    responsiveInd(find(nancheck))=[]; 
                    calcCalcCut(find(nancheck),:)=[]; 
                    diamCalcCut(find(nancheck),:,:)=[]; 
                    diamBSCalcCut(find(nancheck),:,:)=[]; 
                    locoCalcCut(find(nancheck),:)=[]; 
                    stimCalcCut(find(nancheck),:)=[]; 
                end
                
                %if some are responsive then process the responsive trials
                %remove non-responsive calc events
                calcCalcCut_ttt = zeros(sum(responsiveInd), ...
                    size(calcCalcCut,2));
                locoCalcCut_ttt = zeros(sum(responsiveInd), ...
                    size(locoCalcCut,2));
                stimCalcCut_ttt = zeros(sum(responsiveInd), ...
                    size(stimCalcCut,2));
                diamCalcCut_ttt = zeros(sum(responsiveInd), ...
                    size(diamCalcCut,2), size(diamCalcCut,3));
                diamBSCalcCut_ttt = zeros(sum(responsiveInd), ...
                    size(diamBSCalcCut,2), size(diamBSCalcCut,3));
                for c = 1:size(responsiveInd,2) %loop the calcium events
                    %look for non responsive trials
                    %need to check if it is the first responsive trial, if
                    %so set counter to 1, and put that data first in matrix
                    %if it is a subsequent index, add 1 to the counter -
                    %then this is the place of the new data in the matrix
                    if responsiveInd(c) == 1 %check if responsive ind is 1
                        if c == find(responsiveInd==1,1) %make counter
                            counter = 1; 
                        else
                            counter = counter + 1; 
                        end
                        calcCalcCut_ttt(counter,:) = calcCalcCut(c,:);
                        locoCalcCut_ttt(counter,:) = locoCalcCut(c,:);
                        stimCalcCut_ttt(counter,:) = stimCalcCut(c,:);
                        diamCalcCut_ttt(counter,:,:) = diamCalcCut(c,:,:);
                        diamBSCalcCut_ttt(counter,:,:) = diamBSCalcCut(c,:,:);
                    end %end of check if responsive ind is 1
                end %end of looping calc events
                clear calcCalcCut locoCalcCut stimCalcCut diamCalcCut ...
                    diamBSCalcCut;
                calcCalcCut = calcCalcCut_ttt;
                locoCalcCut = locoCalcCut_ttt;
                stimCalcCut = stimCalcCut_ttt;
                diamCalcCut = diamCalcCut_ttt;
                diamBSCalcCut = diamBSCalcCut_ttt;
                clear calcCalcCut_ttt locoCalcCut_ttt stimCalcCut_ttt ...
                    diamCalcCut_ttt diamBSCalcCut_ttt;
                
                %check if any calc epochs have been detected
                if size(calcCalcCut,1)==1
                    %ie may not have any that fit now added baseline and 
                    %post-calc period
                    %inform user if no calc events, and no need to cut rest 
                    %of data
                    disp('no calc events');
                    disp('deleting any old mat files'); 
                    
                    %need to delete any old version of the files if saved
                    %into dir
                    matfile = fullfile([expDir,filesep,...
                        'LONGER_dataCutByCalc_xyFWHM_', ...
                        num2str(prefs.std_multiplier),'std.mat']);
                    delete(matfile);
                    
                else %cut rest of data if calc events detected
                    
                    %error bars for calcium traces (SEM)
                    [calcCalcCut_SEM]=getSEM(calcCalcCut,1);
                    %error bars for loco traces (SEM)
                    [locoCalcCut_SEM]=getSEM(locoCalcCut,1);
                    
                    %mean across branches for diam, and smooth
                    diamCalcCut=squeeze(nanmean(diamCalcCut,2));
                    %error bars for diam
                    [diamCalcCut_SEM]=getSEM(diamCalcCut,1);
                    
%                     for g = 1:size(diamBSCalcCut,2)
%                         %mean across branches for diam bootstrap
%                         diamBSCalcCut_ttt(:,g,:)=squeeze(nanmean( ...
%                             diamBSCalcCut(:,g,:,:),3));
%                     end
%                     clear diamBSCalcCut;
%                     diamBSCalcCut = diamBSCalcCut_ttt; 
%                     clear diamBSCalcCut_ttt;
                    %error bars for diam BS
                    for g = 1:size(diamBSCalcCut,2)
                        [diamBSCalcCut_SEM(g,:)]=getSEM(squeeze ...
                            (diamBSCalcCut(:,g,:)),1);
                    end

                    %save out the data stored in cells as all different 
                    %sizes, but not averaged together yet, just in case 
                    %need raw traces 
                    matfile = fullfile([expDir,filesep,...
                        'LONGER_dataCutByCalc_xyFWHM_', ...
                        num2str(prefs.std_multiplier),'std.mat']);
                    
                    save(matfile,'diamCalcCut','diamBSCalcCut','calcCalcCut', ...
                        'stimCalcCut','locoCalcCut','diamBSCalcCut_SEM',...
                        'diamCalcCut_SEM','calcCalcCut_SEM','calcEvents',...
                        'timeCalcCut','fps','responsiveInd','prefs1',...
                        '-v7.3');
                    
                    %binarise stim and loco - 1=stim/loco on, 0=off
                    for c = 1:size(locoCalcCut,1) %loop calc events
                        %check for 1s before the calc peak until calc peak
                        %if any numbers above 0 then class as loco
                        %According to Patrick Drew, even a flicker affects
                        %haem responses, so doesn't matter what size loco
                        %is in this period 
                        ttt = sum(locoCalcCut(c,find(timeCalcCut>=4,1): ...
                            find(timeCalcCut>=5,1))); 
                        if ttt > 0 %check if loco sum is above 0
                            calcSort.loco(c) = 1;
                        else
                            calcSort.loco(c) = 0;
                        end
                        clear ttt; 
                        ttt = sum(stimCalcCut(c,find(timeCalcCut>=4,1): ...
                            find(timeCalcCut>=5,1))); 
                        if ttt > 0 %check if stim sum is above 0
                            calcSort.stim(c) = 1;
                        else
                            calcSort.stim(c) = 0;
                        end
                    end %end of looping calc events

                    %%%% WARNING LAZY CODING
                    %CBA EDITING ALL SUBSEQUENT CODES YET, SO JUST RENAMING
                    %VARS INTO WHAT I CALL LATER... WILL NEED TO TIDY AT
                    %SOME POINT BEFORE SEND TO OTHERS
                    calcSort.time = timeCalcCut;
                    calcSort.calc = calcCalcCut; 
                    calcSort.diam = diamCalcCut;
                    calcSort.diamBS = diamBSCalcCut;
                    calcSort.loco_all = locoCalcCut;
                    calcSort.stim_all = stimCalcCut; 
                    calcSort.calc_SEM = calcCalcCut_SEM; 
                    calcSort.loco_SEM = locoCalcCut_SEM; 
                    calcSort.diam_SEM = diamCalcCut_SEM; 
                    calcSort.diamBS_SEM = diamBSCalcCut_SEM; 

                    %save data:
                    matfile = fullfile([expDir,filesep,'LONGER_calcSort_xyFWHM_'...
                        num2str(prefs.std_multiplier),'std.mat']);
                    save(matfile,'calcSort','fps','prefs1','responsiveInd');
                    
                    %% plots
                    
                    %average traces:
                    figure;
                    %make fig size of screen
                    set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                        screenSz(3)/2 screenSz(4)]);
                    %calc
                    subplot(2,1,1);
                    %plot raw Calc trace
                    plot(calcSort.time,nanmean(calcSort.calc,1),'g');
                    hold on;
                    if size(calcSort.calc,1)>2
                        shadedErrorBar(calcSort.time,nanmean( ...
                            calcSort.calc,1),calcSort. ...
                            calc_SEM,'lineProps','g');
                    end
                    %plot base
                    plot([calcSort.time(1) calcSort.time(end)],[0 0],'k');
                    %plot Calc peak, and norm period
                    clear yLimits; yLimits = get(gca, 'YLim');
                    plot([prefs.post(1) prefs.post(1)],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %calc peak
                    plot([2 2],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                    plot([3 3],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                    xlim([calcSort.time(1),calcSort.time(end)]);
                    title(['Calc Onset, nCalcEvents=',num2str(size ...
                        (calcSort.calc,1))]);
                    xlabel('Time(s)'); ylabel('\Delta F/F'); 
                    
                    %main title
                    mtit(gcf, 'Cut Vessel Diam by Calc Events',...
                        'fontsize',14,'color', [1 0 0],'xoff',0,'yoff',0.1);
                    
                    %diam:
                    subplot(2,1,2);
                    %plot diam
                    plot(calcSort.time,nanmean(calcSort.diam,1),'r');
                    hold on;
                    if size(calcSort.calc,1)>2
                        shadedErrorBar(calcSort.time, nanmean( ...
                            calcSort.diam,1), calcSort.diam_SEM, ...
                             'lineProps','r');
                    end
                    %plot base
                    plot([calcSort.time(1) calcSort.time(end)],[0 0],'k');
                    %plot calc onset line
                    clear yLimits; yLimits = get(gca, 'YLim');
                    plot([prefs.post(1) prefs.post(1)],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %calc peak
                    plot([2 2],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                    plot([3 3],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                    xlim([calcSort.time(1),calcSort.time(end)]);
                    title('Diameter');
                    xlabel('Time(s)'); ylabel('\Delta d/d'); 
                    
                    saveas(gcf, fullfile([expDir,filesep,...
                        'LONGER_cutByCalcPlot_', num2str(prefs.std_multiplier),...
                        'std.png']));
                    close;
                    
                    %individual calc increase and diam increase
                    figure;
                    %make fig size of screen
                    set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                        screenSz(3) screenSz(4)]);
                    
                    for c = 1:size(calcSort.calc,1)
                        clear calc_ttt diam_ttt; 
                        [calc_ttt]=calcSort.calc(c,:);
                        [diam_ttt]=calcSort.diam(c,:);
                        %plot
                        subplot(ceil(size(calcSort.calc,1)/4),4,c)
                        plot(calcSort.time, calc_ttt,'g');
                        hold on;
                        plot(calcSort.time, diam_ttt,'r');
                        plot([calcSort.time(1) calcSort.time(end)], ...
                            [0 0],'k');
                        yLimits=get(gca, 'YLim');
                        plot([prefs.post(1) prefs.post(1)],...
                        [yLimits(1) yLimits(2)],'b--','LineWidth',2); %calc peak loc
                        plot([-2 -2],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                        plot([-3 -3],...
                        [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                        xlim([calcSort.time(1) calcSort.time(end)]);
                        xlabel('time(s)'); 
                       ylabel('\DeltaF/F or \DeltaD/D');                      
                        if c==1
                           legend('calcium','vessel','AutoUpdate','off'); 
                        end
                    end
                    saveas(gcf, fullfile([expDir,filesep, ...
                        'LONGER_calcDiamTracesOverlaid_',num2str( ...
                        prefs.std_multiplier),'std.png']));
                    close;
                    
                    %also display as two imagescales, so can see if calc
                    %peak is earlier than diam, and if there is even a diam
                    %peak
                    %find time to subtract to align peaks to zero 
                    ttt = find(calcSort.time<=5);
                    timSubtract = calcSort.time(ttt(end)); 
                    %imagesc plot:
                    figure;
                    %make fig size of screen
                    set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                        screenSz(3)/2 screenSz(4)]);
                    subplot(211);
                    imagesc(calcSort.time - timSubtract, [], ...
                        calcSort.calc);
                    xlabel('time(s)');ylabel('calcEvent#'); 
                    title('calcium, peaks aligned to 0'); 
                    subplot(212);
                    imagesc(calcSort.time - timSubtract, [], ...
                        calcSort.diam);
                    hold on; 
                    plot([0 0],[1 size(calcSort.calc,1)],'k','LineWidth',2);
                    xlabel('time(s)');ylabel('calcEvent#'); 
                    title('corresponding diam'); 
                    saveas(gcf, fullfile([expDir,filesep, ...
                        'LONGER_calcDiamImagesc_',num2str( ...
                        prefs.std_multiplier),'std.png']));
                    close;
                    
                end %check if any calc epochs have been detected
            else %check if any responsive trials
                
                disp('no responsive calc trials');
                
            end %end of check for responsive trials
        end %check if any calcium events
        
    end %check if contData_xyFWHM.mat file exists
    
end %end of looping expDirs

end %end of function

