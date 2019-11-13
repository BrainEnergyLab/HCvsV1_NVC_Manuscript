
function cutByCMRO2Events_probe(fname, prefs)

%function by Kira, May 2018, updated May 2019
%written to detect all CMRO2 events deemed as 'responsive' (>2std above
%pre-CMRO2 peak baseline), and then also cut the haem during this cmro2
%increase - to see if cmro2 increase causes any vascular/haem changes

%INPUTS-
%fname = top directory which contains calcAutoRoi.mat (net calc) file and
%contData.mat (cont probe data traces) file, NB/ the data should
%be separated by group (e.g. brain region, or APOE genotype), and you
%should be within the group directory (e.g. V1)

%OUTPUTS-
%saves figures to show average calc peaks and corresponding diam, and saves
%mat file with data in inside individual expDirs - run later code to
%compare across groups

%Useful functions needed for code to run:
%findFolders, findLocoEvents, sortDataByResponsiveTrials

if nargin < 2
    %set which haem param you want to cut around
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
    %for responsive trial sorter-
    %time in seconds to average either side of the max peak, i.e. to see if
    %the peak is above the std dev of the baseline - this may need to be
    %larger for slower responses (e.g. haem vs calc)
    prefs.avgPeak = 0.1; %seconds
end

%find all the data which has been run through net calc extractor
[findmatfile]=findFolders(fname,'contData.mat');
%smooth the data? 0 = no, 1=yes
smoothVessFlag = 0;

%get screen size for plots
screenSz=get(0,'Screensize');

for n = 1:size(findmatfile,2)
    
    clearvars -except fname prefs findmatfile screenSz n smoothVessFlag
    
    %find each session experimental directory containing net calc data
    [expDir,~]=fileparts(findmatfile{1,n});
    
    %inform user of progress
    disp(['processing file ',num2str(n),'/', ...
        num2str(size(findmatfile,2))]);
    disp(expDir);
    
    %load the data
    disp('loading data...');
    load([expDir, filesep, 'contData.mat']);
    
    %extract calcium from data channel 7
    %NB need to smooth cmro2 data or cannot detect peaks
    calcium = smooth(data(7,:),0.01,'loess')';
    %take all haem parameters
    haem = data([1:6],:); 
    
    %check if the experiment is a stimulation experiment
    %NB stim needs to be in the name of the exp dir folder
    stimExpFlag=regexpi(expDir,'stim');
    
    %convert some of the preferences to frames, using the now loaded
    %fps
    prefs1.minDist = round(prefs.minDist * fps);
    prefs1.base = round(prefs.base(2) * fps);
    prefs1.post = round(prefs.post(2) * fps);
    prefs1.plotFlag = prefs.plotFlag;
    prefs1.flickerFlag = prefs.flickerFlag;
    
    %% Bootstrap linear model
    %%%%%%%%% bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %also create multiple reshuffled traces of the vasc diam
    %reshuffle the haem data - keep the calc as is, so have same calc
    %spikes detected
    nIts        = 100;  % Number of repeats of randomization
    %Initialize  variable to save traces in
    shuffledHaemTrace = nan([size(haem,1),nIts,size(haem,2)]);
    
    %shuffle data over specified iterations
    for n = 1:size(haem,1)
        for i = 1:nIts
            haemTraceTemp = haem(n,:);
            %Shuffle the rows (using randperm), and then stitch back together
            shuffledHaemTrace(n,i,:) = ...
                reshape(haemTraceTemp(:,randperm(size(haemTraceTemp,2))),[],1)';
        end % Number of iterations
    end
    
    %% Detect peaks in CMRO2 signal
    
    %normalise calcium trace - remove drift etc. - so can detect peaks
    calcium_norm = calcium;
    calcium_norm = calcium_norm - nanmin(calcium_norm);
    calcium_norm = calcium_norm / nanmax(calcium_norm);
    calcium_norm = calcium_norm - smooth(calcium_norm,100)';
    
    %look for peaks in calc trace from normalised trace (otherwise
    %trend drifts will effect peak detection)
    [calcEvents] = findSignalSpikesInTrace(calcium_norm, fps);
    
    %plot RAW trace with detected spikes
    figure;
    
    %make fig 1/2 size of screen
    set(gcf, 'Position', [screenSz(1) screenSz(2) ...
        screenSz(3)/2 screenSz(4)]);
    
    plot(time, calcium, 'g');
    hold on;
    %       plot(time(:,calcEvents(1,:)),calcium(:,calcEvents(1,:)), 'ko'); %start of peak
    plot(time(:,calcEvents(2,:)),calcium(:,calcEvents(2,:)), 'bo'); %peak
    title('calcium trace with peaks');
    legend('calc','peak');
    xlabel('time(s)'); %ylabel('\Delta/F/F');
    saveas(gcf, fullfile([expDir,filesep,'detectedCMRO2peaks_all.png']));
    close;
    
    %check if any CMRO2 events were actually detected
    if isempty(calcEvents) == 1
        
        %inform user that there are no loco events in this specific exp
        disp(['no calc events: ', expDir]);
        
    else %if there are CMRO2 events detected, proceed...
        
        %check there are no zeros in calcEvents, if so remove
        zeroIndSearch=find(calcEvents(1,:)==0);
        if ~isempty(zeroIndSearch)
            calcEvents(:,zeroIndSearch)=[];
        end
        
        %% plot detected calc events again, but show all haem cont traces 
        figure;
        %make fig size of screen
        set(gcf, 'Position', [screenSz(1) screenSz(2) ...
            screenSz(3) screenSz(4)]);
        
        ax1 = subplot(611);
        plot(calcium, 'g');
        hold on;
        %plot calc events detected on calc trace
        %find the value of the calc data at the time point of the calc
        %event
        plot(calcEvents(2,:), calcium(1,calcEvents(2,:)), ...
            'ko', 'LineWidth',2);
        title('CMRO2, with all detected events');
        xlabel('frames');
        xlim([0 size(calcium,2)]);
        
        ax2 = subplot(612);
        plot(haem(1,:),'m');
        hold on;
        %plot calc events detected on haem trace
        plot(calcEvents(2,:), haem(1,calcEvents(2,:)),'ko', 'LineWidth',2);
        title('flux');
        xlabel('frames'); ylabel('AU');
        xlim([0 size(calcium,2)]);
        
        ax3 = subplot(613);
        plot(haem(2,:),'k');
        hold on;
        %plot calc events detected on haem trace
        plot(calcEvents(2,:), haem(2,calcEvents(2,:)),'ro', 'LineWidth',2);
        title('speed');
        xlabel('frames'); ylabel('AU');
        xlim([0 size(calcium,2)]);
        
        ax4 = subplot(614);
        plot(haem(3,:),'k');
        hold on;
        %plot calc events detected on haem trace
        plot(calcEvents(2,:), haem(3,calcEvents(2,:)),'ro', 'LineWidth',2);
        title('so2');
        xlabel('frames'); ylabel('%');
        xlim([0 size(calcium,2)]);
        
        ax5 = subplot(615);
        plot(haem(4,:),'r');
        hold on;
        plot(haem(5,:),'b');
        plot(haem(6,:),'g');
        %plot calc events detected on haem trace
        plot(calcEvents(2,:), haem(6,calcEvents(2,:)),'ko', 'LineWidth',2);
        title('haem');
        xlabel('frames'); ylabel('AU');
        xlim([0 size(calcium,2)]);
        
        ax6 = subplot(616);
        plot(locomotion,'b');
        xlabel('frames'); ylabel('AU');
        xlim([0 size(calcium,2)]);
        
        linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x');
        
        saveas(gcf, fullfile([expDir,filesep, ...
            'contHaemPlots_detectedCMRO2Events.png']));
        close;
        
        %% check for CMRO2 events having enough time before and after
        %i.e. for cutting around them
        enoughTimInd = 0; %preallocate the index with a zero
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
        
        %% cut around CMRO2 peaks - i.e. so peaks are all aligned
        %also cut around corresponding haem and haem BS data
        %as well as loco and stim
        for c = 1:size(calcEvents,2)
            %calcium trace
            calcCalcCut(c,:) = calcium(:,calcEvents(2,c)-round( ...
                prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                prefs.post(2)*fps)- (prefs.post(1)*fps))+1 );
            %diameter trace
            haemCalcCut(c,:,:) = haem(:,calcEvents(2,c)-round( ...
                prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                prefs.post(2)*fps)- (prefs.post(1)*fps))+1 );
            %bootstrap diam
            for n = 1:size(shuffledHaemTrace,1)
                haemBSCalcCut(c,n,:,:) = squeeze(shuffledHaemTrace(n,:,calcEvents(2,c)-round( ...
                    prefs.base(2)*fps)+1 : calcEvents(2,c)+round(( ...
                    prefs.post(2)*fps)- (prefs.post(1)*fps))+1));
            end
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
        %create corresponding time vector
        timeCalcCut = [1:size(calcCalcCut,2)]/fps;
        
        %% normalise data for first 1s -i.e. way before peak of calc:
        %CMRO2
        for c = 1:size(calcCalcCut,1)
            [calcCalcCut(c,:)]=findDeltaF(calcCalcCut(c,:), ...
                [round(2*fps) round(3*fps)]);
        end
        %haem
        for c = 1:size(haemCalcCut,1)
            [haemCalcCut(c,:,:)]=findDeltaF(haemCalcCut(c,:,:), ...
                [round(2*fps) round(3*fps)]);
        end
        %haem bootstrap
        for c = 1:size(haemBSCalcCut,1)
            for n = 1:size(haemBSCalcCut,2)
                [haemBSCalcCut(c,n,:,:)]=findDeltaF(squeeze(haemBSCalcCut(c,n,:,:)), ...
                    [round(2*fps) round(3*fps)]);
            end
        end
        
        %% get most responsive index
        %NB do this on the normalised trace,
        %as sloping trend could mean some calc events are deemed not
        %sig due to decreasing baseline
        
        %call function:
        clear responsiveInd;
        for c = 1:size(calcCalcCut,1)
            [responsiveInd(c)] = sortDataByResponsiveTrials ...
                (calcCalcCut(c,:),time,prefs,'calc');
        end
        
        %% plot the trials which are categorised as responsive and non
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
            'categoriseResponsiveCMRO2Events_', ...
            num2str(prefs.std_multiplier),'std.png']));
        close;
        
        %% check if any of the trials have been classed as responsive
        if nansum(responsiveInd) > 0
            
            %redo the calcium and diameter cont traces
            %Now with only the responsive events marked:
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
            title('calcium, responsive peaks only');
            xlabel('frames'); xlim([0 size(calcium,2)]);
            ylabel('AU');
            
            ax2=subplot(212);
            plot(haem(1,:),'m');
            hold on;
            %plot calc events detected on diam trace
            plot(calcEvents(2,responsiveInd==1), haem(1, ...
                calcEvents(2,responsiveInd==1)),'ko', 'LineWidth',2);
            title('flux, for EG');
            xlabel('frames'); xlim([0 size(calcium,2)]);
            ylabel('AU');
            
            linkaxes([ax1,ax2],'x');
            saveas(gcf, fullfile([expDir,filesep, ...
                'contPlots_CMRO2Events_fluxEG_responsiveOnly.png']));
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
            
            %% Process the responsive CMRO2 trials only
            %remove non-responsive CMRO2 events - i.e. dnt care about
            %them if non-resp (as looking at haem resp to sig CMRO2
            %events)
            calcCalcCut_ttt = zeros(sum(responsiveInd), ...
                size(calcCalcCut,2));
            locoCalcCut_ttt = zeros(sum(responsiveInd), ...
                size(locoCalcCut,2));
            stimCalcCut_ttt = zeros(sum(responsiveInd), ...
                size(stimCalcCut,2));
            haemCalcCut_ttt = zeros(sum(responsiveInd), ...
                size(haemCalcCut,2), size(haemCalcCut,3));
            haemBSCalcCut_ttt = zeros(sum(responsiveInd), ...
                size(haemBSCalcCut,2), size(haemBSCalcCut,3), ...
                size(haemBSCalcCut,4));
            for c = 1:size(responsiveInd,2) %loop the CMRO2 events
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
                    haemCalcCut_ttt(counter,:,:) = haemCalcCut(c,:,:);
                    haemBSCalcCut_ttt(counter,:,:,:) = haemBSCalcCut(c,:,:,:);
                end %end of check if responsive ind is 1
            end %end of looping calc events
            clear calcCalcCut locoCalcCut stimCalcCut haemCalcCut ...
                haemBSCalcCut;
            calcCalcCut = calcCalcCut_ttt;
            locoCalcCut = locoCalcCut_ttt;
            stimCalcCut = stimCalcCut_ttt;
            haemCalcCut = haemCalcCut_ttt;
            haemBSCalcCut = haemBSCalcCut_ttt;
            clear calcCalcCut_ttt locoCalcCut_ttt stimCalcCut_ttt ...
                haemCalcCut_ttt haemBSCalcCut_ttt;
            
            %% check if any CMRO2 epochs have been detected now
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
                    'dataCutByCMRO2peaks.mat']);
                delete(matfile);
                
            else %cut rest of data if CMRO2 events detected
                
                %error bars for calcium traces (SEM)
                [calcCalcCut_SEM]=getSEM(calcCalcCut,1);
                %error bars for loco traces (SEM)
                [locoCalcCut_SEM]=getSEM(locoCalcCut,1);
                %error bars for haem
                for i = 1:size(haemCalcCut,2)
                    [haemCalcCut_SEM(i,:)]=getSEM(squeeze( ...
                        haemCalcCut(:,i,:)),1);
                end
                %errorbars for bootstrap haem
                for i = 1:size(haemBSCalcCut,2)
                    for j = 1:size(haemBSCalcCut,3)
                        [haemBSCalcCut_SEM(i,j,:)]=getSEM(squeeze ...
                            (haemBSCalcCut(:,i,j,:)),1);
                    end
                end
                
%                 %save out the data stored in cells as all different
%                 %sizes, but not averaged together yet, just in case
%                 %need raw traces
%                 matfile = fullfile([expDir,filesep,...
%                     'dataCutByCMRO2peaks_indCells.mat']);
%                 
%                 save(matfile,'haemCalcCut','haemBSCalcCut','calcCalcCut', ...
%                     'stimCalcCut','locoCalcCut','haemBSCalcCut_SEM',...
%                     'haemCalcCut_SEM','calcCalcCut_SEM','calcEvents',...
%                     'timeCalcCut','fps','responsiveInd','prefs1',...
%                     '-v7.3');
                
                %% binarise stim and loco - 1=stim/loco on, 0=off
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
                calcSort.haem = haemCalcCut;
                calcSort.haemBS = haemBSCalcCut;
                calcSort.loco_all = locoCalcCut;
                calcSort.stim_all = stimCalcCut;
                calcSort.calc_SEM = calcCalcCut_SEM;
                calcSort.loco_SEM = locoCalcCut_SEM;
                calcSort.haem_SEM = haemCalcCut_SEM;
                calcSort.haemBS_SEM = haemBSCalcCut_SEM;
                
                %save data:
                matfile = fullfile([expDir,filesep,'dataCutByCMRO2peaks']);
                save(matfile,'calcSort','fps','prefs1');
                
                %% plots
                
                %average traces:
                figure;
                %make fig size of screen
                set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                    screenSz(3)/2 screenSz(4)]);
                %calc
                subplot(4,1,1);
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
                mtit(gcf, 'Cut Haem by Calc Events',...
                    'fontsize',14,'color', [1 0 0],'xoff',0,'yoff',0.1);
                
                %so2:
                subplot(4,1,2);
                %plot haem
                plot(calcSort.time,squeeze(nanmean(calcSort.haem(:,3,:),1)),'k');
                hold on;
                if size(calcSort.calc,1)>2
                    shadedErrorBar(calcSort.time, squeeze(nanmean( ...
                        calcSort.haem(:,3,:),1)), calcSort.haem_SEM(3,:), ...
                        'lineProps','k');
                end
                %plot base
                plot([calcSort.time(1) calcSort.time(end)],[0 0],'k');
                %plot calc onset line
                clear yLimits; yLimits = get(gca, 'YLim');
                plot([prefs.post(1) prefs.post(1)],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %cmro2 peak
                plot([2 2],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                plot([3 3],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                xlim([calcSort.time(1),calcSort.time(end)]);
                title('so2');
                xlabel('Time(s)'); ylabel('\Delta d/d');
                
                %flux:
                subplot(4,1,3);
                %plot haem
                plot(calcSort.time,squeeze(nanmean(calcSort.haem(:,1,:),1)),'m');
                hold on;
                if size(calcSort.calc,1)>2
                    shadedErrorBar(calcSort.time, squeeze(nanmean( ...
                        calcSort.haem(:,1,:),1)), calcSort.haem_SEM(1,:), ...
                        'lineProps','m');
                end
                %plot base
                plot([calcSort.time(1) calcSort.time(end)],[0 0],'k');
                %plot calc onset line
                clear yLimits; yLimits = get(gca, 'YLim');
                plot([prefs.post(1) prefs.post(1)],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %cmro2 peak
                plot([2 2],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                plot([3 3],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                xlim([calcSort.time(1),calcSort.time(end)]);
                title('Flux');
                xlabel('Time(s)'); ylabel('\Delta d/d');
                
                %flux:
                subplot(4,1,4);
                %plot haem
                plot(calcSort.time,squeeze(nanmean(calcSort.haem(:,6,:),1)),'r');
                hold on;
                if size(calcSort.calc,1)>2
                    shadedErrorBar(calcSort.time, squeeze(nanmean( ...
                        calcSort.haem(:,6,:),1)), calcSort.haem_SEM(6,:), ...
                        'lineProps','r');
                end
                %plot base
                plot([calcSort.time(1) calcSort.time(end)],[0 0],'k');
                %plot calc onset line
                clear yLimits; yLimits = get(gca, 'YLim');
                plot([prefs.post(1) prefs.post(1)],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %cmro2 peak
                plot([2 2],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                plot([3 3],...
                    [yLimits(1) yLimits(2)],'k--','LineWidth',2); %norm period
                xlim([calcSort.time(1),calcSort.time(end)]);
                title('Hbt');
                xlabel('Time(s)'); ylabel('\Delta d/d');
                
                saveas(gcf, fullfile([expDir,filesep,'cutByCMRO2Plot_haems.png']));
                close;
                
                %individual calc increase and haem increase
                figure;
                %make fig size of screen
                set(gcf, 'Position', [screenSz(1) screenSz(2) ...
                    screenSz(3) screenSz(4)]);
                
                for c = 1:size(calcSort.calc,1)
                    clear calc_ttt haem_ttt;
                    [calc_ttt]=calcSort.calc(c,:);
                    [haem_ttt]=squeeze(calcSort.haem(c,:,:));
                    %plot
                    subplot(ceil(size(calcSort.calc,1)/4),4,c)
                    plot(calcSort.time, calc_ttt,'g');
                    hold on;
                    plot(calcSort.time, haem_ttt(1,:),'m'); %flux
                    plot(calcSort.time, haem_ttt(3,:),'k'); %so2
                    plot(calcSort.time, haem_ttt(6,:),'r'); %hbt
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
                        legend('cmro2','flux','so2','Hbt','AutoUpdate','off');
                    end
                end
                saveas(gcf, fullfile([expDir,filesep, ...
                    'CMRO2HaemTracesOverlaid.png']));
                close;

                
            end %check if any cmro2 epochs have been detected
        else %check if any responsive trials
            
            disp('no responsive calc trials');
            
        end %end of check for responsive trials
    end %check if any calcium events
    
    
end %end of looping expDirs

end %end of function

