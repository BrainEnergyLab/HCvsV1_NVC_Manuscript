
function cutByLScalc(topDir, prefs)

%must have run extractLScalcium and linescanVelocityAnalysis first - i.e.
%so have the mat files 'contData_ls_calc.mat' and 'contData_ls_RBCV.mat'
%in the experimental folders

if nargin < 2
    %threshold for removing calcium signal below (just for detecting
    %spikes)
    prefs.calcThresh = 1.5; %1.5 %std
    %minimum time between calcium peaks, will check each peak, and where
    %multiple are within min time of each other, will select the largest
    %peak - this number should at least exceed the norm period before calc
    %peak, i.e. so a responsive pt is not being included in the 'baseline'
    %period
    prefs.timeBtwnPks = 2; %3 %seconds
    %minimum number of detected calcium peaks required per vessel for it to
    %be included (i.e. if only has 1 or 2 may not be worth it)
    %this is done after all the catches, so may have some peaks removed for
    %not having enough time etc. 
    prefs.minCalcEvents = 3; 
    %time (in secs) before and after calc peak to take RBCV trace and for
    %averaging
    prefs.timeAvgRBCV = 3; %seconds
    %minimum amount of time within which the RBCV data in the period
    %before/after calc pk needs to contain data (i.e. in case too much of
    %it is filled with NaNs) - e.g. if there where only 3 data pts out of
    %100 which had actual data and not nans, the nanmean would still give
    %you a value, but this would have been sample from like 0.1 seconds, so
    %isn't really a fair value
    prefs.minTime4Data = 1.5; %seconds
    % amount of time before calc peak to normalise baseline during, so must
    % be no larger than amount of time in baseline (i.e. prefs.timeAvgRBCV)
    prefs.normb4calc = 1; %seconds
%     %amount of time to cut before and after calc peak for looking at RBCV
%     %traces (and getting max pk etc.) 
%     prefs.cutRBCVtrace = [1 3]; %seconds [b4calc aftercalc]
    %amount of time before and after peak dilation to average
    prefs.timeAroundPk = 0.2; %0.2 %seconds
    %threshold for determining if RBCV trace is responsive vs it's own
    %normalised baseline (i.e. whether peak is >= pref * std(baseline)) 
    prefs.respThresh = 1.5; %1.5 %std 
end

%get screen sz for plots
screenSz = get(0,'Screensize');

%% find all the subfolders - i.e. groups - in the top dir
listing=dir(topDir);
for a = 1:size(listing,1)-2
    if size(listing(2+a).name,2) <= 5
        grpNames{a} = listing(2+a).name;
    else
        grpNames{a} = [];
    end
end
clear a;
%remove the empty cells
grpNames = grpNames(~cellfun('isempty',grpNames));

%loop the groups, extract the data
for n = 1:size(grpNames,2)
    
    clearvars -except n grpNames topDir grp screenSz prefs
    
    %inform user of progress
    disp(grpNames{n});
    %label group 
    grp{n}.label = grpNames{n}; 
    
    %check for existence of the extracted ls calc file
    find_tif_file = findFolders([topDir,filesep,grpNames{n}], 'contData_ls_calc.mat');
    
    %loop through all directories with tif file
    for a = 1:size(find_tif_file,2)
        
        %inform user of progress
        disp(['Looping experiments: ', num2str(a), '/', ...
            num2str(size(find_tif_file,2))]);
        
        %find the local folder containing the tif file, get exp dir
        [expDir,~] = fileparts(find_tif_file{1,a});
        %search for RBCV file within exp dir, if doesn't exist inform user
        matFile = findFolders(expDir, 'contData_ls_RBCV.mat');

        if ~isempty(matFile) %check if RBCV mat file exists
            
            %% load data:
            load(matFile{1,1});
            load(find_tif_file{1,a});
            
            %rough diam not working for a sig portion of vessels - will
            %exclude for now
%             %load rough diam of vessel
%             load([expDir,filesep,'roughdiam.mat']);

            
            %catch to check the calc and RBCV have been extracted same:
            if size(velocity,2) == size(calcium,2)
                
                % plots with all time series if want to look at them
                if false
                    figure;
                    subplot(311);
                    plot(time,velocity,'r');
                    title('RBCV');
                    subplot(312);
                    plot(time,calcium,'g');
                    title('calcium');
                    subplot(313);
                    plot(time,locomotion,'b');
                    title('locomotion');
                    xlabel('time(s)');
                end
                
                %% clean the data, for detecting spikes in signal:
                
                %smooth calc trace
                calcium = smooth(calcium); 
                if size(calcium,1) > size(calcium,2)
                    calcium = calcium';
                end
                
                %normalise calc trace
                calcium_norm = calcium - nanmin(calcium);
                calcium_norm = calcium_norm / nanmax(calcium_norm);
                calcium_norm = calcium_norm - smooth(calcium_norm,100)';
                %remove any negative values
                calcium_norm(:, find(calcium_norm < 0)) = 0;
                %find any values below thresh and set them to zero - so just
                %left with spikes
                calcium_norm(:, find(calcium_norm < ...
                    nanstd(calcium_norm)*prefs.calcThresh)) = 0;
                
                % smooth the RBCV trace
                velocity = smoothKnit(velocity); 
                if size(velocity,1) > size(velocity,2)
                    velocity = velocity';
                end
                
                %% find final calc pks by: minimum spacing & largest pk
                
                %find all the values which aren't zero, aka the calc peaks
                calcPkInd_all = find(calcium_norm);
                %get minimum distance (in frames) need between calcium peaks
                %i.e. based on prefs specified (in seconds)
                minDistFrames = find(time>=prefs.timeBtwnPks,1);
                
                %create matrix with indices for calc periods to search for
                %largest peak within - NB this number is how many pks you have
                %found
                newCalcPeriodInd = [1, find(diff(calcPkInd_all)>=minDistFrames)+1];
                %loop number of calc pks - to find largest pk per group
                calcPkInd = [];
                for b = 1:size(newCalcPeriodInd,2)
                    if b < size(newCalcPeriodInd,2) %check if last IDd pk
                        %get out each index within the same calc pk grp
                        checkPkSzInd = calcPkInd_all(newCalcPeriodInd(b): ...
                            newCalcPeriodInd(b+1)-1);
                    else %last pk
                        %get out each index within the same calc pk grp
                        checkPkSzInd = calcPkInd_all(newCalcPeriodInd(b): ...
                            size(calcPkInd_all,2));
                    end %end of check if last pk
                    %find largest calc pk within this grp
                    [~,calcInd] = max(calcium_norm(:,checkPkSzInd));
%                     %TEST: just take 1st calc peak in the group - as looks
%                     %like RBCV already changing too early otherwise??
%                     calcInd = 1;
                    calcPkInd(b) = checkPkSzInd(calcInd);
                    clear checkPkSzInd calcInd;
                end %end of loop number of calc pks
                clear b calcPkInd_all newCalcPeriodInd;
                
                % original data and filtered data, with detected calcium spikes
                % overlaid
                figure;
                %make fig size of screen
                set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
                    screenSz(4)]);
                subplot(311);
                plot(time,calcium, 'g');
                hold on;
                plot(time(:,calcPkInd),calcium(:,calcPkInd),'ko');
                title('calcium');
                subplot(312);
                plot(time,calcium_norm,'k');
                hold on;
                plot(time(:,calcPkInd),calcium_norm(:,calcPkInd),'ro');
                title('calcium filtered');
                subplot(313);
                plot(time,velocity,'r');
                hold on;
                plot(time(:,calcPkInd),velocity(:,calcPkInd),'ko');
                xlabel('time(s)');
                %save figure into expdir
                figSave = 'dataTraces_calcPksMarked.png';
                saveas(gcf, fullfile([expDir,filesep,figSave]));
                close(gcf);
                
                %% find the avg RBCV before and after calc peak:
                
                %convert the specified time to average RBCV within (before and
                %after calc peak) into frames
                frames2avg = find(time>=prefs.timeAvgRBCV,1);
                frames4data = find(time>=prefs.minTime4Data,1);
                
                %loop calcium peaks, and get average RBCV in time before and
                %after
                removePk = []; 
                for b = 1:size(calcPkInd,2)
                    %check if have enough time for RBCV averaging
                    if calcPkInd(b) >= frames2avg && ...
                            calcPkInd(b)+frames2avg <= size(velocity,2)
                        %create temp vars with extracted RBCV data
                        %as need to check that it isn't filled with too many
                        %NaNs, so just 1 sample pt affects data
                        velB4 = velocity(calcPkInd(b)-frames2avg:calcPkInd(b)-1);
                        velAfter = velocity(calcPkInd(b)+1:calcPkInd(b)+frames2avg);
                        %check not too many NaNs inside either data set
                        %if too many will need to exclude this pk
                        if size(find(isnan(velB4)),2) < frames4data && ...
                                size(find(isnan(velAfter)),2) < frames4data
                            %take actual data traces
                            grp{n}.vess{a}.RBCV(b,:) = velocity(:, ...
                                calcPkInd(b)-frames2avg: ...
                                calcPkInd(b)+frames2avg); 
                            grp{n}.vess{a}.calcium(b,:) = calcium(:, ...
                                calcPkInd(b)-frames2avg: ...
                                calcPkInd(b)+frames2avg);
                            % average the RBCV data in time before and after
                            % calc
                            %data has passed the checks
                            %get mean RBCV:
                            grp{n}.vess{a}.velB4_avg(b) = nanmean(velB4);
                            grp{n}.vess{a}.velAfter_avg(b) = nanmean(velAfter);
                            %get CV of RBCV
                            grp{n}.vess{a}.velB4_cv(b) = nanstd(velB4) ...
                                ./ nanmean(velB4);
                            grp{n}.vess{a}.velAfter_cv(b) = ...
                                nanstd(velAfter) ./ nanmean(velAfter);
                        else %too many nans
                            disp(['skipping calc pk ', num2str(b), ...
                                ', too many nans']);
                            %fill traces with NaNs
                            grp{n}.vess{a}.RBCV(b,:) = nan(1,size( ...
                                calcPkInd(b)-frames2avg : ...
                                calcPkInd(b)+frames2avg,2));
                            grp{n}.vess{a}.calcium(b,:) = nan(1,size( ...
                                calcPkInd(b)-frames2avg: ...
                                calcPkInd(b)+frames2avg,2));
                            %fill average RBCV val with NaNs
                            grp{n}.vess{a}.velB4_avg(b) = NaN;
                            grp{n}.vess{a}.velAfter_avg(b) = NaN;
                            grp{n}.vess{a}.velB4_cv(b) = NaN;
                            grp{n}.vess{a}.velAfter_cv(b) = NaN;
                            %save out index of which pks to remove
                            if isempty(removePk) %check if already have remove pk ind
                                removePk = b;
                            else %already have a remove pk in
                                removePk = [removePk; b]; %put next pk to remove on end
                            end %end of check remove pk ind
                        end %end of check if too many nans
                    else %calc pk not enough time before or after for avg
                        %inform user which calc pk is being skipped
                        disp(['skipping calc pk ', num2str(b), ...
                            ', not enough time']);
                        %fill traces with NaNs
                        grp{n}.vess{a}.RBCV(b,:) = nan(1,size( ...
                            calcPkInd(b)-frames2avg : ...
                            calcPkInd(b)+frames2avg,2));
                        grp{n}.vess{a}.calcium(b,:) = nan(1,size( ...
                            calcPkInd(b)-frames2avg: ...
                            calcPkInd(b)+frames2avg,2));
                        %fill average RBCV val with NaNs
                        grp{n}.vess{a}.velB4_avg(b) = NaN;
                        grp{n}.vess{a}.velAfter_avg(b) = NaN;
                        grp{n}.vess{a}.velB4_cv(b) = NaN;
                        grp{n}.vess{a}.velAfter_cv(b) = NaN;
                        %save out index of which pks to remove
                        if isempty(removePk) %check if already have remove pk ind
                            removePk = b;
                        else %already have a remove pk in
                            removePk = [removePk; b]; %put next pk to remove on end
                        end %end of check remove pk ind
                    end %end of check if have enough time for RBCV averaging
                    clear velB4 velAfter; %clear temp vars created for averaging
                end %end of loop calc pks
                
                %% remove calc pk inds which have no vel data:
                calcPkInd(:,removePk)=[]; 
                grp{n}.vess{a}.calcPkInd = calcPkInd; clear calcPkInd; 
                grp{n}.vess{a}.velB4_avg(:,removePk)=[]; 
                grp{n}.vess{a}.velAfter_avg(:,removePk)=[]; 
                grp{n}.vess{a}.velB4_cv(:,removePk)=[]; 
                grp{n}.vess{a}.velAfter_cv(:,removePk)=[]; 
                grp{n}.vess{a}.RBCV(removePk,:)=[]; 
                grp{n}.vess{a}.calcium(removePk,:)=[]; 
                
                %remove any vessels which dnt contribute at least 5
                %responses
                if size(grp{n}.vess{a}.calcPkInd,2)>=prefs.minCalcEvents
                    
                    %label exp dir and animal for vessels which make it
                    %through all checks
                    grp{n}.vess{a}.anLabel = extractBefore(extractAfter ...
                        (expDir,[topDir,filesep,grpNames{n},'\']),'\');
                    grp{n}.vess{a}.expDir = extractAfter(expDir, ...
                        [topDir,filesep,grpNames{n},filesep,...
                        grp{n}.vess{a}.anLabel,'\']);
                    %             %store vess diam info
                    %             grp{n}.vess{a}.roughdiam = roughdiam;
                    
                    %create time vector
                    grp{n}.vess{a}.time = time(:,1:size(grp{n}.vess{a}.RBCV,2));
                    
                    %% normalise avg RBCV val after to avg RBCV val before
                    %(i.e. so can see how it changes)
                    grp{n}.vess{a}.velNorm = grp{n}.vess{a}.velB4_avg ./ ...
                        grp{n}.vess{a}.velAfter_avg;
                    
                    %% work with the traces - normalise RBCV trace and get pk val
                    %then can plot to see if any diam changes
                    %for averaging pk dilation (i.e. so dnt just take a random
                    %noise spike) - will take time either side of max pk and
                    %average this
                    framesAroundPk = find(time>=prefs.timeAroundPk,1);
                    for b = 1:size(grp{n}.vess{a}.RBCV,1) %loop calc pk RBCV traces
                        %normalise the data to 1s before calc peak
                        [grp{n}.vess{a}.RBCV_norm(b,:)] = findDeltaF( ...
                            grp{n}.vess{a}.RBCV(b,:), ...
                            [(frames2avg-find(time>=prefs.normb4calc,1)) frames2avg]);
                        %smoothing whole RBCV trace at start instead of in
                        %individual calc pk traces
%                         %smooth the data as RBCV is very noisy
%                         grp{n}.vess{a}.RBCV_norm(b,:) = smoothKnit(grp{n}.vess{a}.RBCV_norm(b,:)); 
                        %get the avg peak val - but only search during post
                        %calc time
                        [~,maxInd] = max(grp{n}.vess{a}.RBCV_norm ...
                            (b,frames2avg+1:end));
                        %need to find the index of the peak across the whole
                        %trace
                        pkInd = (frames2avg+1) + maxInd; clear maxInd;
                        %check averaging data for pk can fit
                        if pkInd + framesAroundPk < size(grp{n}.vess{a}.RBCV_norm,2) && ...
                                pkInd - framesAroundPk > 0
                            grp{n}.vess{a}.RBCVpk(b) = nanmean( ...
                                grp{n}.vess{a}.RBCV_norm ...
                                (b,pkInd-framesAroundPk:pkInd+framesAroundPk),2);
                            %deem whether RBCV peak is responsive or not resp
                            if grp{n}.vess{a}.RBCVpk(b) >= ...
                                    prefs.respThresh * nanstd( ...
                                    grp{n}.vess{a}.RBCV_norm(b,1:frames2avg))
                                %responsive
                                grp{n}.vess{a}.respInd(b) = 1;
                            else %not responsive
                                grp{n}.vess{a}.respInd(b) = 0;
                            end
                        else %cannot average for max pk, not enough time
                            grp{n}.vess{a}.RBCVpk(b) = NaN;
                            grp{n}.vess{a}.respInd(b) = NaN;
                        end %end of check if enough time for averaging data pk
                    end %end of loop calc pk RBCV traces
                    
                    %remove any peaks which didn't have enough time for
                    %averaging peak so couldnt get peak val
                    removePk = find(isnan(grp{n}.vess{a}.RBCVpk));
                    grp{n}.vess{a}.calcPkInd(removePk) = [];
                    grp{n}.vess{a}.RBCVpk(removePk) = [];
                    grp{n}.vess{a}.respInd(removePk) = [];
                    grp{n}.vess{a}.velB4_avg(:,removePk)=[];
                    grp{n}.vess{a}.velAfter_avg(:,removePk)=[];
                    grp{n}.vess{a}.velB4_cv(:,removePk)=[];
                    grp{n}.vess{a}.velAfter_cv(:,removePk)=[];
                    grp{n}.vess{a}.velNorm(:,removePk)=[];
                    grp{n}.vess{a}.RBCV(removePk,:)=[];
                    grp{n}.vess{a}.RBCV_norm(removePk,:)=[];
                    grp{n}.vess{a}.calcium(removePk,:)=[];
                    
                    %% get out average responsiveness % per vessel
                    grp{n}.vess{a}.respPcnt = ...
                        (size(find(grp{n}.vess{a}.respInd),2) ...
                        / size(grp{n}.vess{a}.respInd,2))*100;
                    
                    %% plots
                    
                    %plot of average RBCV trace before and after calc
                    %calc marked with black line (at prefs.timeAvgRBCV)
                    figure;
                    %make fig 1/2 size of screen
                    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3)/2 ...
                        screenSz(4)]);
                    plot(grp{n}.vess{a}.time, nanmean(grp{n}.vess{a}.RBCV_norm), 'r');
                    hold on;
                    %check if enough data for an errorbar
                    if size(grp{n}.vess{a}.RBCV_norm,1)>=3
                        shadedErrorBar(grp{n}.vess{a}.time,nanmean(grp{n}.vess{a}.RBCV_norm), ...
                            getSEM(grp{n}.vess{a}.RBCV_norm,1),'lineProps', 'r');
                    end
                    plot([0 grp{n}.vess{a}.time(end)],[0 0], 'k'); %plot x axis
                    plot([prefs.timeAvgRBCV, prefs.timeAvgRBCV],[-0.15, 0.15], 'k');
                    xlim([0 grp{n}.vess{a}.time(end)]);
                    title('avg RBCV trace');
                    %save figure into expdir
                    figSave = 'RBCV_cutByCalc.png';
                    saveas(gcf, fullfile([expDir,filesep,figSave]));
                    close(gcf);
                    
                    %plot individual RBCV peaks and resp index
                    figure;
                    %make fig size of screen
                    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
                        screenSz(4)]);
                    for b = 1:size(grp{n}.vess{a}.RBCVpk,2) %loop RBCV traces
                        subplot(ceil(size(grp{n}.vess{a}.RBCVpk,2)/4),4,b);
                        plot(grp{n}.vess{a}.time,grp{n}.vess{a}.RBCV_norm(b,:), 'r');
                        hold on;
                        plot([0 grp{n}.vess{a}.time(end)],[0 0], 'k'); %plot x axis
                        plot([prefs.timeAvgRBCV, prefs.timeAvgRBCV],[-0.15, 0.15], 'k');
                        xlim([0 grp{n}.vess{a}.time(end)]);
                        if grp{n}.vess{a}.respInd(b) == 1
                            title([num2str(b),': resp']);
                        else
                            title([num2str(b),': non-resp']);
                        end
                    end
                    %save figure into expdir
                    figSave = 'RBCV_respInd.png';
                    saveas(gcf, fullfile([expDir,filesep,figSave]));
                    close(gcf);
                    
                else % after removing bad pks / pks with not enough time
                    
                    %check the vessel recording is contributing at least 5
                    %peaks/RBCV recordings - if not clear it 
                    grp{n}.vess{a} = []; 
                    disp('No calc pks detected, skipping...')
                    
                end %end of checking if any pks left
                
            else %RBCV and calc traces not same size, inform user
                
                disp('RBCV and calcium traces not same size, skipping...');
                
            end %end of check if calc and RBCV traces same size
            
        else %no RBCV mat file
            
            disp('No RBCV mat file in dir, skipping...');
            
        end %end of check if RBCV file exists
        
    end %end of loop dirs with tif file
    
end %end of loop grps

save('LScutByCalc_RBCV','grp','prefs');

end %end of function
