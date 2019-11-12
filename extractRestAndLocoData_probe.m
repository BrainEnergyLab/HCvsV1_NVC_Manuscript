function extractRestAndLocoData_probe(fname, prefs)

%function written by Kira, May 2018
%
%function takes all rest and loco periods, and extracts average haem
%responses across all these epochs for each recording
%this saves the data out into the top dir in a cell containing the multiple
%'trials', i.e. for each rest/loco epoch and for each recording
%can be used to compare two groups, e.g. V1 vs HC, or APOE3 vs APOE4 - the
%function must be run separately on each group directory to store them as
%separate mat files
%
%requires user to have previously run 'extractProbeData' which outputs
%contData.mat
%
%INPUTS-
%fname : top directory with contData.mat files within
%prefs : preferences for whether want plots or to remove locomotion periods
%lasting less than 1 second, etc.
%if user specifies the prefs themselves, all must be included
%
%OUTPUTS-
%nothing is outputted from function into workspace, but will save a matfile
%with restperiod data inside
%
%other functions needed for this to run:
%findFolders, cleanLoco, findLocoEvents

if nargin < 2
    %plot the loco trace with onsets and offsets marked
    %automatically set to not export these
    prefs.plotFlag = 1;
    %flag for if want to remove the flickers in the signal (not loco as
    %shorter than 1s)
    %automatically set to remove flickers (=0)
    prefs.flickerFlag = 0;
    %minimum distance between locomotion events - in seconds
    %will be converted to frames wthin code
    prefs.minDistSec = 3; %seconds
    %minimum duration for rest or loco event, must be >=2
    %automatically set at 2s, as 1s will be taken from beginning when haem
    %may be ramping up or down in response to change in loco
    prefs.minDurEvent = 2; %seconds
end

%look for all data files to be analysed
findmatfile=findFolders(fname, 'contData.mat');

%loop through each experimental directory containing the mat file
for a = 1:size(findmatfile,2)
    
    [expDir,matFile]=fileparts(findmatfile{1,a});
    
    %inform user of progress
    disp(['processing file ',num2str(a),'/',num2str(size( ...
        findmatfile,2))]);
    
    %load the data
    load(fullfile(expDir, matFile));
    
    %process the locomotion data
    [locomotion]=cleanLoco(movement);
    
    %update the prefs.minDist (for min dist between loco events) to convert
    %it into frames
    prefs.minDist = prefs.minDistSec * fps;
    %find all the locomotion events
    %make sure the locoEvents and restEvents vars are cleared from last
    %loop
    clear locoEvents restEvents;
    [locoEvents]=findLocoEvents(locomotion,fps,prefs);
    %locoEvents has 3 dimensions: 1 = onset, 2=offset, 3=duration - all in
    %frames
    
    if ~isempty(locoEvents) %check if any loco events have been detected
        
        %% extract loco periods
        
        %remove any loco periods that last for less than 2 seconds
        %as want to remove the 1st second post-rest to give haem responses
        %time to react to onset of loco
        for b = 1:size(locoEvents,2) %loop loco events
            %check if they last less than 2s
            if locoEvents(3,b) < prefs.minDurEvent*fps
                %set too short loco events to NaN
                locoEvents(1,b)=NaN;
                locoEvents(2,b)=NaN;
                locoEvents(3,b)=NaN;
            end %end of duration check
        end %end of loop loco events
        
        %remove NaNs from loco events index:
        locoEvents_t=locoEvents(1,:);
        locoEvents_t(isnan(locoEvents_t))=[];
        locoEvents_tt=locoEvents(2,:);
        locoEvents_tt(isnan(locoEvents_tt))=[];
        locoEvents_ttt=locoEvents(3,:);
        locoEvents_ttt(isnan(locoEvents_ttt))=[];
        clear locoEvents;
        locoEvents(1,:)=locoEvents_t; clear locoEvents_t;
        locoEvents(2,:)=locoEvents_tt; clear locoEvents_tt;
        locoEvents(3,:)=locoEvents_ttt; clear locoEvents_ttt;

        %% extract rest periods
        
        %check if have any loco periods left after removed too short ones
        if ~isempty(locoEvents) %check if any loco events have been detected
            
            %find periods when there is no locomotion, i.e. 'rest periods':
            for b = 1:size(locoEvents,2) %loop loco events
                if b == 1 %for first rest period
                    if locoEvents(1,b) > 2 %check enough room for a base
                        %start of rest event
                        restEvents(1,b) = 1;
                        %start of loco offset - 1 frame
                        restEvents(2,b) = locoEvents(1,b)-1;
                        %calculate duration of rest event and put in dim 3
                        restEvents(3,b) = restEvents(2,b)-restEvents(1,b);
                    else
                        %if there is not enough room for a rest period b4
                        %the first loco event, set rest event 1 as NaN
                        restEvents(1,b) = NaN;
                        restEvents(2,b) = NaN;
                        restEvents(3,b) = NaN;
                    end %end of check enough base time for rest event
                else %for rest of rest periods
                    if locoEvents(b) < size(data,2) %check data fits on end
                        %take offset of previous loco event and add 1 to get
                        %start of rest event:
                        restEvents(1,b) = locoEvents(2,b-1)+1;
                        %take offset of current loco event and -1 frame
                        %of rest event:
                        restEvents(2,b) = locoEvents(1,b)-1;
                        %calculate duration of rest event
                        restEvents(3,b) = restEvents(2,b)-restEvents(1,b);
                    else
                        %if there is not enough time for the final rest
                        %event to fit on the end, then add a nan
                        restEvents(1,b) = NaN;
                        restEvents(2,b) = NaN;
                        restEvents(3,b) = NaN;
                    end %end of check if data fits on end
                end %end of check if 1st rest period
            end %end of loop loco events
            
            %remove any rest periods lasting less than 2 seconds
            for b = 1:size(restEvents,2) %loop rest events
                if restEvents(3,b) < prefs.minDurEvent*fps %check duration
                    restEvents(1,b)=NaN;
                    restEvents(2,b)=NaN;
                    restEvents(3,b)=NaN;
                end %end of check duration
            end %end of loop rest events
            
            %remove NaNs - figure out a better way!!!
            %rest
            restEvents_t=restEvents(1,:);
            restEvents_t(isnan(restEvents_t))=[];
            restEvents_tt=restEvents(2,:);
            restEvents_tt(isnan(restEvents_tt))=[];
            restEvents_ttt=restEvents(3,:);
            restEvents_ttt(isnan(restEvents_ttt))=[];
            clear restEvents;
            restEvents(1,:)=restEvents_t; clear restEvents_t;
            restEvents(2,:)=restEvents_tt; clear restEvents_tt;
            restEvents(3,:)=restEvents_ttt; clear restEvents_ttt;
            
            %take 1 sec after loco starts - i.e. to give haem responses time to
            %ramp up (if they are doing)
            %all events will have enough time, as any <2s in duration have been
            %removed
            for b = 1:size(locoEvents,2) %loop loco events
                locoEvents(1,b) = locoEvents(1,b)+(1*fps); %add 1s to onset
                locoEvents(3,b) = locoEvents(2,b)-locoEvents(1,b); %update dur
            end %end of loop loco events
            
            %remove 1st second of rest period, as loco may still be
            %affecting haem responses from previous loco events
            for b = 1:size(restEvents,2) %loop rest events
                restEvents(1,b)=restEvents(1,b)+(1*fps); %add 1s to start
                restEvents(3,b) = restEvents(2,b)-restEvents(1,b); %duration
            end
            
            %% average data during loco or rest epochs
            %take the average 'data' during each rest period - effectively
            %giving multiple trials
            for b = 1:size(restEvents,2) %loop rest periods
                %extract data during each individual rest period and average
                dataRestTrials{a}(:,b) = nanmean(data(:,restEvents(1,b) ...
                    :restEvents(2,b)),2);
            end %end of loop rest periods
            %also extract data during individual loco events
            for b = 1:size(locoEvents,2) %loop loco events
                dataLocoTrials{a}(:,b) = nanmean(data(:,locoEvents(1,b) ...
                    :locoEvents(2,b)),2);
            end %end of loop loco events
            
        else
            %if no loco events, then the entire recording is a 'rest trial'
            %so can average across it for rest trials
            dataRestTrials{a} = nanmean(data,2);
            %no loco events
            dataLocoTrial{a} = NaN;
            
        end %end of check for any loco events (after remove too short ones)
        
    else %check if any loco events
        %if no loco events, then the entire recording is a 'rest trial' so
        %can average across it for rest trials
        dataRestTrials{a} = nanmean(data,2);
        %no loco events
        dataLocoTrial{a} = NaN;
    end %end of check if any loco events have been detected
    
end %end of looping mat files

%save the data which contains all 'trials' i.e. the average for each rest
%or loco epoch, and all sessions
matfile = fullfile(fname, 'Avg4RestnLocoPeriods_allTrials_allSess.mat');
save(matfile, 'dataLocoTrials','dataRestTrials','findmatfile','-v7.3');

end %end of function

