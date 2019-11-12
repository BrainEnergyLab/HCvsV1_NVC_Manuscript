
function [spikeEvents] = findSignalSpikesInTrace(data, fps, baseTim)

if nargin < 3
    baseTim = 3; %seconds
end

prePeakTime = round(baseTim*fps);

frames = [1:size(data,2)];

%% find all pts of data which are above a set threshold:

%find spikes
thresh = 0.1;
spikePeakLocs = find(data>=thresh*max(data));

if ~isempty(spikePeakLocs)
    %plot trace with detected spikes
    if false
        figure;
        plot(frames, data, 'g');
        hold on;
        plot(frames(:,spikePeakLocs),data(:,spikePeakLocs), 'ko');
    end
    
    %find start of peak
    %take mean of signal - 10% of mean, and find last value preceeding spike
    %that is less than mean threshold - this gives a ROUGH start pt for peak
    meanThresh = nanmean(data)-(thresh*nanmean(data));
    for n = 1:size(spikePeakLocs,2)
        ttt=find(data(:,1:spikePeakLocs(n))<=meanThresh);
        if ~isempty(ttt)
            spikeStartLocs(n) = ttt(end);
        else
            spikeStartLocs(n) = NaN;
        end
    end
    %remove NaNs
    nanInd=find(isnan(spikeStartLocs));
    spikeStartLocs(nanInd)=[];
    spikePeakLocs(nanInd)=[];
    
    % plot(frames(:,spikeStartLocs),data(:,spikeStartLocs), 'bo');
    
    %% now find only correct peak values, i.e. remove repeats
    
    %now we need to find all spikes which share a start pt, and remove the
    %'repeat peak' values - i.e. by only taking the peak with the largest value
    %from the cluster of peaks
    %look for the diff between all the start locs - so can ID which share start
    %pts and which have their own
    %need to check if first val has own start pt or not, if it is same as val
    %proceeding it - allocate it a zero
    if spikeStartLocs(2) - spikeStartLocs(1) == 0
        repeatSpikeInd = [0 diff(spikeStartLocs)];
    else %not same as val preceding it, so allocate it an arbitrary 5
        %put arb num to show we can take this spike - i.e. show it has its own
        %start pt
        repeatSpikeInd = [5 diff(spikeStartLocs)];
    end
    %make a switch index to see where need to search between to find highest
    %peak - this would pick up first val if i'd set it to 5
    [~,switchPtInd] = find(repeatSpikeInd);
    %next need to find largest peaks to correspond with start vals
    %add first value on
    % switchPtInd = [1 switchPtInd];
    %add last value on
    % switchPtInd = [switchPtInd size(repeatSpikeInd,2)];
    
    %retake the start values using only the non repeats
    vals2use = switchPtInd(find([1 diff(spikeStartLocs(switchPtInd))]));
    spikeStartLocs = spikeStartLocs(vals2use);
    
    %check for max detected peak for all peaks which share a start pt
    for n = 1:size(switchPtInd,2)
        if n == 1 %first peak 
            %go from first spike loc to last one before switch
            datapts = spikePeakLocs(1:switchPtInd(n+1)-1);
            %find which is the largest
            [~,maxInd]=max(data(:,datapts));
        elseif n+1 > size(switchPtInd,2) %last peak - think this is ok
            %go from final peak, until size of peaks
            datapts = spikePeakLocs(switchPtInd(n):size(spikePeakLocs,2));
            %find which is the largest
            [~,maxInd]=max(data(:,datapts));
        else
            %take all spike locs with same start pt
            datapts = spikePeakLocs(switchPtInd(n):switchPtInd(n+1)-1);
            %find which is the largest
            [~,maxInd]=max(data(:,datapts));
        end
        if size(maxInd,2) == 1
            spikeLoc_t(n) = datapts(maxInd);
        else
            maxInd = maxInd(1);
            spikeLoc_t(n) = datapts(maxInd);
        end
    end
    %take new peaks - with repeats deleted
    clear spikePeakLocs;
    spikePeakLocs = spikeLoc_t;
    clear spikeLoc_t;
    
    %% now have got the correct peaks - need to find correct START of peak:
    
    %for each 'start of peak' location - search trace before for nearest pt
    %preceding start which is higher than start, then take pt after that
    
    %plot trace with detected spikes
    if false
        figure;
        plot(frames, data', 'g');
        hold on;
        plot(frames(:,spikePeakLocs),data(:,spikePeakLocs), 'ko');
        plot(frames(:,spikeStartLocs),data(:,spikeStartLocs), 'bo');
        plot([frames(1) frames(end)],[meanThresh meanThresh],'k')
    end
    
    %check there is less than 3s between start of peak and peak
    %create index with 1 if the data is okay, and 0 if it needs removing
    for n = 1:size(spikeStartLocs,2)
        if spikePeakLocs(n) - spikeStartLocs(n) <= prePeakTime
            baseTimeInd(n) = 1;
        else
            baseTimeInd(n) = 0;
        end
    end
    %check if any have been classed as being to long between onset and peak
    %if so, remove those start of peak and peak locations from index
    if sum(baseTimeInd) < size(spikeStartLocs,2)
        spikeStartLocs(find(baseTimeInd == 0)) = [];
        spikePeakLocs(find(baseTimeInd == 0)) = [];
    end
    
    %check the min dist between peaks is longer than this base period (e.g. of
    %3s) - if not, only take largest peak
    minDistInd = [prePeakTime, diff(spikePeakLocs)];
    %any that are lower than the require min dist - need to test pos in ind and
    %position prior in ind, see which val is bigger
    ttt=find(minDistInd < prePeakTime);
    for n = 1:size(ttt,2) %loop the spikes with not enough prepeak time
        %check if the data for the peak detected is bigger than the previous
        %peak
        if data(:,spikePeakLocs(ttt(n))) > data(:,spikePeakLocs(ttt(n)-1))
            %if it is bigger, then update the previous peak to be equiv to this
            spikePeakLocs(ttt(n)-1) = spikePeakLocs(ttt(n));
        else
            %if it is smaller, then make this current peak the previous peak
            spikePeakLocs(ttt(n)) = spikePeakLocs(ttt(n)-1);
        end
    end
    %remove repeat values
    spikeStartLocs(find(([1 diff(spikePeakLocs)])==0))=[];
    spikePeakLocs(find(([1 diff(spikePeakLocs)])==0))=[];
    
    spikeEvents(1,:) = spikeStartLocs;
    spikeEvents(2,:) = spikePeakLocs;
    
else
    
    spikeEvents = []; 
    
end

end

