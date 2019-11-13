
function [peakLocs] = findPeaksInTrace(expDir,data,fps,prefs)

% function created by Kira, April 2019 

% takes a continuous data trace, and finds all the spikes in it, NB if you
% send in a trace with multiple ROIs, each ROI has its own threshold 
% requires expDir to save an important fig - i.e. so user can check they
% are happy with the peaks which it has detected 

%INPUTS-
%expDir       :   exp dir to save figs into which show peaks found
%data         :   the continuous data trace with spikes in the signal, i.e.
%                 typically calcium, this can have multiple ROIs in the
%                 first dim and it will loop these
%fps          :   frames per second - used for timing info 
%prefs        :   these are to set the threshold for classifying a peak as
%                 responsive or not responsive, and also for deciding
%                 whether to normalise the data trace when looking for
%                 peaks, i.e. so it is detrended and spikes will all be
%                 above the set threshold - may not need to do this if the
%                 data sent in is ALREADY normalised 

%OUTPUTS-
%peakLocs     :   this is a structure where each struct level is a separate 
%                 ROI - inside each struct level are the peak locations
%                 across the cont trace for each ROI

%% set prefs
if nargin < 4
    
    %set std multiplier for thresh for ppm calc
    prefs.stdMultiplier = 2; %std
    %normalise calc trace for detecting peaks only
    prefs.normaliseTrace = 1; %0 dnt norm calc trace, 1 do norm
    %time to search either side of the calc peak to find the max
    %i.e. WM calc is slower than excitatory, so may need to increase 
    prefs.searchPer = 0.5; %seconds, 0.5
    
end

%% catches to exit function
if ndims(data)>2
    disp('Data must be 2D');
    disp('Exiting function');
    return;
end
if size(data,1) > size(data,2)
    disp('ROIs must be in first dim');
    disp('Exiting function');
    return;
end


%% detect peaks over time

%create a time vector
time = [0:size(data,2)-1]/fps;

%open figure to plot peak detector results - so can eyeball if it is
%working ok
figure;
screenSz=get(0,'Screensize');
%make fig size of screen
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);

for b = 1:size(data,1) %loop ROIs
    
    %% find the peak locations across the whole data trace
    
    %temporarily normalise the data trace - remove drift etc.
    %just for finding peaks in the signal
    if prefs.normaliseTrace == 1
        data_ttt = data(b,:);
        data_ttt = data_ttt - nanmin(data_ttt);
        data_ttt = data_ttt / nanmax(data_ttt);
        data_ttt = data_ttt - smooth(data_ttt,100)';
    else
        data_ttt = data(b,:);
    end
    
    %set peak threshold invidually for each ROI trace
    thresh = nanmean(data_ttt) + (prefs.stdMultiplier ...
        *nanstd(data_ttt));
    
    %find where the data trace exceeds the set threshold
    locs = find(data_ttt>thresh);
    %only take locations separated by 0.5s, as don't want to take same
    %peak repeatedly
    locs = locs([find(diff(locs)>fps/2),size(locs,2)]);
    
    %search 0.5s either side of peak, and find the actual max val
    %i.e. so it is actually located at peak and not on the way up/down
    %from the peak
    for c = 1:size(locs,2) %loop peak locations
        %find time 0.5s before and 0.5s after peak
        start = locs(c)-round(fps*prefs.searchPer);
        %catches to check new start and stop pts dont exceed size of
        %data
        if start <= 0
            start = 1;
        end
        stop = locs(c)+round(fps/2);
        if stop > size(data_ttt,2)
            stop = size(data_ttt,2);
        end
        %look within this search period and find the max
        [~,maxLoc] = max(data_ttt(:,start:stop));
        %replace peak location with new one that is the max val within
        %this search period
        %create temp var with new peak location in
        locs(c) = (start+maxLoc)-1;
    end %end of loop peak locs
    
    % plot peak finder if only 1 ROI/average trace
    if size(data,1) == 1
        title('normalised data trace with thresh');
        hold on;
        %normalised data trace
        plot(time,data_ttt,'k');
        %threshold
        plot([time(1) time(end)],[thresh thresh], 'b', 'LineWidth', 2);
        %peak locations
        plot(time(:,locs),data_ttt(:,locs), 'ro', 'LineWidth', 2);
        xlim([time(1) time(end)]);
        xlabel('Time(s)');
        ylabel('\Delta/ f/f');
    else
        %plot the peak finder and trace for first 20 ROIs, so can check
        %they look right
        if b <= 10
            subplot(5,2,b);
            title(['normalised calc with thresh, ROI: ', num2str(b)]);
            hold on;
            %normalised data trace
            plot(time, data_ttt,'k');
            %threshold
            plot([time(1) time(end)],[thresh thresh], 'b', ...
                'LineWidth', 2);
            %peak locations
            plot(time(:,locs), data_ttt(:,locs), 'ro', ...
                'LineWidth', 2);
            xlim([time(1) time(end)]);
            xlabel('Time(s)');
            ylabel('\Delta/ f/f');
        end
    end 
    
    %save out all the peak locations detected (to output from func):
    %separate in a structure by ROI 
    peakLocs(b).locs = locs; 
    
end %end of looping ROIs

%save figure
saveas(gcf, fullfile([expDir,filesep,'calcPeaksDetected.png']));
close;

end %end of func
