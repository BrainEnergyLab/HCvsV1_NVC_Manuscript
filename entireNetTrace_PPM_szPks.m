
clear all; close all; 
%'D:\Dropbox (Brain Energy Lab)\Everything\Kira\Projects\HCvsV1\data_sortedByAnalysis\zoomedOutCalc'
fname=cd; 

%count peaks prefs:
%set std multiplier for thresh for ppm calc
prefs.stdMultiplier = 2; %2, std
%normalise calc trace for detecting peaks
prefs.normaliseTrace = 1;
%time to search either side of the calc peak to find the max
%i.e. WM calc is slower than excitatory, so may need to increase
prefs.searchPer = 0.5; %seconds, 0.5 for excitatory, 1 for WM

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
        'cell_sig_norm.mat');
    %inform user of progress
    disp(['loading data for ', grpNames{1,n}, ' grp...']); 
    
    clear a; 
    for a = 1:size(findmatfiles,2) %loop through expDirs with data
        
        clear cell_sig fps peakLocs time; 
        %find individual exp dir
        [expDir,matfile]=fileparts(findmatfiles{1,a});
        %load data
        load(fullfile([expDir,filesep,matfile]));
        %added next line to get fps
        cd(expDir); 
        if exist([expDir, filesep, 'stimLocoTraces.mat'],'file')
            load('stimLocoTraces.mat', 'fps'); 
        else
            disp('running func to extract and process loco and stim chs');
            %call function
            extractLocoStimCh(expDir);
            load('stimLocoTraces.mat', 'fps'); 
        end
        
        cd(fname); 
        calcium = nanmean(cell_sig,1); 
        time = [0:size(calcium,2)-1]/fps;
        
        [peakLocs] = findPeaksInTrace(expDir,calcium,fps,prefs);
        
        grp{n}.expDir{a} = expDir; 
        
        %average PPM
        grp{n}.ppm_netTrace(a) = size(peakLocs.locs,2) / (time(end)/60);
        
        %average size of calc peak
        grp{n}.calcPeakVals{a} = calcium(:,peakLocs.locs);
        grp{n}.calcPeakVals_mean(a) = nanmean(calcium(:,peakLocs.locs));
        
        %correlation between ROIs
        r = corrcoef(cell_sig');
        %r will be # ROIs x # ROIs in size
        %put a NaN where the cell correlates with itself
        r(r==1) = NaN;
        %save the correlation between cells into the struct
        grp{n}.corr_all{a} = r;
        %get the mean correlation
        grp{n}.corr_avg(a) = nanmean(r(:)); 
        clear r; 

    end
    
end

%PPM:
disp('PPM:');
[~,p] = ttest2(grp{1}.ppm_netTrace,grp{2}.ppm_netTrace)

%peak sizes:
disp('Pk szs:');
[~,p] = ttest2(grp{1}.calcPeakVals_mean,grp{2}.calcPeakVals_mean)

%peak sizes:
disp('Correlation between cells:');
[~,p] = ttest2(abs(grp{1}.corr_avg),abs(grp{2}.corr_avg))

save('NEW_zoomOutCalc_PPM_corr','grp'); 








