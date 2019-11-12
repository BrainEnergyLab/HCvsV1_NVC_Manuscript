function extractLineScanData(fname, prefs)

%top function to extract line scan RBCV and/or diameter data
%written by Kira & Orla, updated by Kira April 2019
%data should be saved as *RBCV.tif or *diam.tif, and folder should have
%spont or stim or VR in name - to categorise experiment

%INPUTS-
%fname = directory which contains RBCV.tif(s) and/or diam.tif(s)
%NB linescans require preprocessing to create RBCV.tif and/or diam.tif
%vesselCh = this is the PMT channel that the vessel is in, 1 for FITC, or 2
%for TR, it automatically assumes TR if user doesn't specify
%borderSz = this will crop the edges of the data, as there may be noise
%from image registration of the diam.tif, if not specified it will assume
%no border crop (this value is in pixels, and should be an integer,
%typically between 10-50, it will be removed from every edge of the image)

%OUTPUTS-
%no variables are outputted from the function, but .mat files and figures
%will be saved
%
%functions called inside this function for the line scan analysis:
%linescanVelocityAnalysis, GetVelocityRadon,
%linescanDiamAnalysis
%general functions which are called:
%findFolders, pixel4ls, avgDataOverSlidingWindow, findLocoEvents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set preferences if not all arguments are sent in
if nargin<2
    
    %specify if you want to extract ls RBCV, ls diam or both:
    prefs.analysis = 'RBCV'; %'both', 'RBCV', 'diam'
    
    %automatically assumes vessel is in channel 2
    prefs.vesselCh = 2; %TR
    %how many pixels to crop around the edges of the inputted data frame
    prefs.borderSz = 0; %pixels
    
    %spike remover:
    %calls a function to clean the trace for RBCV and diam
    prefs.removeSpikes = 1; %1 - will run it through code to remove 
    %plots on or off-
    prefs.plotRemovedSpikes = 1; %1 - will run it through code to remove 
    %the outlier detection method-
    prefs.method = 'mean'; % 'mean', 'median', 'gsed', 'quartile'
    %for more info follow this link:
    %https://uk.mathworks.com/help/matlab/ref/isoutlier.html
    
    %ms used to calc the window sz for looking at RBCV and diam - keep
    %consitent between both: 
    prefs.windowSz = 40; %ms
    
    %this is for the velocity only:
    %remove bad frames from the calculate velocity trace
    %call on findLocoEvents function to find large jitters in the signal
    %these large noise spikes will be removed
    %set preferences before call function, as these will be different in
    %this case to when want to detect loco events
    prefs.minDist = 100; %min dist between noise spikes
    %whether you want to see the plots for how it has found the noise
    %will set this to zero, as you will see a plot of the velocity before
    %and after noise removal anyway
    prefs.plotFlag = 0;
    %for locomotion spikes in the signal that last less than 1 sec are
    %removed, as these are not true loco events, however in this case we
    %want to detect ALL signal spikes to remove them so turn this off
    prefs.flickerFlag = 0;
    
end %end of check how many arguments sent into func

%inform user of border crop that will be used
disp(['Border crop of ', num2str(prefs.borderSz), ' pixels']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% catches to check what user wants and whether there is any data 
%check whether any relevant tif files in top dir inputted
if strcmp(prefs.analysis,'both') 
    find_RBCVtif_file = findFolders(fname, '*RBCV.tif');
    find_diamtif_file = findFolders(fname, '*diam.tif');
    %loop higher number of expDirs, i.e. in case you have folders where
    %some have diam.tif and not RBCV.tif 
    if size(find_RBCVtif_file,2) ~= size(find_diamtif_file,2)
        disp(['Warning: you dont have same num of RBCV and ', ...
        'diam tif files in your inputted dir']); 
    end
elseif strcmp(prefs.analysis,'RBCV')  
    find_RBCVtif_file=findFolders(fname, '*RBCV.tif');
elseif strcmp(prefs.analysis,'diam')
    find_diamtif_file=findFolders(fname, '*diam.tif');
else
    disp('You need to specify if want RBCV, diam or both');
    disp('Exiting function...');
    return;
end
   
    %% RED BLOOD CELL VELOCITY DATA:

    %check if user requested ls RBCV extraction
    if strcmp(prefs.analysis,'both') || strcmp(prefs.analysis,'RBCV')
        if isempty(find_RBCVtif_file) %check if any RBCV tif files in top dir
            disp('No RBCV tif files detected');
        else
            %loop through tif files found
            for a = 1:size(find_RBCVtif_file,2) %loop exp dirs with tif files
                %find individual exp directories
                [expDir,~]=fileparts(find_RBCVtif_file{1,a});
                %check if the data has already been extracted
                %this is quite a slow function, so don't want to call it if it has
                %already been run
                if exist([expDir, filesep,'contData_ls_RBCV.mat'])
                    disp('line scan velocity data created...');
                else
                    %if data has not been extracted, call function to look for velocity
                    disp('running line scan velocity finder...');
                    linescanVelocityAnalysis(expDir, prefs);
                end %end of check if the velocity code has already been run
            end %end of looping through tif files for RBCV
        end %end of check if any RBCV tif files found
    end %end of check if user requested ls RBCV extraction
    
    %% LINE SCAN DIAMETER DATA:

    %check if user requested ls diam extraction
    if strcmp(prefs.analysis,'both') || strcmp(prefs.analysis,'diam')
        if isempty(find_diamtif_file) %check if any RBCV tif files in top dir
            disp('No diam tif files detected');
        else
            %loop through tif files found
            for a = 1:size(find_diamtif_file,2) %loop exp dirs with tif files
                
                %find individual exp directories
                [expDir,~]=fileparts(find_diamtif_file{1,a});
                
                %check if the data has already been extracted
                %this is quite a slow function, so don't want to call it if it has
                %already been run
                if size(find_diamtif_file,2)>0
                    if exist([expDir, filesep,'contData_ls_diam.mat'])
                        disp('line scan diam data created...');
                    else
                        disp('running line scan diameter finder...');
                        linescanDiamAnalysis(expDir, vesselCh, prefs);
                    end
                end
                
            end %end of loop through mat files
        end
    end %end of check if user requested ls diam extraction

end %end of function
