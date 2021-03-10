
function extractLScalcium(topDir, prefs)

%function to find a tif file 'calcium.tif' with the linescan calc trace
%if you cross the BV (diam) within calc trace, draw a box over this and
%clear it in imageJ - the code will turn this space into NaNs so won't
%affect results 
%get an average per frame using the sliding time window (set to same as
%extracted RBV trace for) - i.e. so the time series are the same size

%send into top dir, as this purely extracts data and saves in local exp dir
%- doesn't call any other functions (other than findFolder)

%if user doesn't input enough arguments
if nargin<2
    %ms used to calc the window sz for looking at RBCV
    %NB match this to the RBCV - so on same time scale
    prefs.windowSz = 80; %ms
    %thresh to multiply std by for binarising the calc image ch
    prefs.thresh = 0.5; 
end

%check for existence of the RBCV.tif file
find_tif_file = findFolders(topDir, '*calcium.tif');

screenSz=get(0,'Screensize');

%loop through all directories with tif file
for a = 1:size(find_tif_file,2)
    
    clearvars -except topDir prefs find_tif_file a screenSz; 
    
    %inform user of progress
    disp(['Looping experiments: ', num2str(a), '/', ...
        num2str(size(find_tif_file,2))]);
    
    %find the local folder containing the tif file
    [expDir,~] = fileparts(find_tif_file{1,a});
    
    %looking for .ini file, which contains parameter info
    find_ini_file = findFolders(expDir, '*.ini');
    x = cell2mat(find_ini_file);
    ini_file = ini2struct(x);
    clear x find_ini_file
    
    %get the number of frames, the width and height of the RBCV ch
    %this has been cropped around the correct part of the scan path, so
    %width and height will be different
    info = imfinfo(find_tif_file{1,a});
    width = info(1).Width;
    height = info(1).Height;
    num_images = numel(info);
    mspline = str2num(ini_file.x_.ms0x2ep0x2eline); %ms per line
    clear inf_file;
    
    %load the tif files for each parameter, across all frames
    %initialise an empty array to store the tiff file in as large dataset
    A = zeros(width,height,num_images);
    
    % Load in all the frames
    % inform user as process can be quite slow
    disp('loading images for all channels...');
    for k = 1:num_images %loop frames
        %disp every 100 frames to show progress
        if (mod(k,100) == 0)
            disp(['Frame ', num2str(k), '/', num2str(num_images)])
        end
        %load the line scan data
        A(:,:,k) = imread(find_tif_file{1,a}, k)';
    end %end of loop frames when loading calc ch
    
    %reorder the data:
    %put all frames on top of each other, so it becomes a 2D image
    %get the time x space line scan info
    rawLine = reshape(A, [width, height*num_images]);
    
    %the image will have zeros in it where the linescan went through a
    %vessel (for diam) and it has been cleared
    %replace the zeros with NaNs
    rawLine(find(rawLine==0))=NaN;
    
    %normalise the image so the values are between 0 & 1
    %then every line scan (i.e. across multiple sessions) will be on the 
    %same scale, and comparible 
    rawLine = rawLine-min(min(rawLine));
    rawLine = rawLine/max(max(rawLine));
    
    if false
        %save part of the linescan as a figure, so the user can see if this
        %normalising has worked okay
        %take first 100 points as would expect beginning of linescan to be
        %good, as it is how the user decided to attempt a line scan
        figure;
        imagesc(rawLine(:,1:100));
        title('first 100 points of linescan after normalising data');
    end
    
    %% step through the data in specified time window (e.g. 80ms), and 
    % average calcium trace
    
    data = rawLine'; 
    
    % NB Drew uses 40ms window but he uses this for arterioles, we needed a
    % larger window for our capillaries to capture more of the streaks
    % this was tested quite extensively, but user can change the window
    % size if it doesn't work for their needs
    % NOTE the window size must be divisible by 4! 
    %NB for HC vs V1 paper - all extracted at 80
    windowsize = round((prefs.windowSz/mspline)/4)*4; %pixels
    
    %define step size for sliding window
    %this is because the data will be averaged across the window, then step 1/4
    %of the window and average again - this means data is overlapped in
    %averaging by multiple windows, and provides a smoothing effect
    %divide the no of pixels in the linescan by 4
    stepsize=round(.25*windowsize); %pixels
    nlines = size(data,1); %time, i.e. no Lines
    npoints = size(data,2); %space, i.e. Pixels
    %calculate how many steps can be taken with the step size specified and
    %size of the data - this will be used for looping the data
    nsteps = floor(nlines/stepsize)-3;
    
    %loop for every step
    for k = 1:nsteps
        
        %some sort of time vector? (used to create time later) - use to get
        %sample rate
        time_pts(k) = 1+(k-1)*stepsize+windowsize/2;
        
        % take info from Continuous data
        % loop through data in steps of step size * 4
        % this will overlap by one time box on each loop
        % 0:T/4 time box = stepsize+windowsize/2 - i.e. first box
        % 4 steps later (i.e. 4 time boxes on, 4*stepsize)
        % - (k-1)*stepsize+windowsize
        data_hold = data(1+(k-1)*stepsize:(k-1)*stepsize+windowsize,:);
        
        %% Threshold image
        %binarise data - so any calcium signal is given a 1, and background
        %is given a 0
        % 1. set threshold: 
        %find mean intensity value of image
        meanim = nanmean(reshape(data_hold, [1, size(data_hold, 1)*size(data_hold,2)]));
        %find std intensity value of image
        stdim = nanstd(reshape(data_hold, [1, size(data_hold, 1)*size(data_hold,2)]));
        % 2. binarise pixels based on whether above or below thresh
        %set any background pixels (less than mean+thresh*std) to 0
        data_hold(data_hold<(meanim+(prefs.thresh*stdim))) = 0;
        %set any pixels which are nans to zero
        data_hold(isnan(data_hold)) = 0;
        %set any signal pixels (>= than mean+thresh*std) to 1
        data_hold(data_hold>=(meanim+(prefs.thresh*stdim))) = 1;
        %plot binarised calc image, if want to check 
        if false
            figure;
            subplot(211);
            imagesc(data_hold); 
            title('Binarised image'); 
        end
        %remove small noise pixels - just leave blocks of signal
        data_hold = bwareaopen(data_hold, 50);
        %plot binarised calc image, with noise speckles removed
         if false
            figure;
            subplot(212);
            imagesc(data_hold); 
            title('Binarised image, noise removed'); 
        end
        %nanmean over the window to get an average value for calc intensity
        %per sample pt 
        calcium(k) = nanmean(nanmean(data_hold));
    end %end of stepping through time windows 
    
    %create time vector for plots (in seconds) 
    time =(time_pts*mspline)/1000; 
    
    figure;
    %make fig size of screen
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
        screenSz(4)/2]);
    plot(time, calcium);
    title('linescan calcium trace');
    xlabel('time(s)');
    %save figure into expdir
    figSave = 'calciumTrace.png';
    saveas(gcf, fullfile([expDir,filesep,figSave]));
    close(gcf);
    
    %save variables
    prefs2output = prefs;
    matfile = fullfile(expDir, 'contData_ls_calc.mat');
    save(matfile,'calcium','time','prefs2output');
    
end %end of looping tif files

end %end of function