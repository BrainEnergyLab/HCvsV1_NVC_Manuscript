
%written by Kira, Nov 2017

%function to find diameter from linescans
%send in: experimental directory, and name of tif file with data
%NB// tif file with vessel diam needs to be cropped around the vessel
%(remove artefacts from background), and remove bad frames

%plots the intensity curve for each line of scan (i.e. along x), and each
%frame, and finds the FWHM of this curve
%saves extracted variables in exp dir

function findDiam_lineScan(exp_dir, vesselCh, tifName)

close all; 

%if dont specify tif file - looks for cleanLine_diam.tif
if nargin<2
    tifName='cleanLine_diam.tif';
    vesselCh=2; %automatically assumes TR, specify 1 if used FITC
elseif nargin<3
    tifName='cleanLine_diam.tif';
end

%look for the tif file specified within the experimental directory
files = findFolders(exp_dir, tifName);
%load other channels for loco, calc, stim
if vesselCh==2
    find_calc=findFolders(exp_dir, 'ch1_diam_crop*.tif');
else
    find_calc=findFolders(exp_dir, 'ch2_diam_crop*.tif');
end
find_loco=findFolders(exp_dir, 'ch3_diam_crop*.tif');
find_stim=findFolders(exp_dir, 'ch4_diam_crop*.tif');

%notify user if no diam tif files found, and exit code
if size(files,2)==0
    disp('no line scan diam tif file detected...')
    return
end

%loop through diam tif files found
for l = 1:size(files,2)
    
    %get info from sciscan outputs, i.e. num frames, and dimensions of image
    fname = files{l};
    [path,name,~] = fileparts(fname);
    info = imfinfo(fname);
    num_frames = numel(info);
    dims = [info(1).Height, info(1).Width];
    
    %change directory to experimental directory (i.e. for saving)
    cd(path);
    
    %open ini file to find one pixel size in metres
    find_ini_file = findFolders(cd, '*.ini');
    x= cell2mat(find_ini_file);
    ini_file = ini2struct(x);
    clear x find_ini_file
    %get size info:
    %convert pixel size from metres to um (multiply by a million)
    pxsz_um = str2num(ini_file.x_.x0x2epixel0x2esz)*1000000;
    %get timing info:
    fps = str2num(ini_file.x_.frames0x2ep0x2esec);
    mspline=str2num(ini_file.x_.ms0x2ep0x2eline);
    lps=str2num(ini_file.x_.lines0x2eper0x2esecond);
    
    %load other channels
    calcium = find_calc{1,l};
    loco = find_loco{1,l};
    stim = find_stim{1,l};
    
    %load diam tif file data into matlab (var A)
    rawIm = zeros(num_frames, dims(1), dims(2));
    for i = 1:num_frames
        rawIm(i,:, :) = imread(fname, i);
        %load other channels
        calc_ch(i,:,:)=imread(calcium,i)';
        loco_ch(i,:,:)=imread(loco,i)';
        stim_ch(i,:,:)=imread(stim,i)';
    end
    clear linescan calcium loco stim;
    
    %read mean intensity from each channel (calc, loco, stim)
    calcium=zeros([1,num_frames*dims(1)]); 
    movement=zeros([1,num_frames*dims(1)]); 
    stim=zeros([1,num_frames*dims(1)]); 
    for h = 1:num_frames
            calcium(1+((h-1)*dims(1)):h*dims(1))=repmat(nanmean(nanmean(squeeze(calc_ch(h,:,:)),1),2),[1, dims(1)]);
            movement(1+((h-1)*dims(1)):h*dims(1))=repmat(nanmean(nanmean(squeeze(loco_ch(h,:,:)),1),2),[1, dims(1)]);
            stim(1+((h-1)*dims(1)):h*dims(1))=repmat(nanmean(nanmean(squeeze(stim_ch(h,:,:)),1),2),[1, dims(1)]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get a rough diameter for the vessel
    
    %threshold the image (i.e. to separate it from the background)
    meanim = mean(rawIm(:));
    stdim = std(rawIm(:));
    thresh=meanim+stdim*0.1;
    
    %binarise image
    bwIm=zeros(size(rawIm));
    for i=1:size(rawIm,1)
        for j=1:size(rawIm,2)
            for k=1:size(rawIm,3)
                if rawIm(i,j,k)>=thresh
                    bwIm(i,j,k)=1;
                    %vessel is labelled 1
                else
                    bwIm(i,j,k)=0;
                    %background is labelled 0
                end
            end
        end
    end
    
    %filters to fill holes, smooth edges, etc.
    clear i;
    for i=1:num_frames
        %fill any holes
        bwIm(i,:,:) = imfill(squeeze(bwIm(i, :, :)),'holes');
        %remove parts of image which are too small (i.e. below certain
        %pixel size) -- set at 300 pixels for now...
        bwIm(i,:,:) = bwareaopen(squeeze(bwIm(i,:,:)), 300);
        %any gaps inside two pos pixels, are also made pos (Morphological closing)
        bwIm(i,:,:) = bwmorph(squeeze(bwIm(i, :, :)),'close',100);
        %apply different Gaussian smoothing filters to images using imgaussfilt
        %Gaussian smoothing filters are commonly used to reduce noise
        %         bwIm(i,:,:)=imgaussfilt(squeeze(bwIm(i,:,:)),0.75);
    end
    
    %create skeleton
    closure = bwmorph(squeeze(mode(bwIm,1)),'close',Inf);
    closure = bwmorph(closure,'spur',5);
    closure = bwmorph(closure,'thin',Inf);
    skeleton = bwmorph(closure,'clean',5);
    %get length of skeleton
    lengthskel=sum(sum(skeleton)); %i.e. find all positive values
    
    %get a rough diameter of the vessel
    %use 1st frame
    %NB all vessel areas labelled as 3, therefore need to divide ans by 3
    %     vesselArea=sum(sum(squeeze(bwIm(1,:,:))));
    vesselArea=sum(sum(squeeze(mode(bwIm,1))));
    %Divide the area by the length to give you the mean width
    %in pixels
    roughdiameter_px=vesselArea/lengthskel;
    
    clear bwIm closure skeleton lengthskel vesselArea;
    
    %convert proxy diam changes (now in pixels) to um
    roughdiameter_um=roughdiameter_px*pxsz_um;
    
    if roughdiameter_um<11
        vessellabel='capillary';
    else
        vessellabel='arteriole';
    end
    
    %save rough diameter out - can use to workout if it is an
    %arteriole/capillary etc.
    save('roughdiameter','roughdiameter_um','roughdiameter_px','vessellabel');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %loop through frames/rows of the image to find full width half max
    clear i j;
    for i=1:num_frames %loop through frames
        for j=1:size(rawIm,2) %i.e. traveling down the vessel in linescan
            
            %extract each row of the linescan, for each frame
            %Smooth the data with the loess method using a span of 50% of the total number of data points
            data=smooth(squeeze(rawIm(i,j,:)),0.4,'loess');
            %             data=smooth(squeeze(rawIm(i,j,:)),0.5,'loess');
            
            % Find the half max value.
            halfMax = (min(data) + max(data)) / 2;
            maxIndex= find(data==max(data),1);
            % Find where the data last reaches the half max in the rise
            index1 = 1+find(data(1:maxIndex,:) <= halfMax, 1, 'last');
            % Find where the data last rises above half the max in the fall
            index2 = maxIndex + find(data(maxIndex+1:end,:) <= halfMax, 1, 'first');
            if ~isempty(index1) && ~isempty(index2)
                fwhm_px(i,j) = index2-index1 + 1; % FWHM in indexes.
            else
                fwhm_px(i,j)=NaN;
            end
            
            %output plots every 50 frames, i.e. to check it is working
            if (mod(i,100)==0) && (mod(j,50)==0)
                %plot the distribution of the data, and the indices for fwhm found
                figure;
                plot(squeeze(rawIm(i,j,:)),'k');
                hold on;
                plot(data,'b', 'LineWidth',2);
                plot([index1,index1],[halfMax,halfMax],'go');
                plot([index2,index2],[halfMax,halfMax],'ro');
                title(['frame:',num2str(i),' skelpt:',num2str(j)]);
            end
            
            
        end
    end
    
    % Put all frames on top of each other, so it becomes a 2D image.
    % get the time x space line scan info
    fwhm_px = reshape(fwhm_px, [1, dims(1)*num_frames]);
%     calcium = reshape(calc_int, [1, dims(1)*num_frames]);
%     movement = reshape(loco_int, [1, dims(1)*num_frames]);
%     stim = reshape(stim_int, [1, dims(1)*num_frames]);
    clear calc_ch loco_ch stim_ch;
    
    %width represents the length of the line you drew
    %multiply height x num_images to recreate the line scan
    
    %width (i.e. x) is the number of pixels in the line - i.e. the size of the line you drew
    %height (i.e. y) - each row represents a scan through the line
    %the num_images (i.e. frames) multiplied by length of the y axis is equivalent to how many time
    %the line was scanned -- this is just how sciscan stores the info, i.e. can put it back into 2 dimensions
    
    %smooth the FWHM pixel data
    disp('smoothing data...');
    fwhm_px_smooth=smooth(fwhm_px);
    disp('finished smoothing');
    
    if size(fwhm_px_smooth,1)>size(fwhm_px_smooth,2)
        fwhm_px_smooth=fwhm_px_smooth';
    end 
    
    %apply low pass filter to data
    d =  fdesign.lowpass('Fp,Fst,Ap,Ast', 0.1, 2, 0.5, 70, lps);
    Hd = design(d, 'equiripple');
    fwhm_px_lpf = filter(Hd, fwhm_px_smooth);
    
    %medfilt to remove offtrend spikes
    fwhm_px_lpf = medfilt1(fwhm_px_lpf,3);
    %remove outliers
    fwhm_px_lpf=deleteoutliers(fwhm_px_lpf);
    
    %remove outliers
    % data = the variable name of your array
    stdDev = std(fwhm_px_lpf(:)); % Compute standard deviation
    meanValue = mean(fwhm_px_lpf(:)); % Compute mean
    zFactor = 1.5; % or whatever you want.
    % Create a binary map of where outliers live.
    outliers = abs(fwhm_px_lpf-meanValue) > (zFactor * stdDev);_dir 
    
    %convert pixels to microns
    %NB use smoothed data
    fwhm_um=fwhm_px_smooth*pxsz_um;
    
    %%%%% process stim and loco
    %process stim so it is in binary
    stimOnIndex=find(stim(:,50:(end-100))>mean(stim)+std(stim))+50; %find when stim is on
    binaryData=zeros(size(stim));
    binaryData(:,stimOnIndex)=1; %convert to binary (i.e. 1)
    clear stim;
    stim=binaryData;
    clear binaryData stimOnIndex;
    %loco
    %smooth movement when extract locomotion for cutting (but not for raw)
    disp('smoothing loco...');
    [locomotion]=cleanLoco(smooth(movement,0.1,'loess'));
    disp('finished smoothing loco...');
    locomotion_raw=abs(diff(movement));
    locomotion_raw=[locomotion_raw,0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% TO FIX: why are these times slightly different??? - BUT not always
    %create a time vector
    time_frames=[1:num_frames]/fps;
    disp(['time(s) using frames: ' num2str(time_frames(end))]);
    time_lines=([1:size(fwhm_px,2)]*mspline)/1000;
    %/1000 to convert from ms to seconds
    disp(['time(s) using lines: ' num2str(time_lines(end))]);
    %out by fps-1.1s
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    %find size of screen (for figure sizing)
    screenSz=get(0,'Screensize');
    %make fig 1/2 size of screen
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);
    subplot(411);
    plot(time_lines, fwhm_px_smooth,'r');
    title('diam calculated from fwhm at each row of linescan');
    xlabel('time(s)');
    ylabel('pixels');
    xlim([time_lines(1) time_lines(end)]);
    subplot(412);
    plot(time_lines, locomotion_raw,'b');
    xlim([time_lines(1) time_lines(end)]);
    title('locomotion');
    xlabel('time(s)'); ylabel('locomotion,A.U.');
    subplot(413);
    plot(time_lines, stim,'k');
    xlim([time_lines(1) time_lines(end)]);
    title('stim');
    xlabel('time(s)');
    subplot(414);
    plot(time_lines, calcium,'g');
    xlim([time_lines(1) time_lines(end)]);
    title('calcium');
    xlabel('time(s)');
    saveas(gcf,'cont_linescan_diam.png');
    
    %save data
    save('contData_linescanDiam','fwhm_px','fwhm_px_smooth','fwhm_um','time_frames','time_lines',...
        'movement','locomotion','locomotion_raw','stim','calcium','fps','mspline','lps');
    
    
end %end of looping through linescan files

end %end of function