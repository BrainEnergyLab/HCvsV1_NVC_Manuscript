
function extractRBCflux(topDir, prefs)

% Function to extract RBC flux 
% Must have a binarised tif file of RBCV in exp
% dir called 'RBCV_binary 
% INPUTS: 
% topDir = top directory containing all experiments with the required files
% inside, NB will search specifically for the RBCV_binary.tif file 
% prefs = optional, if these aren't sent in the prefs below will be used,
% pref inputs are the window sz, the plot print prefs, and the threshold
% for counting a RBC as present or not (see below...) 

if nargin<2
   %NB window sz - keep this as 1000, as RBC flux unit is RBC/s, 
   %and this is 1000ms, i.e. 1s 
   prefs.windowSz_ms = 250; %1000; %1000ms = 1s
   prefs.plots = 1; %1 plots on, 0 plots off 
   %image binarised between 0 and 255 - threshold set just below half way 
   %(125.5) - have checked this by running plots in examples 
   prefs.thresh = 175; %127.5; %175 
end

screenSz=get(0,'Screensize');

%call function to find the tif file across the dirs in the fname
find_tif_file = findFolders(topDir, 'RBCV_binary.tif');

for n = 1:size(find_tif_file,2) %loop all dirs with the data in

    clearvars -except topDir prefs find_tif_file n screenSz; 
    
    disp(['Looping expDir ', num2str(n), '/', num2str(size(find_tif_file,2))]);
    
    %find individual exp dirs for loading and saving data
    [expDir,~] = fileparts(find_tif_file{1,n});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% need a binarise tif file from imageJ for flux and haematocrit
%%%%% crop these lots so dnt just send in crap - need quality control!

%looking for .ini file, which contains parameter info
find_ini_file = findFolders(expDir, '*.ini');
x = cell2mat(find_ini_file);
ini_file = ini2struct(x);
clear x find_ini_file

%timing info
mspline = str2num(ini_file.x_.ms0x2ep0x2eline); %ms per line

%work out how many data points to take for 1 second (i.e. 1000 ms)
windowSz = round(prefs.windowSz_ms/mspline); %in frames
% NB/ could 1 second be too long, as may be counting same cells multiple
% times?? EG Would it be better to do 250ms and then just x answer by 4?? 

%get the number of frames, the width and height of the RBCV ch
%this has been cropped around the correct part of the scan path, so
%width and height will be different
%NB height and num images should be same for loco and RBCV_binary - if
%they're not then the timing wont be the same, so you need to reextract the
%RBCV_binary and not delete any (delete in post-processing)
info = imfinfo(find_tif_file{1,n});
width = info(1).Width;
height = info(1).Height;
num_images = numel(info);
clear info;

%get loco so can look at rest periods only
loco_ch = cell2mat(findFolders(expDir, '*ch_3.tif'));
%loco channel info:
info = imfinfo(loco_ch);
width_loco = info(1).Width;
height_loco = info(1).Height;
num_images_loco = numel(info);
clear info; 

%check loco and RBCV binary channel have same amount of time in them
if height ~= height_loco || num_images ~= num_images_loco
    
    disp(['Exp dir: ', expDir]);
    disp('Check num frames in RBCflux vs loco ch... skipping...');
    
else
    
    % Load in all the frames
    % inform user as process can be quite slow
    disp('loading images for all channels...');
    for k = 1:num_images %loop frames
        %disp every 100 frames to show progress
        if (mod(k,100) == 0)
            disp(['Frame ', num2str(k), '/', num2str(num_images)])
        end
        %load the binary linescan
        A(:,:,k) = imread(find_tif_file{1,n}, k)';
        %load locomotion
        B(:,:,k) = imread(loco_ch, k)';
    end %end of frames loop
    
    
    %reorder the data:
    %put all frames on top of each other, so it becomes a 2D image
    %get the time x space line scan info
    %binary image:
    binaryLine = reshape(A, [width, height*num_images]);
    locoLine = reshape(B, [width_loco, height*num_images]);
    %dnt care about width of loco channel, just take avg for walking or not
    locoLine = nanmean(locoLine,1);
    clear A B;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loop 1 second time windows in data
    %% extract RBC flux (RBC/s), velocity and diameter over these windows
    
    %NB typical RBC/s values from cortex in literature:
    %~50 for WT (Lu et al., 2019), Li et al. (2019) ~40 for WT
    
    %preallocate nans into array outputted by loop to save memory
    RBCflux = nan(1,floor(size(binaryLine,2)/windowSz));
    locomotion = nan(1,floor(size(locoLine,2)/windowSz));
    %loop over 1 second windows...
    for a = 1:floor(size(binaryLine,2)/windowSz) %loop linescan 1s at a time
        
        %% loco
        if a == 1
            locomotion(a) = nanmean(locoLine(:,1:windowSz),2);
        else
            locomotion(a) = nanmean(locoLine(:,((a-1)*windowSz)+1 : a*windowSz),2);
        end
        
        %% RBC flux (RBC/s)
        %extract your window of 1s of binary linescan data
        if a == 1
            IM = binaryLine(:,1:windowSz);
        else
            IM = binaryLine(:,((a-1)*windowSz)+1 : a*windowSz);
        end
        
        %     % smooth image
        %     % https://uk.mathworks.com/matlabcentral/answers/375967-removal-of-shadow-area
        %     IM = imfilter(IM_ttt, fspecial('average',2),'replicate'); %5
        %     clear IM_ttt;
        %     % rebinarise the smoothed image - as it will not be binary now smoothed
        %     for b = 1:size(IM,1) %loop rows of linescan
        %         IM(b,(find(IM(b,:)>0))) = 255;
        %     end %end of loop linescan rows
        
        %invert image (to make RBCs 255 value and lumen 0 value)
        %it'll look like the lumen is black and the RBCs are white -->
        %however, in order to do the intensity profile on the shadows and
        %not the lumen --> it needs to have a value (i.e. not be 0)
        for b = 1:size(IM,1) %loop rows of linescan
            IM(b,(find(IM(b,:)==0))) = 1;
            IM(b,(find(IM(b,:)==255))) = 0;
            IM(b,(find(IM(b,:)==1))) = 255;
        end %end of loop linescan rows
        
        %take the middle few rows of the linescan and average to get the intensity
        %plot of the binary image, i.e. >0 means there is signal
        %scale is between 0 and 255 (as this is how it is binarised, i.e. if
        %all 4 lines were 255 then the average would be 255, but if only 1 line
        %is 255 and the rest are 0, then the average will be below thresh and
        %excluded, i.e. set as no RBC being present
        imIntensityCurve = nanmean(IM(round(size(IM,1)/2)-4:...
            round(size(IM,1)/2)+4,:),1);
        
        % turn off smoothing - lose info
        %save out the number of RBCs per second, by counting the baseline periods
        %(i.e. black shadows)
        % 0s in the intensity curve are black shadows where RBCs are, higher
        % numbers are white for fluorescent plasma
        % find all the zeros to find shadows, then use diff to separate out the
        % shadows which are spaced apart, those spaced apart will have diff >1
        % how many distinct shadow groups can you find (use size) - this is the
        % number of RBCs per second for this window:
        % >1 diff maybe too small, as it means one spec can affect the RBC/s number,
        % need to only count things as RB cells when there a few in a row, set at
        % 5 for now... but this is somewhat arbitrary...
        % NB if take smaller window, multiply by correct number to extrapolate
        % it to be 1s
        
        %find the blocks of 0s (i.e. shadows)
        %create a line with the set threshold
        thresholdPlot = ones(size(imIntensityCurve))*prefs.thresh;
        %find where this threshold intersects with the intensity profile - note
        %that each shadow will intersect the threshold line twice - i.e. shadow
        %on and shadow off - so will need to divide by 2
        L1(1,:) = [1:size(imIntensityCurve,2)];
        L1(2,:) = imIntensityCurve;
        L2(1,:) = [1:size(imIntensityCurve,2)];
        L2(2,:) = thresholdPlot;
        [lineIntersects] = InterX(L1,L2);
        RBCflux(a) = round(size(lineIntersects,2)/2) * 1000/prefs.windowSz_ms; clear L1 L2 lineIntersects;
        
        %     RBCflux(a) = (size(find(diff(find(imIntensityCurve == 0)) >=5 ),2)) * ...
        %         1000/prefs.windowSz_ms; %>1
        
        %check if it looks sensible what has been determined as no RBC vs RBC,
        %anything between the blue lines which is white is no RBC, black is RBC
        %present
        %plots for checking step by step what the code is doing:
        if prefs.plots %check if pref for plots is on
            if a == 1 %only plot for first loop
                % plot the image overlaid with the intensity curve as taken
                % from the middle pixels averaged
                % cnt do average across all as the RBC line is slanted
                figure;
                title(['RBCs=white, ', num2str(prefs.windowSz_ms), 'ms window, ', ...
                    num2str(RBCflux(a)/(1000/prefs.windowSz_ms)), ' nRBCs']);
                yyaxis left; xlim([1 size(IM,2)]);
                imshow(IM);
                hold on;
                yyaxis right; ylim([0 255]); xlim([1 size(IM,2)]);
                plot(imIntensityCurve, 'r', 'LineWidth', 2);
                %plot boundaries either side of threshold
%                 plot([1 size(IM,2)],[prefs.thresh-2 prefs.thresh-2],'y-');
%                 plot([1 size(IM,2)],[prefs.thresh+2 prefs.thresh+2],'y-');
                %plot the threshold line
                plot(thresholdPlot, 'b', 'LineWidth',2);
                %make fig size of screen
                set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
                    screenSz(4)]);
                %save figure into expdir
                figSave = 'RBCflux_binary_countCheck.png';
                saveas(gcf, fullfile([expDir,filesep,figSave]));
                close(gcf);
            end %end of check if first loop
        end %end of check plot prefs
        
%         % set a threshold and remove the intensity curves that are below half
%         % checked this by eye from plot above, and decided on this threshold (as
%         % came from fluorescence which were missing pixels, etc)
%         imIntensityCurve(:,find(imIntensityCurve <= prefs.thresh)) = 0;
%         
%         %can recheck intensity curve with below thres values removed
%         if false
%             figure;
%             plot(imIntensityCurve);
%             title('intensity curve final');
%         end
        
        %clear variables used repeatedly by loop but not outputted
        clear IM imIntensityCurve thresholdPlot;
        
    end %end of loop 1s windows of linescan
    
    [locomotion] = cleanLoco(locomotion);
    [locomotion] = norm_01(locomotion);
    locomotion(:,find(locomotion<=0.01))=0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% plot loco and RBC flux traces into exp dir
    
    figure;
    subplot(211);
    plot([1:size(RBCflux,2)], nanmean(RBCflux,1),'r');
    title('RBC flux');
    xlabel('Time(s)'); ylabel('RBC/s');
    subplot(212);
    plot([1:size(RBCflux,2)], locomotion,'b');
    title('Locomotion');
    xlabel('Time(s)'); ylabel('A.U.');
    %save figure into expdir
    figSave = 'RBCflux_Loco_traces.png';
    saveas(gcf, fullfile([expDir,filesep,figSave]));
    close(gcf);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save variables into expDir
    prefs2output = prefs;
    matfile = fullfile(expDir, 'contData_ls_RBCflux');
    save(matfile,'RBCflux','locomotion','prefs2output');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %end of function
