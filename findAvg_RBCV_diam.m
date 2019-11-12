
function findAvg_RBCV_diam(fname) 
%
%created by Kira, June 2018
%
%Requires user to have previously extracted contData_RBCV and roughdiam 
%Looks for the average RBCV (outside of any stim periods) and the average
%vessel diameter across all sessions within the top dir - allows user to
%look for a relationship between vessel diameter and RBCV
%Note: data should be stored within a group category, e.g. V1 or APOE4, and
%then just run this code within each group separately, later code will
%compare the groups together
%
%INPUTS-
%fname = top directory which has the roughdiam.mat and contData_RBCV.mat
%files inside each individual directory
%
%OUTPUTS-
%a .mat file will be saved into the top dir which contains the average RBCV
%and diam for each session 
%a figure will be saved into the top dir which contains a scatter plot to
%allow the user to look for a pattern 
%
%Additional functions needed for this script to run:
%findFolders, findLocoEvents

%find out which group, e.g. brain region/group
[~,grpNm]=fileparts(fname); 

%get animal names, i.e. so can sort data by animal
listing=dir(fname); 
listing = {listing([listing.isdir]).name}';
listing{1,1} = []; listing{2,1} = []; 
listing = listing(~cellfun('isempty', listing));

%find all the sessions which have been through the code to find the average
%diameter of the vessel 
findSz=findFolders(fname,'contData_ls_RBCV.mat');

%% extract the data for each session
for a = 1:size(findSz,2) %loop all the sessions found 
    
    [expDir,~]=fileparts(findSz{1,a}); %get session exp dir 
    
    %inform user which exp dir
    disp(['Dir: ', expDir]);
 
    %load data
    disp('loading data...'); 
    load(fullfile(expDir, 'roughdiam.mat'));
    load(fullfile(expDir, 'contData_ls_RBCV.mat'));
    
    %provide equivalent of fps info, for line scans
    sps=[(timepts(2)-timepts(1))*mspline]/1000;
    
    %remove any stim periods (+1s to allow return to base)
    %check if it is a stim exp
    stimExpFlag=strfind(expDir,'stim');
    if stimExpFlag>1
        disp('stim experiment, removing stim period before averaging');
        %use loco finder to find stim events 
        [stimEvents]=findLocoEvents(stim,sps);
        %add 1s post stim event to allow RBCV to return to base
        for b = 1:size(stimEvents,2)
            stimEvents(2,b)=stimEvents(2,b)+round(1/sps);
        end
        %recalculate the difference between stim on and off, as updated off
        stimEvents(3,:)=stimEvents(2,:)-stimEvents(1,:); 
        %only take the data outside of these time points 
        %remove these time points, and replace with NaNs
        for b = 1:size(stimEvents,2)
            velocity(:,stimEvents(1,b):stimEvents(2,b))=NaN; 
        end
    else
        disp('not stim experiment... averaging across whole time');
    end
    
    % UPDATE
    %could also remove locomotion periods? However, this may  be
    %unnecessary for RBCV as don't get any signal during loco anyway 
    disp('removing loco periods from data');
    %use loco finder to find stim events
    [locoEvents]=findLocoEvents(locomotion,sps);
    
    if ~isempty(locoEvents)
        %add 1s post stim event to allow RBCV to return to base
        for b = 1:size(locoEvents,2)
            locoEvents(2,b)=locoEvents(2,b)+round(1/sps);
        end
        %recalculate the difference between stim on and off, as updated off
        locoEvents(3,:)=locoEvents(2,:)-locoEvents(1,:);
        %only take the data outside of these time points
        %remove these time points, and replace with NaNs
        for b = 1:size(locoEvents,2)
            velocity(:,locoEvents(1,b):locoEvents(2,b))=NaN;
        end
    end
    
    %rough diam of the vessel
    mean_sz(a)=roughdiam;
    %average and std for velocity outside stim periods
    mean_vel(a)=nanmean(velocity);
    std_vel(a)=nanstd(velocity);
    
    %find animal name so can categorise each data
    for b=1:size(listing,1)
        ttt = strfind(expDir,listing{b,1});
        if ttt>1
            anInd(b)=1;
        else
            anInd(b)=0;
        end
    end
    anLabel(a)=find(anInd);

    %clear the variables that will be extracted, to ensure removed from
    %workspace
    clear roughdiam velocity; 
    
end %end of looping the sessions 

%save the data in top dir
matfile = fullfile(fname,'RBCV_meanstd.mat');
save(matfile,'anLabel','listing','mean_sz','mean_vel','std_vel','grpNm'); 

%% plot the data 

%find how many animals (N) in sample 
nAnimals=max(anLabel); 
%max 7 colours in the colorstring - if N is higher than this, loop back to
%beginning colours 
colorstring={'b','k','m','c','y','g','r'}; 

%check if the loop for plotting is large enough to contain ALL the animals
%in the sample 
%NOTE inefficient way to do this!!
if nAnimals > 3*size(colorstring,2) 
   disp('WARNING: some of the data is not being plotted, as N too high');
   disp('will need to update code to include');
   disp(['code currently incorporates N=',num2str(3*size(colorstring,2))]);
end

figure;
%find size of screen (for figure sizing)
screenSz=get(0,'Screensize');
%make fig size of screen
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);
hold on; 
for a = 1:size(anLabel,2) %loop sessions
        if anLabel(a) <= size(colorstring,2) %loop animal ID
            h(anLabel(a))=errorbar(mean_sz(a),mean_vel(a),std_vel(a), ...
                'Color', colorstring{anLabel(a)}, 'Marker','*');
        elseif anLabel(a) > size(colorstring,2) && anLabel(a) < ...
                2*(size(colorstring,2))
            h(anLabel(a))=errorbar(mean_sz(a),mean_vel(a),std_vel(a), ...
                'Color', colorstring{anLabel(a)-size(colorstring,2)}, ...
                'Marker','o');
        elseif anLabel(a) > 2*(size(colorstring,2)) && anLabel(a) < ...
                3*(size(colorstring,2))
            h(anLabel(a))=errorbar(mean_sz(a),mean_vel(a),std_vel(a), ...
                'Color', colorstring{anLabel(a)-2*(size(colorstring,2))},...
                'Marker','o');
        end %end of looping animal ID 
end %end of looping sessions 
legend([h(:)],listing{:})
xlabel('vessel diameter, um');
ylabel('RBCV, mm/s');
title([num2str(grpNm), ', nAnimals=', num2str(nAnimals), ...
    ', nVessels=', num2str(size(anLabel,2))]);
%save plot into exp_dir as png
saveas(gcf, fullfile(expDir, ['RBCV_diam_scatter_all_',num2str(grpNm),...
    '.png']));
close; %close figure


end %end of function

