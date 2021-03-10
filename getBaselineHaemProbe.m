
function getBaselineHaemProbe(topDir)

clear grp; 

% WRITTEN APRIL 2020 - improved code to get out mean and STD over time of
% haemodynamic measures during rest 

%time in seconds to remove before and after loco event to allow for
%transition time, i.e. so definitely getting rest periods 
prefs.transitionTime = 1; %seconds 

%find all the subfolders - i.e. groups - in the top dir
listing=dir(topDir);
for a = 1:size(listing,1)-2
    if size(listing(2+a).name,2) <= 5
        grpNames{a} = listing(2+a).name;
    else
        grpNames{a} = [];
    end
end
%remove the empty cells 
grpNames = grpNames(~cellfun('isempty',grpNames));

%loop the groups, extract the data
for a = 1:size(grpNames,2)
    
    clearvars -except a grpNames topDir grp prefs
    
    %inform user of progress
    disp(grpNames{a});

    [find_mat_files] = findFolders([topDir,filesep,grpNames{a}], ...
        'contData.mat');
    
    for b = 1:size(find_mat_files,2)
        
        [expDir,~] = fileparts(find_mat_files{1,b}); 
        
        %load data
        %NB/ order of the haem data:
        %1-flux, 2-speed, 3-so2, 4-hbo, 5-hbr, 6-hbt, 7-cmro2
        clear data movement fps time frames; 
        load(find_mat_files{1,b});
        
        %clean loco
        clear locomotion; 
        [locomotion] = cleanLoco(movement);
        [locomotion] = norm_01(locomotion);
        locomotion(:,find(locomotion<=0.01))=0;
        
        %take rest periods only - by detecting loco events, 
        %finding 1s before and after to account for transition time and 
        %setting these haem values to NaNs
        clear locoEvents; 
        [locoEvents] = findLocoEvents(locomotion,fps);
        for c = 1:size(locoEvents,2)
            if locoEvents(1,c)-(round(fps*prefs.transitionTime)) >= 1 && ...
                    locoEvents(2,c)+(round(fps*prefs.transitionTime)) < size(data,2)
                data(:,locoEvents(1,c)-(round(fps*prefs.transitionTime))...
                    :locoEvents(2,c)+(round(fps*prefs.transitionTime)))=NaN;
            elseif locoEvents(1,c)-round(fps*prefs.transitionTime) >= 1 && ...
                    locoEvents(2,c)+round(fps*prefs.transitionTime) >= size(data,2)
                data(:,locoEvents(1,c)-(round(fps*prefs.transitionTime)):end)=NaN;
            else
                data(:,1:locoEvents(2,c)+round(fps*prefs.transitionTime))=NaN;
            end
        end
        
        %now nanmean the remaining hct values to get avg hct during rest
        for c = 1:size(data,1) %loop 7 haem variables
            grp{a}.data_avg(b,c) = nanmean(data(c,:));
            %coefficient of variation - i.e. relative SD, so more
            %comparable between two groups with v diff means
            grp{a}.data_cv(b,c) = (nanstd(data(c,:))./nanmean(data(c,:)))*100;
            grp{a}.expDir{b} = find_mat_files{b};
            %get animal ID out for each recording 
            grp{a}.anLabel{b} = extractBefore(extractAfter ...
                (grp{a}.expDir{b},109),'\');
            if b == 1
                grp{a}.anID(b) = 1;
            elseif strcmp(grp{a}.anLabel(b-1), grp{a}.anLabel(b))
                grp{a}.anID(b) =  grp{a}.anID(b-1);
            else
                grp{a}.anID(b) =  grp{a}.anID(b-1)+1;
            end
        end %end of loop 7 haem variables 

    end %end of mat files within brain region grp
    
end %loop brain region grps

% get out average mean over time and average std over time per animal
for a = 1:size(grp,2)
    for b = 1:max(grp{a}.anID)
        if size(find(grp{a}.anID==b),2)>1
            grp{a}.dataAvg_anAvg(b,:) = nanmean(grp{a}.data_avg ...
                (find(grp{a}.anID==b),:));
            grp{a}.dataCV_anAvg(b,:) = nanmean(grp{a}.data_cv ...
                (find(grp{a}.anID==b),:));
        else
            grp{a}.dataAvg_anAvg(b,:) = grp{a}.data_avg ...
                (find(grp{a}.anID==b),:);
            grp{a}.dataCV_anAvg(b,:) = grp{a}.data_cv ...
                (find(grp{a}.anID==b),:);
        end
    end
end

save('AvgHaem_allSess.mat','grp');

%% plots

haemLabels = {'Flux','Speed','SO2','HbO','Hbr','Hbt','CMRO2'};
unitLabels = {'A.U.','A.U.','%','A.U.','A.U.','A.U.','A.U.'};

figure;
for a = 1 : 7 %loop 7 haem vars 
    if a == 1
        counter = 0;
    else
        counter = counter + 1;
    end 
    
    %plot MEAN over time averaged for each animal:
    %extract data needed for this plot 
    clear x1 y1 
    x1 = [ones(size(grp{1}.dataAvg_anAvg(:,a))); ...
        ones(size(grp{2}.dataAvg_anAvg(:,a)))*2];
    y1 = [grp{1}.dataAvg_anAvg(:,a); grp{2}.dataAvg_anAvg(:,a)];
    %plot mean over time averaged for each animal
    subplot(7,2,a+counter);
    [~,p]=ttest2(y1(find(x1==1)),y1(find(x1==2))); 
    title([haemLabels{a},' Avg, p=', num2str(p)]); clear p; 
    xlim([0 3]);
    hold on;
    ax1 = bar(1,nanmean(y1(find(x1==1),:)));
    ax2 = bar(2,nanmean(y1(find(x1==2),:)));
    scatter(x1,y1,'ko', 'LineWidth',1.5);
    ax1.FaceColor = [.5 0 .5]; %purple
    ax2.FaceColor = [0.91 0.41 0.17]; %orange
    xlabel('Region'); ylabel(unitLabels{a});
    set(gca, 'XTick', []);
    if a == 1
        legend({'HC','V1'});
    end
    
    %plot coefficient of var over time averaged for each animal:
    %extract data needed for this plot - mean 
    clear x1 y1 
    x1 = [ones(size(grp{1}.dataCV_anAvg(:,a))); ...
        ones(size(grp{2}.dataCV_anAvg(:,a)))*2];
    y1 = [grp{1}.dataCV_anAvg(:,a); grp{2}.dataCV_anAvg(:,a)];
    %plot mean over time averaged for each animal
    subplot(7,2,a*2);
    [~,p]=ttest2(y1(find(x1==1)),y1(find(x1==2))); 
    title([haemLabels{a},' CV, p=', num2str(p)]); clear p; 
    xlim([0 3]);
    hold on;
    ax1 = bar(1,nanmean(y1(find(x1==1),:)));
    ax2 = bar(2,nanmean(y1(find(x1==2),:)));
    scatter(x1,y1,'ko', 'LineWidth',1.5);
    ax1.FaceColor = [.5 0 .5]; %purple
    ax2.FaceColor = [0.91 0.41 0.17]; %orange
    xlabel('Region'); ylabel(unitLabels{a});
    set(gca, 'XTick', []);
    
end
saveas(gcf, 'allHaemData_Mean_CV_avgPerAn.png');


end %end of function 

