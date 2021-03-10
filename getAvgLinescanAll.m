
function getAvgLinescanAll(topDir)
%topDir: 'D:\Dropbox (Brain Energy
%Lab)\Everything\Kira\Projects\HCvsV1\data_sortedByAnalysis\Baseline_lsAll'

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
    
    clearvars -except a grpNames topDir grp
    
    %inform user of progress
    disp(grpNames{a});
    %label region
    grp{a}.label = grpNames{a};
    
    [find_tif_files] = findFolders([topDir,filesep,grpNames{a}], ...
        'contData_ls_RBCflux.mat');
    
    for b = 1:size(find_tif_files,2)
        
        disp([num2str(b),'/',num2str(size(find_tif_files,2))]);
        
        %find exp dir
        [expDir,~] = fileparts(find_tif_files{1,b});
        
        %label exp dir
        grp{a}.expDir{b} = find_tif_files{1,b};
        
        %load diameter
        load([expDir,filesep,'roughdiam.mat']);
        grp{a}.rough_diam(b) = roughdiam; clear roughdiam;
        
        %% load flux data:
        load(find_tif_files{1,b});
        %each frame in the data is 1 second, so time is just size of
        %loco
        time = [1:size(locomotion,2)]/(1000/prefs2output.windowSz_ms); 
        %take rest periods only - by detecting loco events and setting
        %these hct values to NaN
        [locoEvents] = findLocoEvents(locomotion,find(time>=1,1));
        for c = 1:size(locoEvents,2)
            RBCflux(:,locoEvents(1,c)+1:locoEvents(2,c))=NaN;
        end
        %now nanmean the remaining flux values to get avg hct during rest
        grp{a}.RBCflux_avg(b) = nanmean(RBCflux);
        grp{a}.RBCflux_CV(b) = nanstd(RBCflux)./nanmean(RBCflux);
        clear RBCflux time locomotion;
        
        %% load velocity
        load([expDir,filesep,'contData_ls_RBCV.mat']);
        %process loco
        [locomotion] = cleanLoco(movement);
        [locomotion] = norm_01(locomotion);
        locomotion(:,find(locomotion<=0.01))=0;
        %take rest periods only - by detecting loco events and setting
        %these hct values to NaN
        [locoEvents] = findLocoEvents(locomotion,find(time>=1,1));
        for c = 1:size(locoEvents,2)
            velocity(:,locoEvents(1,c)+1:locoEvents(2,c))=NaN;
        end
        %now nanmean the remaining vel values to get avg vel during rest
        grp{a}.RBCV_avg(b) = nanmean(velocity);
        grp{a}.RBCV_CV(b) = nanstd(velocity)./nanmean(velocity);
        clear velocity time locomotion;
        
        %% load hct
        load([expDir,filesep,'contData_ls_Hct.mat']);
        %take rest periods only - by detecting loco events and setting
        %these hct values to NaN
        [locoEvents] = findLocoEvents(locomotion,find(time>=1,1));
        for c = 1:size(locoEvents,2)
            hct(:,locoEvents(1,c)+1:locoEvents(2,c))=NaN;
        end
        %now nanmean the remaining hct values to get avg hct during rest
        grp{a}.hct_avg(b) = nanmean(hct);
        grp{a}.hct_CV(b) = nanstd(hct)./nanmean(hct);
        clear hct time locomotion;
        
        clearvars -except grp a b topDir grpNames find_tif_files;
        
    end %end of mat files within brain region grp
    
end %loop brain region grps

save('linescanMeasuresAll_rest','grp'); 

end %end of function

