function avgRestAndLocoData_probe(fname, prefs)
%
%function written by Kira, May 2018
%requires user to have previously run extractRestAndLocoData_probe
%searches for outputted mat file: Avg4RestPeriods_allTrials_allSess
%
%send in topDir, i.e. with multiple groups which are being compared (e.g.
%HC vs V1, or APOE3 vs APOE4) - this will compare the groups
%
%INPUTS -
%fname : top dir with multiple groups in to be compared, each dir should
%have mat file in
%prefs : plot prefs - if want bar charts plotted to compare between groups
%- this is at the bottom of the code, and is currently hard coded for
%Kira's HC vs V1 analysis, so may need to update or turn plots off
%
%OUTPUTS-
%nothing outputted into workspace, but will save data and figures into top
%dir
%
%functions needed to run this code:
%findFolders

if nargin < 2
    %plot bar charts to display comparison between groups
    prefs.plot = 1;
end

%look for all data files to be analysed
finddatafile = findFolders(fname, ...
    'Avg4RestnLocoPeriods_allTrials_allSess.mat');

%find all the subfolders - i.e. groups - in the top dir
listing=dir(fname);
for a = 1:size(listing,1)-2
    if size(listing(2+a).name,2) <= 6
        grpNames{a} = listing(2+a).name;
    else
        grpNames{a} = [];
    end
end
%remove the empty cells 
grpNames = grpNames(~cellfun('isempty',grpNames));

%loop the groups, extract the data
for a = 1:size(grpNames,2)

    %inform user of progress
    disp(['processing file ',num2str(a),'/',num2str(size( ...
        finddatafile,2))]);
    disp(grpNames{a}); 
    
    %load the data
    load(fullfile(fname,filesep,grpNames{a},filesep, ...
        'Avg4RestnLocoPeriods_allTrials_allSess.mat'));
    
    %REMOVE EMPTY CELLS:
    %see if there are any empty cells in loco data:
    emptyInd=cellfun('isempty', dataLocoTrials);
    %can tell by summing index, any empties are labelled 1
    if sum(emptyInd) > 0
        for b = 1:size(emptyInd,2)
            if emptyInd(b)==1
                dataLocoTrials{b}=[];
                
            end
        end
    end
    %remove empty cells from loco var
    dataLocoTrials(cellfun('isempty', dataLocoTrials)) = [];
    %see if there are any empty cells in rest data:
    emptyInd=cellfun('isempty', dataRestTrials);
    %can tell by summing index, any empties are labelled 1
    if sum(emptyInd) > 0
        for b = 1:size(emptyInd,2)
            if emptyInd(b)==1
                dataRestTrials{b}=[];
                
            end
        end
    end
    %remove empty cells from rest var
    dataRestTrials(cellfun('isempty', dataRestTrials)) = [];
    
    %can also put group label into the var - need any folders other than
    %group labels to come after grp labels
    dataAvg{a}.label = grpNames{a};
    
    %get an average for each session (which had data in)
    %dim1 is the diff haem param, dim2 is the data during rest/loco period
    %dataAvg - each cell is diff group; data inside dataAvg is stored
    %sessionsxhaemParam
    %rest:
    for b = 1:size(dataRestTrials,2) %loop sessions with rest data
        for c = 1:size(dataRestTrials{1},1) %loop haem param
            %average across rest period
            dataAvg{a}.rest(b,c)=nanmean(dataRestTrials{1,b}(c,:),2);
        end %end of loop haem param
    end %end of loop sessions
    %loco:
    for b = 1:size(dataLocoTrials,2) %loop sessions with loco data
        for c = 1:size(dataLocoTrials{1},1) %loop haem param
            %average across loco period (dim2)
            dataAvg{a}.loco(b,c)=nanmean(dataLocoTrials{1,b}(c,:),2);
        end
    end %end of loop sessions
    
    %state number of trials
    dataAvg{a}.nTrials_rest = size(dataAvg{a}.rest,1);
    dataAvg{a}.nTrials_loco = size(dataAvg{a}.loco,1);
    
    %take mean, std, and sem across the sessions:
    %NB/ order of the haem data:
    %1-flux, 2-speed, 3-so2, 4-hbo, 5-hbr, 6-hbt, 7-cmro2
    for b = 1:size(dataLocoTrials{1},1) %loop haem param
        %rest mean:
        dataAvg{a}.rest_mean(:,b)=nanmean(dataAvg{a}.rest(:,b));
        %rest sem:
        [dataAvg{a}.rest_sem(:,b)]=getSEM(dataAvg{a}.rest(:,b));
        %rest std:
        dataAvg{a}.rest_std(:,b)=nanstd(dataAvg{a}.rest(:,b));
        %loco mean:
        dataAvg{a}.loco_mean(:,b)=nanmean(dataAvg{a}.loco(:,b));
        %loco sem:
        [dataAvg{a}.loco_sem(:,b)]=getSEM(dataAvg{a}.loco(:,b));
        %loco std:
        dataAvg{a}.loco_std(:,b)=nanstd(dataAvg{a}.loco(:,b));
    end %end of loop haem param
    
end %end of looping folders to extract data

%save the data out into top dir- can be input into SPSS for stats tests
matfile = fullfile(fname, 'Avg4RestnLocoPeriods_allSess.mat');
save(matfile, 'dataAvg','-v7.3');

%% plots to save

if prefs.plot == 1
    
    %which haem parameters to plot (in preferred plotting order)
    data2plot = [7; 1; 3];
    dataLabel = {'cmro2'; 'flux'; 'so2'};
    colorstring1 = {'g','m','k'}';
    colorstring2 = {'g.','m.','k.'}';
    yLabelNm = {'A.U.'; 'A.U.'; '%'};
    
    figure;
    screenSz=get(0,'Screensize');
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
        screenSz(4)/2]);
    for a = 1:size(data2plot,1) %loop data to plot
        subplot(1,3,a);
        bar(1:4,[dataAvg{1}.rest_mean(data2plot(a)), ...
            dataAvg{1}.loco_mean(data2plot(a)), ...
            dataAvg{2}.rest_mean(data2plot(a)), ...
            dataAvg{2}.loco_mean(data2plot(a))], colorstring1{a});
        hold on;
        errorbar(1:4,[dataAvg{1}.rest_mean(data2plot(a)), ...
            dataAvg{1}.loco_mean(data2plot(a)), ...
            dataAvg{2}.rest_mean(data2plot(a)), ...
            dataAvg{2}.loco_mean(data2plot(a))], ...
            [dataAvg{1}.rest_sem(data2plot(a)), ...
            dataAvg{1}.loco_sem(data2plot(a)), ...
            dataAvg{2}.rest_sem(data2plot(a)), ...
            dataAvg{2}.loco_sem(data2plot(a))], colorstring2{a});
        xlim([0 5]);
        title(dataLabel(a));
        %update to use grp label insteadof hardcode region
        set(gca, 'XTickLabel', {[grpNames{1},' rest, n',...
            num2str(dataAvg{1}.nTrials_rest)], ...
            [grpNames{1},' loco, n', num2str(dataAvg{1}.nTrials_loco)], ...
            [grpNames{2},' rest, n', num2str(dataAvg{2}.nTrials_rest)], ...
            [grpNames{2},' loco, n', num2str(dataAvg{2}.nTrials_loco)]});
        ylabel(yLabelNm(a));
    end %end of loop data to plot
    %save plot into exp_dir as png
    saveas(gcf, fullfile(fname, 'barchart_haemParam_compGrps.png'));
    close; %close figure
    
end

end %end of function