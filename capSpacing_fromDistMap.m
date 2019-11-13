clear all; close all; 

cd('D:\Dropbox (Brain Energy Lab)\Everything\Kira\Projects\HCvsV1\data_sortedByAnalysis\zstack_VascularMorphology\invivoZstacksStore');

%top dir 
fname = cd; 

%find all the groups in the top dir and put names into a cell
listing=dir(fname);
for n = 1:size(listing,1)-2
    if size(listing(2+n).name,2) <= 6
        grpNames{n} = listing(2+n).name;
    else
        grpNames{n} = [];
    end
end
%remove the empty cells 
grpNames = grpNames(~cellfun('isempty',grpNames));
clear listing; 

%loop the groups, and extract the distmap data
for n = 1:size(grpNames,2) %loop grps
    
    %inform user of progress
    disp(['processing file ',num2str(n),'/',num2str(size( ...
        grpNames,2))]);
    %disp grp name
    disp(grpNames{n});
    
    %inside grp folder look for dist maps 
    all_dirs = findFolders([fname,filesep,grpNames{n}], 'distMap*');
    
    for a = 1:size(all_dirs,2)
        
        disp(['Processing file ', num2str(a), '/', num2str(size(all_dirs,2))]);
        disp([num2str(all_dirs{1,a})]);
        
        distMap_ch = all_dirs{1,a};
        [expDir,~] = fileparts(distMap_ch);
        
        %output animal exp dir 
        grp{n}.expDir{a} = expDir;
        
        %call function to find the fps and pixel size for this recording
        [~, pxsz] = INIinfo(expDir);
        grp{n}.pxsz_um(a) = pxsz*1000000; % pixel size (um)
        clear pxsz;
        
        disp('Loading distance map...');
        %call func to load tif file with loco in
        %NB distance map gives the distance in pixels, need to convert this to
        %um
        [grp{n}.distMap{a}] = loadTifFileIn2Mat(distMap_ch);
        
        clear distMap_ch expDir;
        
        % Extract the dist map pixel values into a 2d matrix with no zero
        % values 
        %create temp var with distmap values (in um)
        ttt_um = grp{n}.distMap{a};
        %reshape this temp var to be in 2D - i.e. so get all values
        ttt_um = reshape(ttt_um,[1 size(ttt_um,1)*size(ttt_um,2)*size(ttt_um,3)]);
        %remove zeros
        ttt_um(:,find(ttt_um==0))=[];
        %put all values in order from low to high um
        ttt_um = sort(ttt_um, 'ascend'); 
        
        %store into outputted structure
        grp{n}.distMap2D{a} = ttt_um; 
        
        %in case need to check plot of data distribution 
        if false
            %Construct a histogram with a normal distribution fit.
            %and put data into 10 bins
            figure;
            histfit(ttt_um,10,'normal');
        end
 
        %find average dist(um) of furthest away 10% (i.e. most at risk)
        %if go for 5 change /10 to /20
        grp{n}.capspacing_high(a) = prctile(ttt_um, 95);
        grp{n}.capspacing_median(a) = prctile(ttt_um, 50);
        
        grp{n}.capspacing_90pcntile(a) = prctile(ttt_um, 90);
        grp{n}.capspacing_80pcntile(a) = prctile(ttt_um, 80);
        grp{n}.capspacing_75pcntile(a) = prctile(ttt_um, 75);
        grp{n}.capspacing_70pcntile(a) = prctile(ttt_um, 70);
        grp{n}.capspacing_60pcntile(a) = prctile(ttt_um, 60);

        clear ttt_um; %clear temp 2d dist map 
        
    end %end of looping dist maps
    clear all_dirs a; 
    
end %end of looping grps
clear n; 

save('capSpacing_fromDistMaps_top5%','-v7.3'); 

%get the average across z stacks for each cube (4 z stacks per region, and
%5 cubes per zstack)
zstackInd = [1:5:20,21]; 
for n = 1:size(grp,2) %loop grps
    for a = 1:size(grp{n}.capspacing_high,2) %loop distmaps loaded
        for b = 1:size(zstackInd,2)-1
            
            
            grp{n}.capspacing_high_perZstack(b) = ...
                nanmean(grp{n}.capspacing_high(:,zstackInd(b):zstackInd(b+1)-1));
            
            grp{n}.capspacing_median_perZstack(b) = ...
                nanmean(grp{n}.capspacing_median(:,zstackInd(b):zstackInd(b+1)-1));
            
            grp{n}.capspacing_90pcntile_perZstack(b) = ...
                nanmean(grp{n}.capspacing_90pcntile(:,zstackInd(b):zstackInd(b+1)-1));
            
            grp{n}.capspacing_80pcntile_perZstack(b) = ...
                nanmean(grp{n}.capspacing_80pcntile(:,zstackInd(b):zstackInd(b+1)-1));
            
            grp{n}.capspacing_75pcntile_perZstack(b) = ...
                nanmean(grp{n}.capspacing_75pcntile(:,zstackInd(b):zstackInd(b+1)-1));
            
            grp{n}.capspacing_70pcntile_perZstack(b) = ...
                nanmean(grp{n}.capspacing_70pcntile(:,zstackInd(b):zstackInd(b+1)-1));
            
            grp{n}.capspacing_60pcntile_perZstack(b) = ...
                nanmean(grp{n}.capspacing_60pcntile(:,zstackInd(b):zstackInd(b+1)-1));
            
        end
    end
end

%plot the distribution of pixels 
edges = [1:2.5:51]; %[bottom data val, step sz, top data val] - in this eg gives 20 bins

%% get number of pixels per bin and then get mean and error for this across zstacks
for n = 1:size(grp,2) %loop grps
    for a = 1:size(grp{n}.distMap2D,2) %loop distmaps loaded
        for b = 1:size(edges,2)-1
            
%             if b < size(edges,2)-1
                %extract index and data to look at
                edgeInd = [edges(b),edges(b+1)];
%             else
%                 edgeInd = [edges(b),edges(b)+(edges(2)-edges(1))]; 
%             end
            data_ttt = sort(grp{n}.distMap2D{a},'ascend'); 
            %get index for all data which falls within bin limits
            ind_ttt = find(data_ttt>=edgeInd(1) & data_ttt<edgeInd(2));
            
            
            grp{n}.numPxPerBin_all(a,b) = sum(data_ttt(ind_ttt));
            
            
        end
    end
end

%get mean number of pixels per bin for each stack
for n = 1:size(grp,2) %loop grps
    for a = 1:size(zstackInd,2)-1
        
        grp{n}.numPxPerBin_avgPerStack(a,:) = ...
            nanmean(grp{n}.numPxPerBin_all([zstackInd(a):zstackInd(a+1)-1],:),1);
        
    end
end

%get the mean px count and SEM between stacks per bin
for n = 1:size(grp,2) %loop grps
    for a = 1:size(grp{n}.numPxPerBin_avgPerStack,2)
        
        grp{n}.meanPxPerBin(a) = nanmean(grp{n}.numPxPerBin_avgPerStack(:,a),1);
        grp{n}.semPxPerBin(a) = getSEM(grp{n}.numPxPerBin_avgPerStack(:,a),1);
        grp{n}.stdPxPerBin(a) = nanstd(grp{n}.numPxPerBin_avgPerStack(:,a),1);
        
    end
end

fontSz = 30; 

binCentres = [1.25:2.5:51];

%overall plot across all stacks
figure;
hold on; 
set(gca,'LineWidth',2)

bar(binCentres,grp{1}.meanPxPerBin,'faceColor', [0.5 0 0.5], ...
    'FaceAlpha' , 0.5,'BarWidth',1,'LineWidth',1);
bar(binCentres,grp{2}.meanPxPerBin,'faceColor', [ 0.91 0.41 0.17], ...
    'FaceAlpha' , 0.5, 'BarWidth',1,'LineWidth',1);
legend({'HC','V1'},'Autoupdate','off','FontSize',25); 
legend boxoff 
 
errorbar(binCentres,grp{1}.meanPxPerBin,grp{1}.semPxPerBin,'k.', 'LineWidth',2);
errorbar(binCentres,grp{2}.meanPxPerBin,grp{2}.semPxPerBin,'k.', 'LineWidth',2);   

yLimits = get(gca,'YLim'); 

plot([nanmean(grp{1}.capspacing_median_perZstack)+2.5, ...
    nanmean(grp{1}.capspacing_median_perZstack)+2.5], ...
    yLimits, 'Color',[0.5 0 0.5], 'LineStyle','--', 'LineWidth',3);    
plot([nanmean(grp{2}.capspacing_median_perZstack)+2.5, ...
    nanmean(grp{2}.capspacing_median_perZstack)+2.5], ...
    yLimits, 'Color',[ 0.91 0.41 0.17], 'LineStyle','--', 'LineWidth',3); 
plot([nanmean(grp{1}.capspacing_high_perZstack)+2.5, ...
    nanmean(grp{1}.capspacing_high_perZstack)+2.5], ...
    yLimits, 'Color',[0.5 0 0.5], 'LineStyle','--', 'LineWidth',3);    
plot([nanmean(grp{2}.capspacing_high_perZstack)+2.5, ...
    nanmean(grp{2}.capspacing_high_perZstack)+2.5], ...
    yLimits, 'Color',[ 0.91 0.41 0.17], 'LineStyle','--', 'LineWidth',3); 

xlabel('Distance to nearest vessel (um)',...
    'FontName','Arial','FontWeight','Bold','FontSize',fontSz); 
ylabel('Number of pixels','FontName','Arial','FontWeight','Bold','FontSize',fontSz); 
ax = gca;
ax.FontSize = fontSz; 
ax.FontWeight = 'Bold'; 


fontSz = 14; 
%plot for each zstack
for n = 1:size(grp,2)
    for a = 1:size(grp{n}.numPxPerBin_avgPerStack,1)
        
        figure;
        hold on;
        set(gca,'LineWidth',2)
        
        if n == 1
            bar(binCentres,grp{n}.numPxPerBin_avgPerStack(a,:),'faceColor', [0.5 0 0.5], ...
                'FaceAlpha' , 0.5,'BarWidth',1,'LineWidth',1);
        else
            bar(binCentres,grp{n}.numPxPerBin_avgPerStack(a,:),'faceColor', [ 0.91 0.41 0.17], ...
                'FaceAlpha' , 0.5, 'BarWidth',1,'LineWidth',1);
        end
        title(['Grp=,',num2str(n),', stack=',num2str(a)]); 
        
        
        yLimits = get(gca,'YLim');
        
        plot([grp{n}.capspacing_median_perZstack(a)+2.5, ...
            grp{n}.capspacing_median_perZstack(a)+2.5], ...
            yLimits, 'k--', 'LineWidth',3);
        plot([grp{n}.capspacing_high_perZstack(a)+2.5, ...
            grp{n}.capspacing_high_perZstack(a)+2.5], ...
            yLimits, 'k--', 'LineWidth',3);

        
        xlabel('Distance to nearest vessel (um)',...
            'FontName','Arial','FontWeight','Bold','FontSize',fontSz);
        ylabel('Number of pixels','FontName','Arial','FontWeight','Bold','FontSize',fontSz);
        ax = gca;
        ax.FontSize = fontSz;
        ax.FontWeight = 'Bold';
        
    end
end









