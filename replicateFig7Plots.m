
cd('D:\Dropbox (Brain Energy Lab)\Everything\ForCatherine\NVCinHCvsV1_paperPlan\Diffusionmodelling\Kira\');
load('Result021019_50th_95th.mat')

%make time vector:
tlist = [0:0.0001:.1]; %[0:0.0001:.1]; %tspan (minutes)
%specify Km
Km = 0.001; %(=1uM)

%% HC normal cap spacing: 

%extract the line between cap and edge from which to get o2 concs
xInd = find((Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(1,:)>=0) & ...
    Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(1,:)<=Result.capspacingav(1));
yInd = find((Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(2,:)>-0.00005) & ...
    Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(2,:)<=0.00005);
indOfInterest = intersect(xInd,yInd);
%plot arena with line to check looks ok
figure;
title('HC arena, with ROI between cap and edge'); 
scatter(Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(1,:), ...
    Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(2,:), 'k.');
hold on;
scatter(Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(1,indOfInterest),...
    Result.hippocampus.averagespacing{1, 1}.Mesh.Nodes(2,indOfInterest), 'x');

for n = 1:size(Result.vmaxes,2)

%extract the o2 conc along this line: 
HC_avg_o2conc(n,:,:) = Result.hippocampus.averagespacing{1, n}. ...
    Temperature(indOfInterest,:);

%extract the dist (um) along this line:
%i.e. take values from x coordinates and multiply by 1000 to get into um
HC_avg_dist(n,:) = Result.hippocampus.averagespacing{1, n}.Mesh.Nodes(1,indOfInterest)*10000; 
%for some reason x coords are not in order?? so need to reorder in
%ascending
[HC_avg_dist(n,:),ind_ttt] = sort(squeeze(HC_avg_dist(n,:)),'ascend');
%then take this index and reorder the o2 conc, as needs to match cap dist
HC_avg_o2conc(n,:,:) = HC_avg_o2conc(n,ind_ttt,:);

%extract vo2 averaged across time, so get a vo2 val per space
%between caps
for o = 1:size(HC_avg_o2conc,2)
    HC_avg_vo2(n,o,:) = (Result.vmaxes(n) * mean(HC_avg_o2conc(n,o,:),3))/ ...
        (mean(HC_avg_o2conc(n,o,:),3)+Km);
end

clear o ind_ttt; 

end
clear n xInd yInd indOfInterest; 


%% HC max cap spacing: 

%extract the line between cap and edge from which to get o2 concs
xInd = find((Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(1,:)>=0) & ...
    Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(1,:)<=Result.capspacingmax5pc(1));
yInd = find((Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(2,:)>-0.00005) & ...
    Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(2,:)<=0.00005);
indOfInterest = intersect(xInd,yInd);
%plot arena with line to check looks ok
figure;
title('HC arena, with ROI between cap and edge'); 
scatter(Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(1,:), ...
    Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(2,:), 'k.');
hold on;
scatter(Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(1,indOfInterest),...
    Result.hippocampus.maxspacing{1, 1}.Mesh.Nodes(2,indOfInterest), 'x');

for n = 1:size(Result.vmaxes,2)

%extract the o2 conc along this line: 
HC_max_o2conc(n,:,:) = Result.hippocampus.maxspacing{1, n}. ...
    Temperature(indOfInterest,:);

%extract the dist (um) along this line:
%i.e. take values from x coordinates and multiply by 1000 to get into um
HC_max_dist(n,:) = Result.hippocampus.maxspacing{1, n}.Mesh.Nodes(1,indOfInterest)*10000; 
%for some reason x coords are not in order?? so need to reorder in
%ascending
[HC_max_dist(n,:),ind_ttt] = sort(squeeze(HC_max_dist(n,:)),'ascend');
%then take this index and reorder the o2 conc, as needs to match cap dist
HC_max_o2conc(n,:,:) = HC_max_o2conc(n,ind_ttt,:);

%extract vo2 maxd across time, so get a vo2 val per space
%between caps
for o = 1:size(HC_max_o2conc,2)
    HC_max_vo2(n,o,:) = (Result.vmaxes(n) * mean(HC_max_o2conc(n,o,:),3))/ ...
        (mean(HC_max_o2conc(n,o,:),3)+Km);
end

clear o ind_ttt; 

end
clear n xInd yInd indOfInterest;  

%% V1 normal cap spacing: 

%extract the line between cap and edge from which to get o2 concs
xInd = find((Result.cortex.averagespacing{1, 1}.Mesh.Nodes(1,:)>=0) & ...
    Result.cortex.averagespacing{1, 1}.Mesh.Nodes(1,:)<=Result.capspacingav(2));
yInd = find((Result.cortex.averagespacing{1, 1}.Mesh.Nodes(2,:)>-0.00005) & ...
    Result.cortex.averagespacing{1, 1}.Mesh.Nodes(2,:)<=0.00005);
indOfInterest = intersect(xInd,yInd);
%plot arena with line to check looks ok
figure;
title('V1 arena, with ROI between cap and edge'); 
scatter(Result.cortex.averagespacing{1, 1}.Mesh.Nodes(1,:), ...
    Result.cortex.averagespacing{1, 1}.Mesh.Nodes(2,:), 'k.');
hold on;
scatter(Result.cortex.averagespacing{1, 1}.Mesh.Nodes(1,indOfInterest),...
    Result.cortex.averagespacing{1, 1}.Mesh.Nodes(2,indOfInterest), 'x');

for n = 1:size(Result.vmaxes,2)

%extract the o2 conc along this line: 
V1_avg_o2conc(n,:,:) = Result.cortex.averagespacing{1, n}. ...
    Temperature(indOfInterest,:);

%extract the dist (um) along this line:
%i.e. take values from x coordinates and multiply by 1000 to get into um
V1_avg_dist(n,:) = Result.cortex.averagespacing{1, n}.Mesh.Nodes(1,indOfInterest)*10000; 
%for some reason x coords are not in order?? so need to reorder in
%ascending
[V1_avg_dist(n,:),ind_ttt] = sort(squeeze(V1_avg_dist(n,:)),'ascend');
%then take this index and reorder the o2 conc, as needs to match cap dist
V1_avg_o2conc(n,:,:) = V1_avg_o2conc(n,ind_ttt,:);

%extract vo2 averaged across time, so get a vo2 val per space
%between caps
for o = 1:size(V1_avg_o2conc,2)
    V1_avg_vo2(n,o,:) = (Result.vmaxes(n) * mean(V1_avg_o2conc(n,o,:),3))/ ...
        (mean(V1_avg_o2conc(n,o,:),3)+Km);
end

clear o ind_ttt; 

end
clear n xInd yInd indOfInterest; 


%% v1 max cap spacing: 

%extract the line between cap and edge from which to get o2 concs
xInd = find((Result.cortex.maxspacing{1, 1}.Mesh.Nodes(1,:)>=0) & ...
    Result.cortex.maxspacing{1, 1}.Mesh.Nodes(1,:)<=Result.capspacingmax5pc(2));
yInd = find((Result.cortex.maxspacing{1, 1}.Mesh.Nodes(2,:)>-0.00005) & ...
    Result.cortex.maxspacing{1, 1}.Mesh.Nodes(2,:)<=0.00005);
indOfInterest = intersect(xInd,yInd);
%plot arena with line to check looks ok
figure;
title('V1 arena, with ROI between cap and edge'); 
scatter(Result.cortex.maxspacing{1, 1}.Mesh.Nodes(1,:), ...
    Result.cortex.maxspacing{1, 1}.Mesh.Nodes(2,:), 'k.');
hold on;
scatter(Result.cortex.maxspacing{1, 1}.Mesh.Nodes(1,indOfInterest),...
    Result.cortex.maxspacing{1, 1}.Mesh.Nodes(2,indOfInterest), 'x');

for n = 1:size(Result.vmaxes,2)

%extract the o2 conc along this line: 
V1_max_o2conc(n,:,:) = Result.cortex.maxspacing{1, n}. ...
    Temperature(indOfInterest,:);

%extract the dist (um) along this line:
%i.e. take values from x coordinates and multiply by 1000 to get into um
V1_max_dist(n,:) = Result.cortex.maxspacing{1, n}.Mesh.Nodes(1,indOfInterest)*10000; 
%for some reason x coords are not in order?? so need to reorder in
%ascending
[V1_max_dist(n,:),ind_ttt] = sort(squeeze(V1_max_dist(n,:)),'ascend');
%then take this index and reorder the o2 conc, as needs to match cap dist
V1_max_o2conc(n,:,:) = V1_max_o2conc(n,ind_ttt,:);

%extract vo2 maxd across time, so get a vo2 val per space
%between caps
for o = 1:size(V1_max_o2conc,2)
    V1_max_vo2(n,o,:) = (Result.vmaxes(n) * mean(V1_max_o2conc(n,o,:),3))/ ...
        (mean(V1_max_o2conc(n,o,:),3)+Km);
end

clear o ind_ttt; 

end
clear n xInd yInd indOfInterest;  


%time series o2 conc
squeeze(nanmean(HC_avg_o2conc(1,:,:),2))*1000;

%dist and o2 conc
squeeze(nanmean(HC_avg_o2conc(7,:,:),3))*1000; 

%vo2 cap dist
HC_avg_vo2(1,:)/0.5;
ans';

%vmax vs o2 1/2way
HC_avg_o2conc(:,52,size(tlist,2))*1000;
HC_max_o2conc(:,43,size(tlist,2))*1000;
V1_avg_o2conc(:,36,size(tlist,2))*1000;
V1_max_o2conc(:,30,size(tlist,2))*1000;

%vmax vs vo2
% vmax = [0.5,1,1.5,2];
vmax = [2];
% vmax = [2];
for o = 1:size(HC_avg_vo2,1)
    ttt(o) = HC_avg_vo2(o,52)/vmax(o)
end
% ttt';
for o = 1:size(HC_max_vo2,1)
    ttt(o) = HC_max_vo2(o,43)/vmax(o)
end
% ttt';
for o = 1:size(V1_avg_vo2,1)
    ttt(o) = V1_avg_vo2(o,36)/vmax(o)
end
% ttt';
for o = 1:size(V1_max_vo2,1)
    ttt(o) = V1_max_vo2(o,30)/vmax(o)
end
% ttt';


figure; 
hold on; 
for a = 1:9
   plot(tlist,squeeze(HC_avg_o2conc(a,9,:))*1000); 
end

figure; 
hold on; 
for a = 1:9
    
   plot(tlist,squeeze(HC_max_o2conc(a,9,:))*1000); 
end




    