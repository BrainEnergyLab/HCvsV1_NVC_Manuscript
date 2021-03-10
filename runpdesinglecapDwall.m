%% Run Kira's pdecapnetwork3_1cam.m script for different brain areas and capillary spacings


% Brain regions are V1 and HC, normFlag determines whether to use the mean
% or top 5% of capillary spacings

clear all;
close all; 

filenamepath = uigetfile('*.csv');
inputdata = readtable(filenamepath);

medspacing = table2array(inputdata(:,3));
pc95spacing = table2array(inputdata(:,4));

thecapspacings = cat(1,medspacing, pc95spacing)./10000; %convert to cm too

Ncapspacings = length(thecapspacings);

Oxygenstart = [1 2]; %1 - HCOxy, 2 - V1Oxy
Flags = 1:Ncapspacings; % %1 for avg cap space, 2 for furthest cap space
O2cap = [0.014 0.021]; % O2 at cap in mM - 14um is HC estimate, 21um is V1 estimate

% Runs in order HC, av then furthest, then V1 av then furthest


for k = 1:length(O2cap)
    close all
    clear Result
    O2level = O2cap(k);
    kval = num2str(k);
    saveTitle = ['modelResults_varyO2t_' kval '.mat']; 
    
    

for i = 1:2
    
    brainregion = Oxygenstart(i);
    
    for j = 1:Ncapspacings
        
        
        capspacing = thecapspacings(j);
        
        [R, Vmax] = pdesinglecapDwall(O2level, capspacing);
        
        if brainregion == 1
            
            
            Result.hippocampus(:,j) = R;
        
            
            
       else
            
            Result.cortex(:,j) = R;
            
            
            
        end
        
        close all
        
% prompt = 'Have you closed pdetool? Y/N [Y]: ';
% str = input(prompt,'s');
% if isempty(str)
%     str = 'Y';
% end
     end
%    
%         prompt = 'Have you closed pdetool? Y/N [Y]: ';
% str = input(prompt,'s');
% if isempty(str)
%     str = 'Y';
% end
    
Result.capspacings = thecapspacings; % store average capillary spacings for HC (1) and V1 (2)
Result.oxygenlevels = O2level; % store max 5 pc (averaged across stacks) of capillary spacings for HC (1) and V1 (2)
Result.vmaxes = Vmax; %store Vmaxes used
Result.inputtable = inputdata;
save(saveTitle, 'Result', '-v7.3')

end

%O2result(k) = Result;;

end

; 
        