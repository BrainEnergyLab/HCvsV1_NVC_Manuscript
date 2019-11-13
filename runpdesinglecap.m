%% Run Kira's pdecapnetwork3_1cam.m script for different brain areas and capillary spacings


% Brain regions are V1 and HC, normFlag determines whether to use the mean
% or top 5% of capillary spacings

clear all;
close all; 

Regions = [1 2]; %1 - HC, 2 - V1
Flags = [1 2]; % %1 for avg cap space, 2 for furthest cap space

% Runs in order HC, av then furthest, then V1 av then furthest


for i = 1:2
    
    brainregion = Regions(i);
    
    for j = 1%:2
        
        
        normFlag = Flags(j);
        
        [R, Vmax, capspaceav, capspacemax5pc] = pdesinglecap(brainregion, normFlag);
        
        if brainregion == 1
            
            if normFlag == 1
            Result.hippocampus.averagespacing = R;
            
            else 
            Result.hippocampus.maxspacing = R;
            
            end
            
            
       else
            if normFlag == 1
            Result.cortex.averagespacing = R;
            
            else 
            Result.cortex.maxspacing = R;
            end
            
        end
        
prompt = 'Have you closed pdetool? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end
    end
   
        prompt = 'Have you closed pdetool? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end
    
Result.capspacingav = capspaceav; % store average capillary spacings for HC (1) and V1 (2)
Result.capspacingmax5pc = capspacemax5pc; % store max 5 pc (averaged across stacks) of capillary spacings for HC (1) and V1 (2)
Result.vmaxes = Vmax; %store Vmaxes used

end
        