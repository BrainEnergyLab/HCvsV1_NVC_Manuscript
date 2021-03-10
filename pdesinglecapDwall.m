

function [R, Vmax] = pdesinglecapDwall(O2level, capspacings)

%% select region:

% brainregion = 1; %1 - HC, 2 - V1
% normFlag = 1; %1 for avg cap space, 2 for furthest cap space

% capspaceav = [0.00144 0.00103]; %now this is 75th centile not median (from distance map, add 2.5 so from middle of 5um capillary)
% capspacemax5pc = [0.00295 0.002164]; % and won't run as Flags = 1....80th centile not 95th centile


%% set model param 
Vmax = 2; %[0.5 1 1.5 2 2.5 3]; %[1 1.5 2 2.5 3];  % This is a range to cover the Sakzadzic example and the range reported by Gould et al (Kleinfeld and Tsai co authors) (2017)
O2conc = O2level; % HC = 14; V1 = 21 uM;  %in mM - calculated from sO2 to pO2(RBC) using Hill curve as per Lyons and Parpa
Km = 0.001; %(=1uM)
Dtissue = 9.25e-4; %cm2/min O2 Diffusion coeefficient in tissue
Dwall = 0.4*Dtissue; % to reflect up to 40% smaller diffusion coefficient in vessel wall than tissue https://www.ncbi.nlm.nih.gov/m/pubmed/24027235/
numcaps = 1; 
% if normFlag == 1
%     capspacings = capspaceav; % %spacing between caps from 3D distance map (cm) -average spacing
%    
% else
%     capspacings = capspacemax5pc;% max spacing
% end
%capillary radius (cm)
capradius = 0.00025; %2.5e-4cm
outercapradius = 0.00027; %2.7e-4 cm% plus 200 nm for BM i.e. radius goes from 2.5um to 2.7um

%%  create simple arena with 1 cap:
%coords of outer arena:
coords(:,1) = [1,0,0,capspacings]; %1 for circle,x coord, ycoord, radius
%name of outer arena
name{1,:} = 'C1';
%coords of cap1:
%NB radius of cap is 0.00025 for both v1 and hc
coords(:,2) = [1,0,0,capradius]; %1 for circle,x coord, ycoord, radius
%vessel outer wall coordinate
coords(:,3) = [1,1e-9,1e-9,outercapradius]; % for outer edge of capillary

%define capillary wall
name{2,:} = 'lumen';
name{3,:} = 'ablumen';
% gd(1,:) = name(1,:);
% catnames = [name(1,:) '-('];
for n = 1:3
 
% 
% xc = coords(2,:);
% yc = coords(3,:);
% r = coords(4,:); 
% %thename = capilname;


% for n = 1:3
%     pdecirc(xc(n),yc(n),r(n),name{n})
%     
% end
%the default name for the circles will be C(n).
%    if n < numcaps
%        catnames = [catnames thename '+'];
%    else
%        catnames = [catnames thename ')'];
%    end
% end


%sf = catnames;
sf = ['C1' '-' 'lumen'];
ns = name';
% generate geometry for arena
[dl,bt] = decsg(coords,sf,ns); 

%% use arena to get diffusion model:
Diffusionmodel = createpde('thermal','transient');

pg = geometryFromEdges(Diffusionmodel,dl);
pdegplot(Diffusionmodel,'EdgeLabels','on', 'CellLabels', 'on', 'FaceLabels', 'on')
axis equal
xlim('auto')

msh = generateMesh(Diffusionmodel, 'Hmax', 0.0001); %was 0.0001
figure
pdeplot(Diffusionmodel)
axis equal
title 'Model with Finite Element Mesh Displayed'

capconc = O2conc; %(brainregion); % in mmolar
thermalBC(Diffusionmodel,'Edge', [5,6,7,8],'Temperature',capconc);
thermalBC(Diffusionmodel,'Edge',[1,2,3,4], 'HeatFlux',0);
thermalProperties(Diffusionmodel, 'Face', 1, 'ThermalConductivity',Dwall,...
                                'MassDensity',1,...
                                'SpecificHeat',1);   
thermalProperties(Diffusionmodel, 'Face', 2, 'ThermalConductivity',Dtissue,...
                                'MassDensity',1,...
                                'SpecificHeat',1);   
tlist = 0:0.001:.1; %tspan (changed from 0.001)
thermalIC(Diffusionmodel,0);

for a = 1:size(Vmax,2)
    
    disp([num2str(a),'/',num2str(size(Vmax,2))]); 
    
    %internal heat source won't work with 1 cap??? 
    qFunc = @(region,state) findsink(region,state,Vmax(a),Km);
    %Specficy internal heat source
    internalHeatSource(Diffusionmodel,qFunc);

    R{a} = solve(Diffusionmodel,tlist);
    
    figure;
    axis equal
    pdeplot(Diffusionmodel,'XYData',R{a}.Temperature(:,end),'Contour','on','ColorMap','hot');
    axis equal
    title(['Vmax=',num2str(Vmax(a)),', Concentration, Final Time, Transient Solution']);

end

% if brainregion == 1
%     HC = R;
%     clear R;
% else
%     V1 = R;
%     clear R;
% end
          
end




