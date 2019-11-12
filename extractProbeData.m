function extractProbeData(fname)

%function written by Kira, August 2017, updated May 2018
%
%function searches through experimental directories for excel files
%containing oxyprobe/LDF data 
%it will then extract the data as a .mat file, and save into the directory
%
%INPUTS-
%fname = directory which contaisn excel sheets 
%NB/ currently this will search for ALL excel sheets in the folders
%inputted, to not limit name choice, so no other excel sheets should be in
%the folders inputted or it will error 
%can update and force user to have 'spont' or 'stim' in the name if needed
%OUTPUTS-
%no variables are outputted, a .mat file, with the data from the excel now
%as a matlab file will be saved into the local exp dirs
%
%other functions needed in path to run this function:
%findFolders

%search for all excel sheets within the specified directory 
find_excel_file = findFolders(fname, '*.xlsx');

%loop to open each excel file
for i = 1:size(find_excel_file,2)
    
    %inform user of progress
    disp(['processing file ',num2str(i),'/',num2str(size( ...
        find_excel_file,2))]);
    
    %access each individual directory which contains an excel file
    x = find_excel_file{1,i};
    %save the data from the excel file into a variable called data
    data_ttt = xlsread(x);
    %transpose the data so that the haem traces are dim1, and the 
    %time/frames dim2
    %NB/ there should be 9 traces
    data_ttt = data_ttt';
    % Haem traces are 2:7
    %2=flux, 3=so2, 4=speed, 5=hbo, 6=hbr, 7=hbt
    data = data_ttt(2:7,:);
    %calculate CMRO2
    %(flux*hbr)/hbt
    data(7,:) = (data_ttt(2,:)).*(data_ttt(6,:)./data_ttt(7,:));
    
    %loco info is in channel 8
    %this trace will look different depending if readRotEncoder or virmen 
    %was used during recording (this will be sorted in later code) -- 
    %for now just extract channel
    movement = data_ttt(8,:);
    %process raw loco data
    [locomotion] = cleanLoco(movement);
    locomotion_raw = abs(diff(movement));
    locomotion_raw = [locomotion_raw,0];
    
    %stim info is in channel 9
    stim = data_ttt(9,:);
    %find when stim is on, binarise stim channel:
    stim = stim > (mean(stim)+std(stim));
    %set first 50 points to zero, in case of initial spike
    stim(1:50) = 0;
    
    %create frame vector (for plots)
    frames = [1:size(data_ttt,2)];
    
    %fps of Hb-LDF probe is 40
    fps = 40; %40Hz
    %create time vector (for plots)
    time = frames/fps;
    
    %save data into correct exp_dir
    [exp_dir,~] = fileparts(x);
    matfile = fullfile(exp_dir, 'contData');
    save(matfile,'data','movement','locomotion','locomotion_raw',...
        'stim','time','frames','fps');
    
end %end of loop excel files

end %end of function
