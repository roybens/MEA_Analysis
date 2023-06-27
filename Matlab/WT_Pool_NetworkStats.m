clear
close all

%Read the project names from the files and the DIV start date of the
%project
refFile = '/home/mmp/Documents/wt_network_sup.xlsx';


% extract runID info from reference excel sheet
T = readtable(refFile);
ProjectFilePaths = T.(1);
DIV0 = T.(2);
CHIPS = T.(3);
DIVOPT = T.(4);
DIVOPT = datestr(DIVOPT,'yymmdd');
%select the DIV to be scanned. eg : 25 or 27
fileNames= {};
fileFolders ={};
for k = 1 : length(ProjectFilePaths)
    % Iterate over the CHIPID values
    chip_names  = char(CHIPS(k));
    chips = strsplit(chip_names,',');
    for chip = chips
        % Generate the file pattern with the updated chip ID
        filePattern = fullfile(ProjectFilePaths(k), sprintf('**/%s/%s/Network/**/*raw.h5', DIVOPT(k,:), char(chip)));
        filePattern = char(filePattern);
        
        % Get the files that match the pattern
        files = dir(filePattern);
        
        % Append the file names to the fileNames cell array
        fileNames = [fileNames;{files.name}'];
        fileFolders = [fileFolders;{files.folder}'];
    end

    for j = 1 : length(fileNames)
        % reset recording info
        scan_runID = nan;
        scan_chipID = nan;
        meanIBI = nan;
        meanBurstPeak = nan;
        nBursts = nan;
        meanSpikesPerBurst = nan;
        spikesPerBurst = nan;
        hd5Date = nan; 
        scan_div = nan;
    
        baseFileName = fileNames(j).name;
        pathFileNetwork = fullfile(fileFolders(j), baseFileName);
        % extract dir information
        fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
        scan_runID = str2double(fileDirParts{end-1}); % extract runID
        scan_runID_text = fileDirParts{end-1};
        scan_chipID = str2double(fileDirParts{end-3}); % extract chipID
    end



end

%CHIP IDS of interest that too from the file.

