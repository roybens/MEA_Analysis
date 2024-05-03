clear
close all

dataDir ='/mnt/disk15tb/mmpatil/Syngap_Organoid/Syngap_Organoid/';
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(dataDir, '**/Network/**/*raw.h5'); 
theFiles = dir(filePattern);
error_l =[];
% set paths to the Network recording file
%pathFileNetwork =  '/mnt/disk15tb/mmpatil/Syngap_Organoid/Syngap_Organoid/230829/21609/Network/000123/data.raw.h5';
for k = 1 : length(theFiles)
% select which well to analyze (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

baseFileName = theFiles(k).name;
pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
% extract dir informationfileNames
fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
scan_runID = str2double(fileDirParts{end-1}); % extract runID
scan_runID_text = fileDirParts{end-1};
scan_chipID = str2double(fileDirParts{end-3}); % extract chipID
scan_chipID_text = fileDirParts{end-3};

fprintf(1, 'Now reading %s\n', pathFileNetwork);

% create fileManager object for the Network recording
try
    networkData = mxw.fileManager(pathFileNetwork);
catch
    error_l = [error_l string(scan_runID_text)];
    continue
end

relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);

% Let's use the time and channel vectors to visualise our data in a raster plot.

% figure('Color','w','position',[0 0 400 800]);
% subplot(2,1,1);
% plot(relativeSpikeTimes.time,relativeSpikeTimes.channel,'.','MarkerSize',2,'Color','#135ba3')
% ylabel('Channel')
% xlabel('Time [s]')
% title('Raster Plot','fontsize',11)
% xlim([0 round(max(relativeSpikeTimes.time)/4)])
% ylim([1 max(relativeSpikeTimes.channel)])
% box off;

% In a well-interconnected neuronal culture, bursts of activity will often be 
% visible by eye. In order to detect them automatically and to quantify their 
% amplitude, we take a three-step approach. First, we bin all the spike times 
% into small time windows (size adjusted by the parameter _binSize_). Note that 
% from this point on, we disregard channel number, treating all the spikes together 
% as network activity. 

binSize = 0.02;
timeVector = 0:binSize:max(relativeSpikeTimes.time);
[binnedTimes, ~] = histcounts(relativeSpikeTimes.time,timeVector);
binnedTimes(end+1) = 0;
binnedTimes = binnedTimes.';

% Second, we convolve the resulting histogram with a Gaussian kernel to produce 
% a smoothed curve of network activity. The parameter gaussianSigma, which is 
% the kernel standard deviation in seconds, determines how much the curve will 
% be smoothed; with smaller values, finer details of the burst dynamics will be 
% resolved (ie, onset and offset peaks in activity bursts), although noise may 
% be captured as well. We normalise network activity by the number of active electrodes, 
% which may vary widely from culture to culture. 

gaussianSigma = 0.16;
kernel = mxw.util.normpdf(-3*gaussianSigma:binSize:3*gaussianSigma,0,gaussianSigma); 
kernel = kernel*binSize;
firingRate = conv(binnedTimes,kernel,'same');
firingRate = firingRate/binSize;
firingRateNorm = firingRate/length(unique(networkData.rawMap.spikes.channel)); 
%save('./NetworkAct16795.mat',"firingRateNorm")
% Let's plot the Network Activity below the raster. 

% subplot(2,1,2);
% plot(timeVector,firingRateNorm,'Color','#135ba3')
% xlim([0 round(max(relativeSpikeTimes.time)/4)])
% ylabel('Firing Rate [Hz]')
% xlabel('Time [s]')
% title('Network Activity','fontsize',11)
% hold on;
% Extract the desired portion of the path


%% way to save in the pipeline.
% parts = strsplit(pathFileNetwork, '/');
% extractedPath = strjoin(parts(end-6:end-1), '/');
% newPath = ['../AnalyzedData/' extractedPath];

% Check if the directory exists, and if not, create it
% if ~exist(newPath, 'dir')
%     mkdir(newPath);
% end
path = sprintf('./TemplateMats/%s_%s_NetworkTemplates.mat', scan_runID_text , scan_chipID_text);

save(path, "firingRateNorm"); % Replace 'YourData' with your actual data

end

if ~isempty(error_l)
    error_str = strjoin(error_l,', ');
    fprintf('Unable to read file with runID: %s, file(s) skipped.\n',error_str);
end