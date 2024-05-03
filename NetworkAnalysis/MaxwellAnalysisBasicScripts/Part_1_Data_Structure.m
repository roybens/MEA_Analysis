%% MaxWell Matlab Toolboxes Tutorials
% by Elena Gronskaya [updated 25.11.2020]
% 
% PART I: Data Structure

% Welcome! This tutorial will take you through the structure of the recording 
% files. We will see how to access and extract the raw data and spike information 
% and do some basic plotting of  electrode configurations. Two types of recordings 
% will be considered in parallel: the* Activity Scan* and the *Network* recording. 
% These tutorials will address data of the new (non-legacy) format recorded with 
% either MaxOne or MaxTwo. 
% 
% Sample data can be downloaded from here: <https://share.mxwbio.com/d/dfcccbd3e5a74d2a8ae8/ 
% https://share.mxwbio.com/d/dfcccbd3e5a74d2a8ae8/> (the folders 'scan' and 'network')
% 
% For this tutorial, please make sure to use the scan_data_full file. (For the Network recording, the same file is always used).  

%% Load the data
% The recording files are stored in the HDF5 (.h5) file format. First, set the 
% path to the Activity Scan and Network files that you would like to analyze. 

clear
close all

pathFileActivityScan = '~/Documents/MATLAB/mmpatil/scan_data_full.raw.h5'; 
pathFileNetwork = '~/Documents/MATLAB/mmpatil/network_data_full.raw.h5'; 

%% Access the data
% An HDF5 file can be accessed directly in the following way.

% information about the file
networkDataInfo = h5info(pathFileNetwork);
rawpath = '/recordings/rec0000/well000/groups/routed';
% h5read(pathFileNetwork, [rawpath '/raw'],[start_samples start_channels],[nr_samples nr_channels]);
% the raw signal (e.g., first 5000 samples for first 100 channels)
networkRawData = h5read(pathFileNetwork,[rawpath '/raw'],[1 1],[5000 100]); 
% a smaller chunk of data (e.g., samples 4000-5000 for first 100 channels)
networkRawData2 = h5read(pathFileNetwork,[rawpath '/raw'],[4000 1],[1000 100]); 
% the spikes that were detected and stored by the MaxLab Live software
networkSpikes = h5read(pathFileNetwork,'/recordings/rec0000/well000/spikes');

%% The fileManager 

% There is a dedicated script that extracts the data for you (using the functions above) and structures it in an accessible way. 
% It is called fileManager.m (found in the +mxw folder of the MaxWell Matlab Toolboxes) and it is the starting point for any further data analysis with the Toolboxes. 
% To load data using the fileManager, set the file path and specify the well to be analysed:
% activityScanData = mxw.fileManager(pathFileActivityScan, wellID);
% networkData = mxw.fileManager(pathFileNetwork, wellID);
% WellID represents the selected well (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

% create fileManager object for the Activity Scan
activityScanData = mxw.fileManager(pathFileActivityScan,wellID);
% create fileManager object for the Network recording
networkData = mxw.fileManager(pathFileNetwork,wellID);

%% The fileManager data structure
% The fileManager object contains some key info (the file path, the number 
% of recordings in the file, filter information, etc.) as well as the data itself, 
% organized in four sub-structures:
% - fileObj
% - rawMap
% - processedMap
% - extractedSpikes

% Note that networkData has a single entry for all of the sub-structures 
% above, while activityScanData *has multiple entries for each of the Activity 
% Scan files. The exception is processedMap, which combines the electrode mapping 
% info across all of the recordings. We can use this mapping data to visualise 
% the location of the electrodes on the MEA.

% xy coordinates of electrodes recorded across all of the Activity Scan files
x = activityScanData.processedMap.xpos;
y = activityScanData.processedMap.ypos;

% plot xy
figure('Color','w');
subplot(2,1,1);
plot(x,y,'.','Color','#e6e7e8')
hold on

% The recording number should be set according to its order in the 'fileObj" array. 
recNr = 3;
% plot the electrodes of just one of the Activity Scan recordings
x = activityScanData.rawMap(recNr).map.x;
y = activityScanData.rawMap(recNr).map.y;
plot(x,y,'.','Color','#43b2ef');

% flip y-axis to match the electrode arrangement on the MEA
axis ij;
% fix aspect ratio between x and y to avoid image distortion
axis equal;
xlim([0 4000]);
ylim([0 2000]);
box off
xlabel('\mum');
ylabel('\mum')
title('Activity Scan Electrode Configuration')

% plot electrode configuration of the Network recording
x = networkData.processedMap.xpos;
y = networkData.processedMap.ypos;

% plot xy
subplot(2,1,2);
plot(x,y,'.','Color','#135ba3');
hold on

% find location of a given channel
chNr = 15;
ind = (networkData.rawMap.map.channel==chNr);
xCh = networkData.rawMap.map.x(ind);
yCh = networkData.rawMap.map.y(ind);
plot(xCh,yCh,'or');

% flip y-axis to match the electrode arrangement on the MEA
axis ij;
% fix aspect ratio between x and y to avoid image distortion
axis equal;
xlim([0 4000]);ylim([0 2000]);
box off
xlabel('\mum');ylabel('\mum')
title('Network Recording Electrode Configuration')

%% fileManager.processedMap
% Let's go back to the four data sub-structures within the *fileManager *object 
% and look at how they are interrelated:
% - fileObj
% - rawMap
% - processedMap
% - extractedSpikes

% For _networkData_, because there is just one recording file, _networkData.processedMap_ 
% contains the same mapping information as _networkData.rawMap.map_. The only 
% difference is that _processedMap_ orders the x and y position data according 
% to electrode number while _rawMap_ orders the data according to channel number. 

rawMap_resorted = sort(networkData.rawMap.map.electrode,'ascend');
isequal(networkData.processedMap.electrode,rawMap_resorted)

% For _activityScanData, processedMap_ combines electrode info across all of 
% the files, sorted by electrode number in ascending order.

%% fileManager.fileObj and fileManager.rawMap
% Now let's look at_ fileObj. _Most of the fields within _fileObj _are single 
% pieces of information about the recording file (like the _firstFrameNum _or_ 
% _the_ samplingFreq_). The electrode and spike data is stored in just two sub-structures: 
% _*map* _and_ *spikes*. _Note that this data is simply a copy of _rawMap._

isequal(networkData.fileObj.spikes,networkData.rawMap.spikes)
isequal(networkData.fileObj.map,networkData.rawMap.map)

% the same holds true for each of the files in the activityScanData
isequal(activityScanData.fileObj(3,1).spikes,activityScanData.rawMap(3).spikes)
isequal(activityScanData.fileObj(3,1).map,activityScanData.rawMap(3).map)

% In the spikes sub-structure, you will find three fields of matching lengths: 
% - frameno
% - channel
% - amplitude 
% These arrays pool data from all the spikes detected across all electrodes 
% in a given recording. The spike time stamps are stored as frame numbers,
% and each spike time is associated with a corresponding channel number and 
% an amplitude (in uV). 
% 
% In extractedSpikes, the spike frame numbers and amplitudes from rawMap 
% are reorganized according to ascending electrode number (same electrode order 
% as in processedMap.electrode). This is why there are only 1020 rows in
% networkData.extractedSpikes but 6600 rows in activityScanData.extractedSpikes 
% (or more, depending on the scan resolution): more electrodes were recorded 
% across all of the Activity Scan files.
% 
% To summarize:
% - fileObj contains a full copy of rawMap data
% - processedMap is the same as rawMap.map, except it combines mapping data 
% across multiple files (when present) and sorts everything according to electrode 
% number rather than channel number
% - extractedSpikes groups the spike information according to electrode number  
%% 
% Now you should be able to easily navigate the fileManager data structure, 
% identifying the electrode configuration of a recording as well as the times 
% and amplitudes of all the recorded spikes. In the next tutorial (PART II: Spikes), 
% we will go into some basic analysis and plotting of the spike data. You may 
% also be wondering, what happened to the raw signal that we learned how to extract 
% with the h5read function? Indeed, the raw signal is not a part of the fileManager 
% object. However, using the fileManager methods, you can easily extract and 
% plot the raw and bandpass-filtered data. These methods will be covered in  
% Part III: Traces.