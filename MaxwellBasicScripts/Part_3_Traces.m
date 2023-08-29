%% MaxWell Matlab Toolboxes Tutorials
% by Elena Gronskaya [updated 25.11.2020]
% 
%% PART III: Traces
% In this tutorial, we will see how to load and visualise the raw and bandpass-filtered 
% voltage traces with reference to the spike times data in _fileManager._
% 
% Sample data can be downloaded from here: <https://share.mxwbio.com/d/dfcccbd3e5a74d2a8ae8/ 
% https://share.mxwbio.com/d/dfcccbd3e5a74d2a8ae8/> (the folders 'scan' and 'network')
% 
% For this tutorial, please make sure to use the scan_data_full file. 
% (For the Network recording, the same file is always used).  
%
% As usual, we start with loading the data and creating fileManager objects. 

clear
close all

% set paths pathFileActivityScan = '~/Documents/MATLAB/mmpatil/scan_data_full.raw.h5'; 
pathFileActivityScan = '/mnt/disk15tb/mmpatil/Spikesorting/Data/activity_16757_1min_fullsccan/data.raw.h5'; 
pathFileNetwork = '/mnt/disk15tb/mmpatil/Spikesorting/Data/Mandar/Trace_20230607_13_53_45_block1.raw.h5'; 
% select which well to analyze (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

% create fileManager object for the Activity Scan
activityScanData = mxw.fileManager(pathFileActivityScan,wellID);
% create fileManager object for the Network recording
networkData = mxw.fileManager(pathFileNetwork,wellID);
%% Loading the Voltage Traces Data: all electrodes
% As you may recall, fileManager objects contain only spike data. We need 
% to use a separate fieManager method to load the raw traces and the bandpass-filtered 
% traces. Given the high sampling rate and the large number of electrodes, the 
% dataset can be very large. Thus, it is advisable to either load a small chunk 
% of data for all of the electrodes (a few seconds) or to load the entire recording 
% for just a few electrodes. Let's start with the first approach: loading 5 s 
% across all recorded electrodes.  

startSample = 1; % sample number from which to start data extraction
dataChunkSize = 5; % length of the data (s) you want to extract

%  get the sampling rate
fsNetwork = networkData.fileObj.samplingFreq;

% extract unfiltered Network data[ axonTraces, electrodeGroups, ts, w, s ] = computeAxonTraces( fileManagerObj, axonTrackElec, varargin )
    % COMPUTEAXONTRACES computes the axon traces from the neurons on top of
    % the electrodes in 'axonTrackElec'. The first step to compute the 
    % traces is 
[ntwTracesRaw, ~, ~] = networkData.extractRawData(startSample,dataChunkSize*fsNetwork);

% extract bandpass filtered Network data
[ntwTracesBP, ~, ntwElectrodesArray] = networkData.extractBPFData(startSample,dataChunkSize*fsNetwork);
 
% To load Activity Scan data, we need to use the additional
% argument 'files' to select which recording to extract. The recording 
% number should be set according to its order in the 'fileObj" array. 
recNr = 3;
% get the sampling rate
fsActivity = activityScanData.fileObj.samplingFreq;

% extract unfiltered Activity Scan data
[scanTracesRaw, ~, ~] = activityScanData.extractRawData(startSample,dataChunkSize*fsActivity,...
    'files',recNr);

% extract bandpass filtered Activity Scan data
[scanTracesBP, ~, scanElectrodesArray] = activityScanData.extractBPFData(startSample,dataChunkSize*fsActivity,...
     'files',recNr);

% The output Traces data is organised in the following way:
% - every row represents a data point 
% - every column represents an electrode in 'electrodesArray',
% which is a list of electrode numbers from which the data was extracted.
% Thus, the dimensions of the traces data are: [(data_chunk_size*fsActivity) x total_electrodes]. 
disp(dataChunkSize*fsActivity);
disp(['Total number of electrodes in Network scan: ' num2str(length(networkData.fileObj.map.electrode))]);

    % the electrodes in 'axonTrackElec'. The first step to compute the 
    % traces is in Activity scan: ' num2str(length(activityScanData.fileObj(recNr).map.electrode))]);

%% Plotting the Traces
% Let's plot and compare the raw and the bandpass-filtered signals from the 
% electrode with the highest-amplitude spikes. 

% find electrodes with biggest spike
[~, ind] = min(min(ntwTracesBP));

figure('Color','w');
subplot(2,2,1)
plot(ntwTracesRaw(:,ind),'Color','#135ba3')
xticks((0:1:dataChunkSize)*fsNetwork)
xticklabels(0:1:dataChunkSize)
xlabel('Time [s]')
ylabel('Amplitude [\muV]')
title('Network Data (Raw)')
subplot(2,2,3)
plot(ntwTracesBP(:,ind),'Color','#135ba3')
xticks((0:1:dataChunkSize)*fsNetwork)
xticklabels(0:1:dataChunkSize)
xlabel('Time [s]')
ylabel('Amplitude [\muV]')
title('Network Data (Filtered)')

% find electrodes with biggest spike
[~, ind] = min(min(scanTracesBP));
subplot(2,2,2)
plot(scanTracesRaw(:,ind),'Color','#135ba3')
xticks((0:1:dataChunkSize)*fsActivity)
xticklabels(0:1:dataChunkSize)
xlabel('Time [s]')
ylabel('Amplitude [\muV]')
title('Activity Scan Data (Raw)')
subplot(2,2,4)
plot(scanTracesBP(:,ind),'Color','#135ba3')
xticks((0:1:dataChunkSize)*fsActivity)
xticklabels(0:1:dataChunkSize)
xlabel('Time [s]')
ylabel('Amplitude [\muV]')
title('Activity Scan Data (Filtered)')
%% Loading the Voltage Traces Data: specific electrodes
% Now let's take a look at the second approach, loading a larger chunk of data 
% from just a few specific electrodes. Let's say we wanted to load 20 s, and we 
% wanted to start 10 s into the recording 

%  get the sampling rate
fs = networkData.fileObj.samplingFreq;

startSample = 10*fs; % sample number from which to start data extraction
dataChunkSize = 20; % length of the data (s) you want to extract

% find five electrodes with the strongest spikes and low signal standard
% deviation
numEl = 5;

maxSpikeAmps = abs(min(ntwTracesBP));
stdData = std(ntwTracesBP);
selectedEls = int32.empty(numEl,0);

c = 1;
while c <= numEl+1
    [biggestSpike, ind] = max(maxSpikeAmps);
    maxSpikeAmps(ind) = NaN;
    if stdData(ind)<median(stdData)
        selectedEls(c) = ind;        
        c = c+1;
    end
end

selectedEls = ntwElectrodesArray{1}(selectedEls);

% load the traces by adding the optional input variable 'electrodes', 
% as an array of electrode numbers 
[longerNtwTracesBP, ~, ~] = networkData.extractBPFData(startSample,...
              dataChunkSize*fs,'electrodes',selectedEls);
          
%% Plotting the Traces and Spike Times
% Let's plot the bandpass-filtered signal of the Network recording from the 
% five electrodes we selected. 
% 
% Just for fun, let's then check how well the spikes that were detected and 
% stored during the recording correspond with the spikes apparent in the voltage 
% traces.

figure('Color','w');
sgtitle(sprintf('20 s of Network Data from %i electrodes',numEl))
for i = 1: length(selectedEls)
    subplot(length(selectedEls),1,i)
    plot(longerNtwTracesBP(:,i),'Color','#135ba3')
    hold on
    
    % check how well the spike data in the fileManager object corresponds
    % to the spikes apparent in the voltage traces 
    
    % first, get the spike timestamps from the selected electrode:
    spikeFrames = networkData.extractedSpikes.frameno(networkData.processedMap.electrode==selectedEls(i));
    spikeFrames = spikeFrames{:};

    % next, select the spikes that fall within the chunk of extracted data    
    spikeFrames = spikeFrames(spikeFrames<(dataChunkSize*fs+networkData.fileObj.firstFrameNum+startSample)...
        &spikeFrames>networkData.fileObj.firstFrameNum+startSample);

    % adjust the spike times so that the start sample corresponds to time point 0
    spikeFrames = spikeFrames - startSample;

    % plot the spike times above the voltage traces
    scatter((spikeFrames-networkData.fileObj.firstFrameNum),(ones(length(spikeFrames),1)*max((longerNtwTracesBP(:,i)))),'r.')   
 
    xticks((0:4:dataChunkSize)*fs)
    xticklabels(0:4:dataChunkSize)
    ylabel('\muV')
    if i== length(selectedEls)       
        xlabel('Time [s]')
    end
    legend(['Electrode nr ',num2str(selectedEls(i))],'Location','eastoutside')
    box off
end

% plot the electrode position
figure('Color','w');
plot(networkData.rawMap.map.x,networkData.rawMap.map.y,'.','Color','#43b2ef');
hold on;
ind_selected_els = ismember(networkData.rawMap.map.electrode,selectedEls);
plot(networkData.rawMap.map.x(ind_selected_els),networkData.rawMap.map.y(ind_selected_els),'r.','Markersize',15);
xlabel('\mum');ylabel('\mum');axis equal;axis ij;
title('Electrode Position');
for i = 1:numEl
    text(networkData.rawMap.map.x(find(networkData.rawMap.map.electrode==selectedEls(i)))+70, ...
	networkData.rawMap.map.y(find(networkData.rawMap.map.electrode==selectedEls(i))), num2str(networkData.rawMap.map.electrode(find(networkData.rawMap.map.electrode==selectedEls(i)))) );
end
xlim([-200 4100]);
ylim([-100 2200])

%% Extracting and Plotting Spike Waveforms  
% Next, let's combine the spike information with the extracted traces to take 
% a closer look at the spike waveforms. For the five electrodes selected above, 
% we will take a 30-sample cutout (10 samples before the frame number of the spike's 
% maximum value and 20 samples after) for all spikes detected, and plot these 
% on top of each other for each electrode. If an electrode picks up activity from 
% more than one unit, we should see several distinct waveform shapes. 

figure('Color','w','Position',[0 100 300 1200]);
preSamples = 10;
postSamples = 20;
for i = 1:numEl
    subplot(numEl,1,i)
    hold on
    
    % repeat the procedure above to extract the appropriate spike frame
    % numbers
    spikeFrames = networkData.extractedSpikes.frameno(networkData.processedMap.electrode==selectedEls(i));
    spikeFrames = spikeFrames{:};    
    spikeFrames = spikeFrames(spikeFrames<(dataChunkSize*fs+networkData.fileObj.firstFrameNum+startSample)...
        &spikeFrames>networkData.fileObj.firstFrameNum+startSample);
    spikeFrames = spikeFrames - startSample;

    % for each spike, plot the waveform cut-out
    for j = 1:length(spikeFrames)
        idxStart = spikeFrames(j)-networkData.fileObj.firstFrameNum-preSamples;
        idxStop = spikeFrames(j)-networkData.fileObj.firstFrameNum+postSamples;
        if idxStart>=1 && idxStop<=length(longerNtwTracesBP)
            plot(longerNtwTracesBP((idxStart:idxStop),i),'Color','#135ba3')
        end
    end
    xlim([1 preSamples+postSamples])
    xlabel('Time [samples]')
    ylabel('Amplitude [\muV]')
    legend(['Electrode nr ',num2str(selectedEls(i))],'Location','northoutside')
end
%% Extracting Spike Waveforms using extractCutOuts
% There is a dedicated method in fileManager for extracting waveform cut-outs. 
% As input, fileManager.extractCutOuts takes three main arguments: spikesTimePoints 
% (in frame numbers), prePointsSpike and postPointsSpike (how many samples to 
% cut on the left and on the right of spikesTimePoints, respectively), as well 
% as an optional argument, 'files', in which you specify the recording number 
% you would like to work with when using Activity Scan data. Unlike the code in 
% the previous section, this method extracts the cut-outs across *all* electrodes 
% in the recording for each spike time specified. Even though many of the electrodes 
% will have no activity at a given spike's timestamp, this method is useful for 
% obtaining single unit "footprints", the distribution of a neuron's electrical 
% activity across the MEA chip. Let's take the biggest spike from each of the 
% five electrodes considered above, and see what the activity across the rest 
% of the electrodes looks like at those timestamps. 

% for each of the five selected electrodes, find the frame number of the biggest spike
selectedSpikeTimes = int64.empty(numEl,0);
for i = 1:numEl
    spikeFrames = networkData.extractedSpikes.frameno(networkData.processedMap.electrode==selectedEls(i));
    spikeFrames = spikeFrames{:};    
    spikeAmps = networkData.extractedSpikes.amplitude(networkData.processedMap.electrode==selectedEls(i));
    spikeAmps = spikeAmps{:};    
    [minVal, minIdx] = min(spikeAmps);
    selectedSpikeTimes(i) = spikeFrames(minIdx)-networkData.fileObj.firstFrameNum;
end

% at each of the five selected spike times, extract 30-sample cut-outs from all the recorded electrodes. 
preSamples = 10;
postSamples = 20;
[waveformCutOuts, electrodesArrayCutOuts] = networkData.extractCutOuts(double(selectedSpikeTimes), preSamples, postSamples);
 
% The output variable waveformCutOuts has the following structure: 
% 
% Each column contains cut-outs extracted from all of the electrodes for a given 
% time point. The cut-outs are concatenated so that the first cut-out comes from 
% the first electrode in electrodesArrayCutOuts and so on. The length of each 
% column is therefore equal to the length of each cut-out (30 samples) times the 
% number of electrodes in the file.
% 
% By default, the extractCutOuts method works with bandpass-filetered data. 
% There is also an option to extract raw waveforms, using the method extractRawCutOuts, 
% which works exactly the same way. 

%% Plotting Neural "Footprints" using plot.axonTraces
% Let's plot the five single-neural "footprints" using the cut-out spike waveforms 
% and a dedicated function, mxw.plot.axonTraces. 

% get x and y coordinates of the electrodes
x = networkData.rawMap.map.x;
y = networkData.rawMap.map.y;

% loop through all the five spike time-points
figure('Color','w','Position',[0 100 300 1200]);
for i = 1:numEl
    subplot(numEl,1,i)
    % reshape the cancatenated cut-outs array into a matrix with the
    % dimensions: cut-out length x total electrode number
    wfs = waveformCutOuts(:,i);
    wfs = reshape(wfs,preSamples+postSamples,length(wfs)/(preSamples+postSamples));
    mxw.plot.axonTraces(x,y,wfs,'PointSize', 100 ,'PlotWaveforms',true,...
        'Figure',false,'PlotHeatMap','true','Title',strcat('Neuron ',num2str(i)),...
        'WaveformColor','k','ylabel','\muV')
    axis equal;
    [val, ind] = min(min(wfs));
    xlim([x(ind) - 150 x(ind) + 150]);
    ylim([y(ind) - 150 y(ind) + 150]);
    xlabel('\mum');
    ylabel('\mum');
end
%
% Note that the locations of the footprints correspond to the locations of the 
% five electrodes we plotted above, since we used the times of each electrode's 
% biggest spike to generate our footprints. In this example we used network data, 
% so the distribution of electrodes around the biggest-spike electrode in each 
% footprint is quite sparse and uneven. Thus, for plotting footprints, it is more 
% appropriate to use Activity Scan data (or data from another recording configuration 
% with a higher spatial resolution) that has been spike-sorted (using, for example, 
% one of the spike_sorting tools available in the MaxWell Matlab Toolbox).