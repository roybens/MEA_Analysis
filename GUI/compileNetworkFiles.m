function [] = compileNetworkFiles(data)
%COMPILENETWORKFILES Summary of this function goes here
%   Detailed explanation goes here
    

 

    %ui components.
    fig = data.fig;
    d = uiprogressdlg(fig,'Title','Compiling files',...
        'Message','Start','Cancelable','on');
    drawnow

    % Unpack the data structure
    projectName = data.projectName;
    div0Date = data.div0Date;
    parentFolderPath = data.parentFolderPath;
    refDir = data.refDir;
    opDir = data.opDir;
    gaussianSigma = data.gaussianSigma;
    binSize = data.binSize;
    minPeakDistance = data.minPeakDistance;
    thresholdBurst = data.thresholdBurst;
    use_fixed_threshold = false;
    thresholdStartStop = data.thresholdStartStop;
    % Set Threshold function for later use
    threshold_fn = 'Threshold';
    if use_fixed_threshold 
        threshold_fn = 'FixedThreshold';
    end
    
    % make output folder
    if not(isfolder(append(opDir,'Network_outputs/Raster_BurstActivity/Plot60s')))
        mkdir(append(opDir,'Network_outputs/Raster_BurstActivity/Plot60s/'));
    end
       % make output folder
    if not(isfolder(append(opDir,'Network_outputs/Raster_BurstActivity/Plot120s')))
        mkdir(append(opDir,'Network_outputs/Raster_BurstActivity/Plot120s'));
    end
       % make output folder
    if not(isfolder(append(opDir,'Network_outputs/Raster_BurstActivity/Plot300s')))
        mkdir(append(opDir,'Network_outputs/Raster_BurstActivity/Plot300s'));
    end
           % make output folder
    if not(isfolder(append(opDir,'Network_outputs/Raster_BurstActivity/Plot600s')))
        mkdir(append(opDir,'Network_outputs/Raster_BurstActivity/Plot600s'));
    end
%     %open the log file 
%     logFile = fopen(append(opDir,'/log_file.txt'),'A');
% 
%     % Check if the file was opened successfully
%     if logFile == -1
%         error('Could not open log file.');
%     end
    logFile = data.logFile;
    % extract runID info from reference excel sheet
    T = readtable(refDir);
    run_ids = unique(T.(4));
    run_id_and_type = [T(:,[4,7])];
    
    % defines
    % convert div 0 to date datatype
    div0_date = datetime(div0Date, "InputFormat",'MM/dd/yyyy');
    % create table elements
    Run_ID = [];
    DIV = [];
    Well =[];
    NeuronType={};
    Time = [];
    Chip_ID = {};
    IBI = [];
    Burst_Peak = [];
    Number_Bursts = [];
    Spike_per_Burst = [];
    
    % create a list to catch error runIDs
    error_l = [];
    
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
        % Check for Cancel button press
        if d.CancelRequested
            break
        end
        % Update progress, report current estimate
        d.Value = k/length(theFiles);
        
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
    

        % extract dir informationfileNames
        baseFileName = theFiles(k).name;%pathFileNetwork
        pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
        
        fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
        scan_runID = str2double(fileDirParts{end-1}); % extract runID
        scan_runID_text = fileDirParts{end-1};
        scan_chipID = fileDirParts{end-3}; % extract chipID
    
        if ismember(scan_runID,run_ids)
            fprintf(1, 'Now reading %s\n', pathFileNetwork);
            fprintf(logFile, 'Now reading %s\n', pathFileNetwork);
            % create fileManager object for the Network recording
            
            idx = T.Run_ == scan_runID;

            wellsIDs = T.Wells_Recorded(idx);
            if iscell(wellsIDs)
            if ismember(',', wellsIDs{1})
            wellsIDs = strsplit(wellsIDs{1}, ',');
            wellsIDs = cellfun(@str2double,wellsIDs);
            else
                error('wellsIDs are not comma separated correctly');
            end
            end
            
            neuronTypes = T.NeuronSource(idx);
            if ismember(',', neuronTypes{1})
            neuronTypes = strsplit(neuronTypes{1}, ',');
            end
            for z = 1:length(wellsIDs)
            wellID=wellsIDs(z);
            fprintf(logFile, 'Processing Well %d\n', wellID);
            neuronSourceType = neuronTypes(z);
            try
                networkData = mxw.fileManager(pathFileNetwork,wellID);
            catch
                error_l = [error_l string(scan_runID_text)];
                continue
            end
           
            % get the startTime of the recordings
            hd5_time = networkData.fileObj.stopTime;
            try
                hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
            catch
                hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            end
            scan_div = fix(daysact(div0_date , hd5Date));
        
    
            % compute Network Activity and detect bursts
            relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);
            networkAct = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
            networkStats = computeNetworkStats_JL(networkAct, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
        
            
            %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
            %average IBI
            meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
            %average Burst peak (burst firing rate y-value)
            meanBurstPeak = mean(networkStats.maxAmplitudesValues);
            %Number of bursts
            nBursts = length(networkStats.maxAmplitudesTimes);
            
    
    
            %average spikesPerBurst        
            if length(networkStats.maxAmplitudesTimes)>3
                peakAmps = networkStats.maxAmplitudesValues';
                peakTimes = networkStats.maxAmplitudesTimes;
                
                % get the times of the burst start and stop edges
                edges = double.empty(length(peakAmps),0);
                for i = 1:length(peakAmps)
                   % take a sizeable (±6 s) chunk of the network activity curve 
                   % around each burst peak point
                   idx = networkAct.time>(peakTimes(i)-6) & networkAct.time<(peakTimes(i)+6);
                   t1 = networkAct.time(idx);
                   a1 = networkAct.firingRate(idx)';
                  
                   % get the arunIDstemp = run_id_and_type(:,1);mplitude at the desired peak width
                   peakWidthAmp = (peakAmps(i)-round(peakAmps(i)*thresholdStartStop));
                   
                   % get the indices of the peak edges
                   idx1 = find(a1<peakWidthAmp & t1<peakTimes(i));
                   idx2 = find(a1<peakWidthAmp & t1>peakTimes(i));
                   
                   if ~isempty(idx1)&&~isempty(idx2)       
                       tBefore = t1(idx1(end));
                       tAfter = t1(idx2(1));
                       edges(i,[1 2]) = [tBefore tAfter];
                   end
                end
                
               % identify spikes that fall within the bursts
                ts = ((double(networkData.fileObj.spikes.frameno)...
                    - double(networkData.fileObj.firstFrameNum))/networkData.fileObj.samplingFreq)';
                ch = networkData.fileObj.spikes.channel;
                
                spikesPerBurst = double.empty(length(edges),0);
                tsWithinBurst = [];
                chWithinBurst = [];
                for i = 1:length(edges)
                   idx = (ts>edges(i,1) & ts<edges(i,2));
                   spikesPerBurst(i) = sum(idx); 
                   tsWithinBurst = [tsWithinBurst ts(idx)];
                   chWithinBurst = [chWithinBurst ch(idx)'];
                end
                meanSpikesPerBurst = mean(spikesPerBurst);
            end
       
            % append information to table elements
            Run_ID = [Run_ID scan_runID];
            DIV = [DIV scan_div];
            Well = [Well wellID];
            NeuronType{end+1} = neuronSourceType{1};
            Time = [Time hd5Date];
            Chip_ID{end+1}=scan_chipID;
            IBI = [IBI meanIBI];
            Burst_Peak = [Burst_Peak meanBurstPeak];
            Number_Bursts = [Number_Bursts nBursts];
            Spike_per_Burst = [Spike_per_Burst meanSpikesPerBurst];
%             runIDstemp = run_id_and_type(:,1).Variables;
%             types = run_id_and_type(:,2).Variables; % to do fix columnnames
%             index = find(runIDstemp == scan_runID);
%             targetType = types{index};
            % plot results
                figure('Color','w','Position',[0 0 400 800],'Visible','off');
                subplot(2,1,1);
                mxw.plot.rasterPlot(networkData,'Figure',false);
                box off;
                %xlim([0 round(max(relativeSpikeTimes.time)/4)])
                xlim([0 60])
                ylim([1 max(relativeSpikeTimes.channel)])
                
                subplot(2,1,2);
                mxw.plot.networkActivity(networkAct,'Threshold',thresholdBurst,'Figure',false);
                box off;
                hold on;
                plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
                %xlim([0 round(max(relativeSpikeTimes.time)/4)])
                xlim([0 60])
                ylim([0 20])
                saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/Plot60s/Raster_BurstActivity_',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.png'))
                subplot(2,1,1);
                xlim([0 120])
                ylim([1 max(relativeSpikeTimes.channel)])
                subplot(2,1,2);
                xlim([0 120])
                ylim([0 20])
                saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/Plot120s/Raster_BurstActivity_',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.png'))
                
                %savefig(appen
                subplot(2,1,1);
                xlim([0 300])
                ylim([0 max(relativeSpikeTimes.channel)])
                subplot(2,1,2);
                xlim([0 300])
                ylim([0 20])
                saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/Plot300s/Raster_BurstActivity_',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.png'))
                 subplot(2,1,1);
                  xlim([0 600])
                  ylim([0 max(relativeSpikeTimes.channel)])
                  subplot(2,1,2);
                  xlim([0 600])
                  ylim([0 20])
                  saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/Plot600s/Raster_BurstActivity_',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.png'))
%                 
%                 
                
                %savefig(append(opDir,'Network_outputs/Raster_BurstActivity/Raster_BurstActivity',scan_runID_text,'.fig'))
            end
        end
        if k==2
            a=0;
        end
    
    if d.CancelRequested
        d.close()
        fprintf(logFile,'User Interruption')
        error("User interruption")
    end
    end
    %% construct table
    % convert row list to columns
    Run_ID = Run_ID';
    DIV = DIV';
    Well = Well';
    NeuronType= NeuronType';
    Time = Time';
    Chip_ID = Chip_ID';
    IBI = IBI';
    Burst_Peak = Burst_Peak';
    Number_Bursts = Number_Bursts';
    Spike_per_Burst = Spike_per_Burst';
    % make table
    T = table(Run_ID,DIV,Well,NeuronType,Time,Chip_ID,IBI,Burst_Peak,Number_Bursts,Spike_per_Burst);
    T = sortrows(T,"Run_ID","ascend");
    writetable(T, fullfile(opDir,'Network_outputs/Compiled_Networks.csv'));
    
    if ~isempty(error_l)
        error_str = strjoin(error_l,', ');
        fprintf(logFile,'Unable to read file with runID: %s, file(s) skipped.\n',error_str);
        fprintf('Unable to read file with runID: %s, file(s) skipped.\n',error_str);
    end
    fprintf(logFile,'Network analysis successfully compiled.\n');
    fprintf('Network analysis successfully compiled.\n');




end

