%This code will calculate Power Spectral Density from MEA raw data

clear
close all


pathFileNetwork = '/Users/mandarmp/Documents/mmpatil/networkdata.raw.h5';


networkData = mxw.fileManager(pathFileNetwork);
% rawFrames = double(networkData.rawMap.spikes.frameno);
% rawSpikes = double(networkData.rawMap.spikes.amplitude);

    % sampling frequency
fsNetwork = networkData.fileObj.samplingFreq;

burstEdges=GetBurstEdges(networkData);

tmin = 30.5;

tmax = 32.3;
    
% 
%     % get spike time stamps 
% tsNetwork = double(networkData.fileObj.spikes.frameno - networkData.fileObj.firstFrameNum)/fsNetwork;
%     
    
startSample = tmin*fsNetwork; % sample number from which to start data extraction
endSample = (tmax - tmin) * fsNetwork; % length of the data (s) you want to extract
    

% extract unfiltered Network data
[ntwTracesRaw, ~, ~] = networkData.extractRawData(startSample,endSample);
    

% get the mean firing rate for each electrode
meanFiringRate = mxw.activityMap.computeSpikeRate(networkData);

%getting the electrode with the highest firing rate.

[maxFiringRate, electrode] = max(meanFiringRate);

LFPData = bandpass(ntwTracesRaw(:,electrode),[1 ,100],1000);

figure(1);
subplot(2,1,1)
div_factor = 10;
plot((1:fsNetwork),LFPData(1:fsNetwork,:));
ylim([-50 50])

MUAData = bandpass(ntwTracesRaw(:,electrode),[100,700],1000);

subplot(2,1,2)
plot((1:fsNetwork),MUAData(1:fsNetwork,:));
ylim([-50 50])
% Correlate the electrode with the chip electrode

electrodeOnChip = networkData.processedMap.electrode(electrode);

%neighbouring 7 electrodes:

n = 7;

% execute
x = networkData.rawMap.map.x;
y = networkData.rawMap.map.y;
d = [];

for i = 1:length(x)
    d(i) = sqrt( (x((networkData.fileObj.map.electrode==electrodeOnChip)) - x(i))^2 + (y((networkData.fileObj.map.electrode==electrodeOnChip)) - y(i))^2 );
end
[val, ind] = sort(d);

requiredRawDataChannels= ntwTracesRaw(:,ind(1:n));


%Taking LFP of the surrounding electrodes.

LFPData2 = bandpass(requiredRawDataChannels,[1 ,100],1000);

%average it 
meanLFPData2 = mean(LFPData2,2);

%pxx = pwelch(ntwTracesRaw);

%nfft = 512; % play with this to change your frequency resolution
nfft = 512;
noverlap = round(nfft * 0.25); % 75% overlap
window = hanning(nfft);
[psd,freqs] = pwelch(meanLFPData2, 500, [], [], fsNetwork);
    %notes: F = fsNetwork/nfft
figure;
plot(log(psd))   
xlim ([0 100]);

freqs = freqs';
psd = psd';

% FOOOF settings
settings = struct();  % Use defaults
f_range = [1, 30];

% Run FOOOF
fooof_results = fooof(freqs, psd, f_range, settings);
%ylim([0 0.4]);    
% selData = Pxx(3:end,:);
% meanPSD = mean(Pxx,2);


%     % plot(F, Pxx);
%     %semilogy(F, Pxx);
%     semilogy(F,meanPSD,color);
%     % loglog(F,Pxx);
%     xlabel('Frequency (Hz)')
%     xlim ([0 11000]);
%     ylim([0 2]);
%     title(IDstring +' '+ genoStr + ' '+ DIV);
%     grid on
%     hold on
% 
% 
% end
% 
% exportFile = [opDir 'mean_PSDs.pdf']; %folder and filename for raster figures
% exportgraphics(fig, exportFile ,'Append',true,'Resolution',150)
% toc
% 
% 
% 
% 
% %% PSD code
% %signal = reshape(networkStats_opt.maxAmplitudesTimes,1,[]);
% % signal = networkStats_opt(:).';
% % n = length(signal);
% % x = signal;
% % 
% % pxx = pwelch(x);
% % pwelch(x)
% 
% %https://stackoverflow.com/questions/21750075/matlab-python-power-spectral-density-of-non-uniform-time-series