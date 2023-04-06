%This code will calculate Power Spectral Density from MEA raw data

clear
close all


pathFileNetwork = '/home/mmpatil/Documents/spikesorting/Data/organoid/16866/Network/000006/data.raw.h5';


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
fs = 1000; % Set the sampling rate
LFPData2 = bandpass(requiredRawDataChannels,[1 100],fs);

%average it 
meanLFPData2 = mean(LFPData2,2);

%pxx = pwelch(ntwTracesRaw);

%nfft = 512; % play with this to change your frequency resolution

% nfft = 512;
% noverlap = round(nfft * 0.25); % 75% overlap
% window = hanning(nfft);
% window_length = round(length(meanLFPData2)/4); % Set the window length to 25% of the signal length
% window_overlap = round(window_length/2); % Set the overlap to 50% of the window length
% nfft = 2^(nextpow2(window_length)); % Set the number of FFT points to the next power of 2
%
% % 
%  [p,f] = pwelch(meanLFPData2, window_length, window_overlap, nfft, fs);


[psd,freqs] = pwelch(meanLFPData2, fs,[],[],fs);
    %notes: F = fsNetwork/nfft
figure;
plot(freqs,log10(psd))   
xlim ([0 100]);
ylim([-2 2])
  

% freqs = freqs';
% psd = psd';

% FOOOF settings
% FOOOF settings
settings =struct(...
        'peak_width_limits', [2, 12], ...
        'max_n_peaks', Inf, ...
        'min_peak_height', 0.0, ...
        'peak_threshold', 2.0, ...
        'aperiodic_mode', 'fixed', ...
        'verbose', true ...
    );
f_range = [1 100];


% Run FOOOF
fooof_results = fooof(freqs', psd', f_range, settings, true);

fooof_plot(fooof_results)
% 
% if grp_proc_info_in.fooof_background_mode == 1 %fixed
%     fooofed_psd = fooof_results.peak_params(1,1) - log10(freqs.^fooof_results.peak_params(1,2)); %knee parameter is 0, and not output
% else %knee
%     fooofed_psd = fooof_results.peak_params(1,1) - log10(fooof_results.peak_params(1,2) + freqs.^fooof_results.peak_params(1,3));
%     %fooofed_psd = 10^fooof_results.background_params(1,1) * (1./(fooof_results.background_params(1,2)+freqs.^fooof_results.background_params(1,3)));
% end
% background_fit = fooofed_psd;
% for i=1:size(fooof_results.peak_params,1)
%     fooofed_psd = fooofed_psd + fooof_results.gaussian_params(i,2) * exp(-((freqs - fooof_results.gaussian_params(i,1)).^2) / (2*fooof_results.gaussian_params(i,3)).^2);
% end
% plot(freqs,log10(psd),freqs,fooofed_psd,freqs,background_fit,'--')
% xlabel('Frequency')
% ylabel('Power')
% legend('Original Spectrum', 'Full Model Fit', 'Background Fit')
% saveas(gcf, strcat(filename,'.png'))





% plot_one_foof_results(freqs,psd,f_range,fooof_results)

%https://github.com/vpapadourakis/fooof_mat/blob/master/matlab_wrapper_examples/plot_one_foof_results.m
function [] = plot_one_foof_results(freqs, psd, f_range,fooof_results)

%% get model data 
%calc L
bckgr = fooof_results.peak_params;
if numel(bckgr)<3 %offset, slope
    b = bckgr(1); chi = bckgr(2); k = 0;  
else %offset, knee, slope
    b = bckgr(1); chi = bckgr(3); k = bckgr(2);
end

L = b - log(k+(freqs.^chi));

%calc sumG
gaussians = fooof_results.gaussian_params;% this for plotting. This is used for the fit
% peaks = fooof_results.peak_params; %these are the "true" amplitudes,
% after subtracting the background (and other gaussians?)
nG = size(gaussians,1); nFreqs = numel(freqs); 
G = nan(nG,nFreqs);
for iG = 1:nG
    cf = gaussians(iG,1); a = gaussians(iG,2); w = gaussians(iG,3);
    ee = (-(freqs-cf).^2)./(2*w^2);
    G(iG,:) = a*exp(ee);
end
sumG = sum(G,1); 

%% plot
P = L + sumG; 

figure('color',[1 1 1]); hold on 
plot(freqs,log10(psd),'k'); 
plot(freqs,P,'color',[0.8 0 0.2 .5],'linewidth',2);
plot(freqs,L,'--','color',[0.2 0 0.8 .6],'linewidth',2);

legend('Original Spectrum','Full Model Fit','Background fit')
xlim(f_range); 
grid on
xlabel('Frequency'); ylabel('log10 Power');

% figure; hold on 
% plot(freqs,sumG);
% xlim(f_range);

end%function




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