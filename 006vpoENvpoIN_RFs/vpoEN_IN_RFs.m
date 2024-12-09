%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing Spike-Triggered Averages and Spectrotemporal Stimulus
% Properties driving vpoEN and vpoIN responses
%
% Related to Cell manuscript "Male-Male Interactions Shape Mate Selection
% in Drosophila" by Tom Hindmarsh Sten, Rufei Li, Florian Hollunder, Shade
% Eleazer, and Vanessa Ruta
% 
% % Description:
% This MATLAB script computes spike-triggered averages (STAs) and
% spectrotemporal receptive fields for vpoEN and vpoIN responses to playback of
% acoustic recordings from MF pair or MMF triads. The script processes
% calcium imaging and auditory stimulus recordings to derive key metrics 
% such as inter-pulse intervals (IPIs), carrier frequency (CF) distributions, 
% and their relationship to neural activity. Key steps include: 
%   - Spectrogram Analysis: Perform Short-Time Fourier Transform (STFT) on the 
%    auditory stimulus to compute power spectra and center-of-mass frequencies.
%   - IPI Computation: Calculate inter-pulse intervals (IPIs) using pulse peak 
%    detection within overlapping temporal bins.
%   - Spike-Triggered Averages: Identify calcium response peaks and compute 
%    STAs of IPIs and CFs within a defined time window around each response.
%   - Distribution Analysis: Visualize IPI and CF distributions associated with 
%    neural responses.
%   - Stimulus-Response Mapping: Correlate neural responses with IPI and CF bins 
%    to gain insight into the the multiplexed receptive fields
%   - Spectrotemporal Receptive Fields: Create spectrptermporal receptive
%     fields for vpoEN and vpoIN neurons (spectral only, no IPI considered)
%
% Visualization:
% - Temporal dynamics of IPI and CF around neural responses.
% - Distributions of response-associated IPIs and CFs.
% - Heatmaps visualizing responses based on CF/IPI combinations
% - Heatmaps visualizing spectrotemporal receptive fields.
%
% Input:
% - `vpoEN_IN_exStimResponse.mat`: Contains stimulus and response data from
%   7 animals. Contains 4 structs (vpoEN/IN responses for MF playback or MMF
%   playback). Each struct has four fields: the acoustic stimulus, calcium
%   imaging traces across 7 animals, and synchronized timestamps for both signals. 
%
% Output:
% - Plots of temporal response dynamics, histograms of IPI/CF distributions, 
%   and spectrotemporal receptive field heatmaps.
%
% Dependencies:
% - Signal Processing Toolbox for STFT and peak detection.
% - function: prctfilt.m (https://github.com/flatironinstitute/CaImAn-MATLAB/blob/master/utilities/prctfilt.m)
%
% Note:
% - The script is modular, allowing for analysis of different datasets 
%   (e.g., vpoIN_MF, vpoIN_MMF). Adjust parameters and data loading as needed.
%
%
% Author: Tom Hindmarsh Sten
% Date: 12/2/2024
% Affiliation: The Rockefeller University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start by computing the "Spike" triggered averages, where the spike is the peak responses in calcium

clear
load('vpoEN_IN_exStimResponse.mat')
data = vpoIN_MF; %can swap for vpoIN

responseCtr = 1;

% Define parameters for STFT
fs = 5000; % Sample rate
windowSize = 512; % Size of each segment (in samples)
overlap = 256; % Overlap between segments (in samples)
nfft = 1024; % Number of FFT points

% Compute spectrogram (STFT)
[spec, freq, timeFreq] = spectrogram(data.stim, windowSize, overlap, nfft, fs);

%remove noisy low frequencies
keep = find(freq>60);
spec = spec(keep,:);
freq = freq(keep);

% Convert to power
powerSpectrogram = abs(spec).^2; % or abs(spec) for amplitude

% Initialize array to store center-of-mass frequencies
comFrequency = zeros(1, length(timeFreq));

% Loop through each time point to calculate center-of-mass frequency
for t = 1:length(timeFreq)
    powerT = powerSpectrogram(:, t); % Power distribution at time t
    powerT = (powerT - min(powerT))./(max(powerT) - min(powerT)); %normalize
    powerT(powerT<1/exp(1)) = 0; %remove low-power freqs
    comFrequency(t) = sum(freq .* powerT) / sum(powerT); % Weighted average
end

% Now compute a running average of the IPI. Start by defining
% IPI bin parameters: 200 ms bins with 50% overlap
minPulseDist = 0.005;
binSize = 0.2;
overlap = 0.1;
binStarts = data.stimT(1):overlap:(data.stimT(end) - binSize); % Bin start times
binCentersIPI = binStarts+overlap;
nBins = length(binStarts);

% Initialize an array to store the average IPI for each bin
mnIPI = zeros(1, nBins);

% Find pulses in the sound trace
[~, pulseTimes] = findpeaks(data.stim, data.stimT, 'MinPeakDistance',...
    minPulseDist ,'MinPeakProminence',0.03);

% Loop through each bin and compute the average IPI
for b = 1:nBins
    % Define bin boundaries
    binStart = binStarts(b);
    binEnd = binStart + binSize;

    % Find pulses within this bin
    pulsesInBin = pulseTimes(pulseTimes >= binStart & pulseTimes <= binEnd);

    % If pulses are detected, compute IPIs
    if length(pulsesInBin) > 1
        IPI_inBin = diff(pulsesInBin); % Compute IPIs within the bin
        mnIPI(b) = mean(IPI_inBin); % Average IPI for the bin
    else
        % If no pulses or only one pulse, set IPI to 200 ms
        mnIPI(b) = 0.2;
    end
end


% Now, lets do an STA for IPI and CF
for j = 1:size(data.response,1)

    %find responses
    thresh = 0.3; %in DFF0
    responseData = prctfilt(data.response(j,:),25,100); %steady the baseline
    [responsePeaks, responseTimes] = findpeaks(responseData,data.imT,'MinPeakHeight',...
        thresh,'MinPeakProminence',thresh);

    %remove responses too close to each other
    minDist = 1; %sec
    selectedPeaks = [];
    selectedTimes = [];
    c = 1;
    while c <= length(responseTimes)
        % Select the current peak and its location
        selectedPeaks = [selectedPeaks, responsePeaks(c)];
        selectedTimes = [selectedTimes, responseTimes(c)];

        % Skip peaks that are within the minimum distance from the current peak
        c = c + find(responseTimes(c+1:end) - responseTimes(c) > minDist, 1);
    end

    %upsample responses
    imTfast = min(data.imT):0.02:max(data.imT);
    imRespFast = interp1(data.imT,data.response(j,:),imTfast);

    %upsample IPIs
    IPIBinsFast = min(binCentersIPI):0.02:max(binCentersIPI);
    IPIfast = interp1(binCentersIPI,mnIPI,IPIBinsFast);

    %Compute average IPI and Freq in a window surrounding each response
    for k = 1:length(selectedTimes)

        window = selectedTimes(k)-1:0.05:selectedTimes(k)+0.5;

        for w = 1:length(window)-1

            idxIm = find(imTfast>window(w) & imTfast<window(w+1));
            idxFreq = find(timeFreq>window(w) & timeFreq<window(w+1));
            idxIPI = find(IPIBinsFast>window(w) & IPIBinsFast<window(w+1));

            IPI_STA(w,responseCtr) = mean(IPIfast(idxIPI));
            Freq_STA(w,responseCtr) = mean(comFrequency(idxFreq));
            Response_STA(w,responseCtr) = mean(imRespFast(idxIm));

        end

        %Since the lag is ~450-550ms, lets store the average frequency
        %and IPI in the 400-600ms preceding each response to build a
        %distribution.
        IPIDist(responseCtr,1) = mean(IPI_STA(9:13,responseCtr));
        FreqDist(responseCtr,1) = nanmean(Freq_STA(9:13,responseCtr));

        responseCtr = responseCtr+1;
    end
end


tBasis = -1:0.05:0.45;

%Plot STAs
figure;
subplot(1,3,1);
plot(tBasis,nanmean(IPI_STA,2),'r','LineWidth',2)
xlabel('Time from Response');
ylabel('IPI')
box off
title('Inter-pulse interval'); hold on
axis([-1 0.5 0 .2])

subplot(1,3,2);
plot(tBasis,nanmean(Freq_STA,2),'k','LineWidth',2)
xlabel('Time from Response');
ylabel('CF (Hz)')
box off
title('Frequency'); hold on
axis([-1 0.5 0 400])

subplot(1,3,3);
plot(tBasis,nanmean(Response_STA,2),'g','LineWidth',2)
xlabel('Time from Response');
ylabel('DFF')
box off
title('Response'); hold on
set(gcf,'Color','w');
axis([-1 0.5 0 1])

%Plot distributions
figure;
subplot(1,2,1);
histogram(IPIDist,0:0.02:0.2,'Normalization','probability')
xlabel('IPI (sec)');
ylabel('Probability')
axis([0 0.2 0 0.8])
box off
title('Inter-pulse interval'); hold on

subplot(1,2,2);
histogram(FreqDist,0:20:800,'Normalization','probability')
xlabel('Carrier Freq (Hz)');
ylabel('Probability')
title('Carrier Frequency'); hold on
axis([0 500 0 0.8])

set(gcf,'Color','w'); box off

%% Get the evoked activity for each IPI/Frequency pair

clear;
load('vpoEN_IN_exStimResponse.mat'); % Load aligned data
recordingCtr = 1; % Counter for recordings

data = vpoEN_MMF; % Change to process a different batch of experiments (EN/IN –/MF/MMF)

% Parameters for Short-Time Fourier Transform (STFT)
fs = 5000; % Sampling rate (Hz)
windowSize = 512; % Size of each segment (samples)
overlap = 256; % Overlap between segments (samples)
nfft = 1024; % Number of FFT points

% Compute spectrogram
[spec, freq, timeFreq] = spectrogram(data.stim, windowSize, overlap, nfft, fs);

% Remove noisy low frequencies (below 60 Hz)
keep = freq > 60;
spec = spec(keep, :);
freq = freq(keep);

% Convert to power spectrum
powerSpectrogram = abs(spec).^2; % Use squared magnitude for power

% Initialize center-of-mass frequency array
comFrequency = zeros(1, length(timeFreq));

% Calculate center-of-mass frequency for each time point
for t = 1:length(timeFreq)
    powerT = powerSpectrogram(:, t); % Power distribution at time t
    powerT = (powerT - min(powerT)) / (max(powerT) - min(powerT)); % Normalize
    powerT(powerT < 1/exp(1)) = 0; % Suppress low-power frequencies
    comFrequency(t) = sum(freq .* powerT) / sum(powerT); % Weighted average
end

% Parameters for IPI computation
minPulseDist = 0.005; % Minimum distance between pulses (s)
binSize = 0.2; % Bin size (s)
overlap = 0.1; % Overlap between bins (s)
binStarts = data.stimT(1):overlap:(data.stimT(end) - binSize); % Bin start times
binCentersIPI = binStarts + overlap; % Bin centers
nBins = length(binStarts);

% Initialize array for mean IPI
mnIPI = zeros(1, nBins);

% Detect pulses in the sound trace
[pulsePeaks, pulseTimes] = findpeaks(data.stim, data.stimT, ...
    'MinPeakDistance', minPulseDist, 'MinPeakProminence', 0.03);

% Compute average IPI for each bin
for b = 1:nBins
    binStart = binStarts(b);
    binEnd = binStart + binSize;

    % Find pulses in the current bin
    pulsesInBin = pulseTimes(pulseTimes >= binStart & pulseTimes <= binEnd);

    % Compute IPIs if there are multiple pulses
    if length(pulsesInBin) > 1
        IPI_inBin = diff(pulsesInBin); % Inter-pulse intervals
        mnIPI(b) = mean(IPI_inBin); % Average IPI
    else
        mnIPI(b) = 0.2; % Default IPI for empty bins
    end
end

% Compute response for each IPI/Carrier Frequency bin
for j = 1:size(data.response, 1)

    % Upsample responses to a common time base
    dT = 0.1; % Time step (s)
    commonT = 0:dT:max(data.imT); % Common time vector
    imRespFast = interp1(data.imT, data.response(j, :), commonT); % Upsampled response
    IPIFast = interp1(binCentersIPI, mnIPI, commonT); % Upsampled IPI
    FreqFast = interp1(timeFreq, comFrequency, commonT); % Upsampled frequency

    % Define bins for frequency and IPI
    FreqBins = 60:10:800; % Frequency bins (Hz)
    IPIBins = 0:0.005:0.15; % IPI bins (s)

    % Compute average response for each bin with a 550 ms delay
    for f = 1:length(FreqBins) - 1
        for p = 1:length(IPIBins) - 1
            offset = round(0.5 / dT); % Delay offset
            idx = find(IPIFast > IPIBins(p) & IPIFast <= IPIBins(p + 1) & ...
                FreqFast > FreqBins(f) & FreqFast <= FreqBins(f + 1));
            if ~isempty(idx)
                mnResp(f, p, recordingCtr) = mean(imRespFast(offset + idx)); % Average response
            end
        end
    end

    recordingCtr = recordingCtr + 1; % Increment recording counter
end

% Plot the results
tmp = nanmean(mnResp, 3); % Average across recordings
tmp = fillmissing(tmp, 'constant', 0); % Fill missing values with zeros

figure;
imagesc(FreqBins(2:end), IPIBins(2:end), imgaussfilt(tmp, 1)'); % Smoothed plot
axis xy; % Flip y-axis for standard orientation

% Define colormap for the plot
INMap = [linspace(1, 55/255, 256); linspace(1, 125/255, 256); linspace(1, 113/255, 256)];
ENMap = [linspace(1, 90/255, 256); linspace(1, 52/255, 256); linspace(1, 112/255, 256)];
colormap(ENMap'); box off;
clim([0 .1])

xlabel('Carrier Frequency (Hz)'); % Label x-axis
ylabel('IPI (s)'); % Label y-axis
set(gcf,'Color','w')

%% Compute canonical spectrotemporal receptive fields

clear;
load('vpoEN_IN_exStimResponse.mat'); % Load aligned data
responseCtr = 1; % Counter for recordings

% Load dataset for IN experiments
data = vpoIN_MF; % Adjust to process a different batch (e.g., EN/IN - MF/MMF)

% Parameters for Short-Time Fourier Transform (STFT)
fs = 5000; % Sampling rate (Hz)
windowSize = 256; % Size of each segment (samples)
overlap = 128; % Overlap between segments (samples)
nfft = 512; % Number of FFT points

% Compute spectrogram
[spec, freq, timeFreq] = spectrogram(data.stim, windowSize, overlap, nfft, fs);

% Remove low frequencies (if needed)
keep = freq > 0;
spec = spec(keep, :);
freq = freq(keep);

% Convert spectrogram to power spectrum
powerSpectrogram = abs(spec).^2; % Use squared magnitude for power

% Loop through recordings in the response data
for j = 1:size(data.response, 1)

    % Filter response data using a 25th percentile filter
    thresh = 0.3; % Threshold for peak detection
    responseData = prctfilt(data.response(j, :), 25, 100);
    
    % Find peaks in the response data
    [responsePeaks, responseTimes] = findpeaks(responseData, data.imT, ...
        'MinPeakHeight', thresh, 'MinPeakProminence', thresh);
    warning('off', 'signal:findpeaks:largeMinPeakHeight'); % Suppress warnings

    % Remove responses that are too close to each other
    minDist = 1; % Minimum distance between responses (seconds)
    selectedPeaks = [];
    selectedTimes = [];
    c = 1; % Index for response times
    while c <= length(responseTimes)
        % Add the current peak and its location to the selected list
        selectedPeaks = [selectedPeaks, responsePeaks(c)];
        selectedTimes = [selectedTimes, responseTimes(c)];

        % Skip peaks closer than minDist
        c = c + find(responseTimes(c+1:end) - responseTimes(c) > minDist, 1);
    end

    % Remove times outside the stimulus playback window
    selectedTimes(selectedTimes + 0.5 > max(data.stimT)) = [];

    % Upsample response data for finer time resolution
    imTfast = min(data.imT):0.02:max(data.imT);
    imRespFast = interp1(data.imT, data.response(j, :), imTfast);

    % Compute average frequency spectrum before each response
    for k = 1:length(selectedTimes)
        % Define a time window around the response
        window = selectedTimes(k)-1:0.05:selectedTimes(k)+0.5;

        for w = 1:length(window)-1
            % Get indices for imaging and frequency spectrogram
            idxIm = find(imTfast > window(w) & imTfast < window(w+1));
            idxFreq = find(timeFreq > window(w) & timeFreq < window(w+1));

            % Compute the mean spectrum within the current window
            mnSpectrum = nanmean(powerSpectrogram(:, idxFreq), 2);

            % Store results
            Freq_PSpectrum(w, :, responseCtr) = mnSpectrum;
            Response_Spectrum(w, responseCtr) = mean(imRespFast(idxIm));
        end

        % Increment response counter
        responseCtr = responseCtr + 1;
    end
end

% Remove extreme outliers (values > 5 * standard deviation)
for i = 1:size(Freq_PSpectrum, 1)
    for j = 1:size(Freq_PSpectrum, 2)
        thresh = mean(squeeze(Freq_PSpectrum(i, j, :))) + std(squeeze(Freq_PSpectrum(i, j, :))) * 5;
        throw = find(Freq_PSpectrum(i, j, :) > thresh);
        Freq_PSpectrum(i, j, throw) = 0; % Set outlier values to 0
    end
end

% Plot average frequency spectrum over time
figure;
tBasis = -1:0.05:0.45; % Time axis (relative to response)
idxPlot = find(freq < 800); % Plot frequencies below 800 Hz
imagesc(tBasis, freq(idxPlot), nanmean(Freq_PSpectrum(:, idxPlot, :), 3)');
axis xy; % Ensure y-axis is oriented correctly
colormap(gray); % Use grayscale colormap
xlabel('Time from response (sec)');
ylabel('Frequency (Hz)');
clim([0 0.2]); % Adjust color limits
