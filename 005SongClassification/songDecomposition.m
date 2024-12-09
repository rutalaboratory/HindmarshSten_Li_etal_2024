%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Detecting Ideal Pulse and Sine Segments in Drosophila Song
%
% Related to Cell manuscript "Male-Male Interactions Shape Mate Selection
% in Drosophila" by Tom Hindmarsh Sten, Rufei Li, Florian Hollunder, Shade
% Eleazer, and Vanessa Ruta
% 
% % Description:
% This script processes sound signals from an experimental recording. It computes 
% spectrograms, applies song template filters, scores matches using dynamic time 
% warping (DTW), and classifies song sequences into distinct categories (e.g., 
% pulse or sine songs). Additionally, it performs envelope analysis to estimate 
% the fraction of sound bouts containing pure pulse song. Key steps
% include: 
% - Computing spectrogram of input data
% - Computing approximate spectral densities of idealized pulse, sine, and
%   agonist song. 
% - Deriving filters for pulse and sine song, applying them to real and
%  idealized signals. 
% - Segmenting input data into bouts of sine and pulse song. 
% - Computing fraction of sound bouts containing pure pulse song. 
%
% Visualization:
% - Approximated spectral templates for pulse, sine, and agonist songs.
% - Filtered spectrograms with song templates.
% - Classified pulse and sine song sequences over time.
% - Envelope-based analysis of sound signals.
%
% Input:
% - `exRecording_MF.mat`: Contains sound recordings, sampling frequency, and copulation time.
%   Can be replaced by exRecording_MMF.mat.   
% - `canonicalSongTraces.mat`: Contains canonical song traces for template
%   generation. Contains one struct with examples of each song type
%
% Output:
% - Plots of filtered spectrograms, song templates, and classified sequences.
% - Metrics such as DTW scores and the fraction of sound bouts with pure pulse songs.
%
% Dependencies:
% - Signal Processing Toolbox for STFT and peak detection.
% - Customized function: fitGaussianSongSpectra.m
% - Song statistics: canonicalSongTraces.mat
%
%
% Author: Tom Hindmarsh Sten
% Date: 12/2/2024
% Affiliation: The Rockefeller University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing and loading sound signals

% Load data and clear variables
load('exRecording_MMF.mat');
clear sineScore pulseScore;

% Settings
Fs = exRecording.Fs; % Sampling frequency (Hz)
NDft = 1024; % Number of FFT points for spectrogram
sampleWindow = 500; % Spectrogram window size
sampleOverlap = 400; % Overlap between spectrogram windows
filterFreqRange = [80, 1000]; % Bandpass filter frequency range (Hz)

% Select the dataset to process
i = 30; % Index of the dataset
t = (1/Fs):(1/Fs):(length(exRecording.sound)/Fs); % Time vector (s)

% Find time of copulation
if ~isnan(exRecording.copulationTime)
    idxCop = find(t > exRecording.copulationTime, 1);
else
    idxCop = length(exRecording.sound);
end

% Extract sound sequence up to copulation
sSequence = exRecording.sound(1:idxCop);
sSequence = bandpass(sSequence, filterFreqRange, Fs);
t = t(1:idxCop);

% Compute spectrogram of the tester sequence
[s, f] = spectrogram(sSequence, sampleWindow, sampleOverlap, NDft, Fs);

% Define timepoints for spectrogram columns
allStarts = 0:(sampleWindow-sampleOverlap):(length(sSequence)-sampleWindow);
allEnds = sampleWindow:(sampleWindow-sampleOverlap):length(sSequence);

%% Compute spectral features for pulse, sine, and agonist templates

% Load canonical song traces for pulse, sine, and agonist songs
load('canonicalSongTraces.mat');

[pulseSpectrum, fPulse] = fitGaussianSongSpectra(gtPulseBouts, NDft, Fs, sampleWindow, sampleOverlap);
[sineSpectrum, fSine] = fitGaussianSongSpectra(gtSineBouts, NDft, Fs, sampleWindow, sampleOverlap);
[agonistSpectrum, fAgonist] = fitGaussianSongSpectra(gtAgonistBouts, NDft, Fs, sampleWindow, sampleOverlap);

% Normalize spectra
normalize = @(x) rescale(x, 0, 100);
pulseSpectrum = normalize(pulseSpectrum);
sineSpectrum = normalize(sineSpectrum);
agonistSpectrum = normalize(agonistSpectrum);

% Plot approximated pulse/sine/agonist spectra
figure; 
plot(fPulse(:,1),pulseSpectrum,'b','LineWidth',2); hold on;
plot(fSine(:,1),sineSpectrum,'r','LineWidth',2); 
plot(fAgonist(:,1),agonistSpectrum,'g','LineWidth',2); 
legend('Pulse','Sine','Agonist');
xlabel('Freq (Hz)'); ylabel('Density (au)')

%% Combine filters

% Apply weighting to compute cost functions
sineFilterWeighted = sineSpectrum - (pulseSpectrum + agonistSpectrum);
pulseFilterWeighted = pulseSpectrum - (sineSpectrum + agonistSpectrum);

% Re-normalize to the range [-100, 100]
renormalize = @(x) rescale(x, -100, 100);
sineFilterWeighted = renormalize(sineFilterWeighted);
pulseFilterWeighted = renormalize(pulseFilterWeighted);

% Apply frequency range restrictions
validFreqs = (fPulse(:,1) >= 80 & fPulse(:,1) <= 450);
sineFilterWeighted(~validFreqs) = 0;
pulseFilterWeighted(~validFreqs) = 0;

% Plot filters
figure; 
plot(fPulse(:,1),pulseFilterWeighted,'b','LineWidth',2); hold on;
plot(fSine(:,1),sineFilterWeighted,'r','LineWidth',2); 
legend('Pulse Filter','Sine Filter');
xlabel('Freq (Hz)'); ylabel('Density (au)')

%% Filter spectrogram with sine and pulse templates

% Normalize spectrogram between -1 and 1
normSpectrogram = 2.*(abs(s)-min(abs(s)))./(max(abs(s))-min(abs(s))) - 1;

% Apply filters to spectrogram
sineFilteredSpectra = normSpectrogram .* sineFilterWeighted;
pulseFilteredSpectra = normSpectrogram .* pulseFilterWeighted;

%% Score matches using dynamic time warping (DTW)

% Prepare idealized template for pulse
burstTrain = gtPulseBouts(8).trace(1:sampleWindow);
[sPrimePulse,~] = spectrogram(burstTrain,sampleWindow,sampleOverlap,NDft,5000);
sPrimePulse = 2*((abs(sPrimePulse) - min(abs(sPrimePulse)))./(max(abs(sPrimePulse)) - min(abs(sPrimePulse)))) -1;

% Prepare idealized template for some
sineTrain = gtSineBouts(6).trace(1:sampleWindow);
[sPrimeSine,~] = spectrogram(sineTrain,sampleWindow,sampleOverlap,NDft,5000);
sPrimeSine = 2*((abs(sPrimeSine) - min(abs(sPrimeSine)))./(max(abs(sPrimeSine)) - min(abs(sPrimeSine)))) -1;

% Run idealized copy through same filter as data
sPrimeSineFiltered = sPrimeSine.*sineFilterWeighted;
sPrimePulseFiltered = sPrimePulse.*pulseFilterWeighted;

% Compute DTW scores for each time point
pulseScore = zeros(1, size(pulseFilteredSpectra, 2));
sineScore = zeros(1, size(sineFilteredSpectra, 2));
for k = 1:size(pulseFilteredSpectra, 2)
    pulseScore(k) = dtw(sPrimePulseFiltered(validFreqs), pulseFilteredSpectra(validFreqs, k));
    sineScore(k) = dtw(sPrimeSineFiltered(validFreqs), sineFilteredSpectra(validFreqs, k));
end

% Smooth and z-score the scores
smoothZScore = @(x) smooth(zscore(x * -1));
pulseScore = smoothZScore(pulseScore);
sineScore = smoothZScore(sineScore);

% Plot the results
figure;
plot(sSequence, 'k');
yyaxis right;
plot(allStarts, sineScore, 'r-', 'DisplayName', 'Sine Score');
hold on;
plot(allStarts, pulseScore, 'k-', 'DisplayName', 'Pulse Score');
legend;
title('Spectrogram Filtered Scores');

%% Transform to continuous sequence in time - add mutual exclusion and thresholding

% initialize continous weights
pulseWeight = zeros(length(sSequence),1);
sineWeight = zeros(length(sSequence),1);
pulseCtr = zeros(length(sSequence),1);
sineCtr = zeros(length(sSequence),1);

% transform binned weights to continous weights
for i = 1:length(allStarts)
    pulseWeight(allStarts(i)+1:allEnds(i)) = pulseWeight(allStarts(i)+1:allEnds(i))+pulseScore(i);
    sineWeight(allStarts(i)+1:allEnds(i)) = sineWeight(allStarts(i)+1:allEnds(i))+sineScore(i);

    pulseCtr(allStarts(i)+1:allEnds(i)) = pulseCtr(allStarts(i)+1:allEnds(i))+1;
    sineCtr(allStarts(i)+1:allEnds(i)) = sineCtr(allStarts(i)+1:allEnds(i))+1;
end

% Normalize by number of frames input
pulseWeight = pulseWeight./pulseCtr;
sineWeight = sineWeight./sineCtr;

% Apply mutual exclusion
applyExclusion = @(x, y) (x .* (x > y));
pulseWeight = applyExclusion(pulseWeight, sineWeight);
sineWeight = applyExclusion(sineWeight, pulseWeight);

% Threshold signals
threshold = 0.75;
pulses = pulseWeight > threshold;
sines = sineWeight > threshold;

% Plot pure song sequences
figure;
plot(t,sSequence, 'k');
hold on;

maskedPulses = sSequence; maskedPulses(~pulses) = NaN; 
maskedSines = sSequence; maskedSines(~sines) = NaN; 

plot(t,maskedPulses, 'b', 'DisplayName', 'Pulse');
plot(t,maskedSines, 'r', 'DisplayName', 'Sine');
legend;
title('Classified Song Signals');

%% Envelope analysis

% Compute envelope
[ub,lb] = envelope(sSequence,300,'peak');

% Get baseline values of bounds
baselineUB = mean(ub);
baselineLB = mean(lb);

% Select all sound bouts exceeding threshold
thresh = 0.01;
contUB = regionprops(ub>thresh,'PixelList');
contLB = regionprops(lb<-thresh,'PixelList');

% Extend peaks to valleys
ubEnvelopeIndices = [];
for i = 1:length(contUB)
    if length(contUB(i).PixelList)>1e3
        ubEnvelopeIndices = cat(1,ubEnvelopeIndices,contUB(i).PixelList(:,2));
    end
end

lbEnvelopeIndices = [];
for i = 1:length(contLB)
    if length(contLB(i).PixelList)>1e3
        lbEnvelopeIndices = cat(1,lbEnvelopeIndices,contLB(i).PixelList(:,2));
    end
end

% Combine
allEnvelopes = unique(cat(1,lbEnvelopeIndices,ubEnvelopeIndices));

% Compute fraction of envelopes containing pure pulse song. 
fractionPurePulse = mean(pulses(allEnvelopes));
disp([num2str(fractionPurePulse*100),'% of sound bouts consist of pure song'])