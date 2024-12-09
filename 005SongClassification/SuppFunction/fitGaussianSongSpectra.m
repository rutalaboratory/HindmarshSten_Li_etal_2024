function [songSpectrum, f] = fitGaussianSongSpectra(gtSongs, NDft, Fs, sampleWindow, sampleOverlap)

%%
%
% Helper function to compute the average spectrogram of a series of
% acoustic signals. Signals should be a structure with a single field
% (trace). NDft is the number of FFT points for spectrogram. Fs is sample rate.
% sampleWindow is the spectrogram window size. sampleOverlap give s the overlap between
% spectrogram windows.
%
%%

%initialize signal and frequency
s = []; f = [];

%loop through all traces
for k = 1:length(gtSongs)
    trace = gtSongs(k).trace;

    %if shorter than sample window, pad with zeros
    if length(trace) < sampleWindow
        trace(end:sampleWindow) = 0;
    end
    
    %if longer than sample window, remove tail. 
    trace = trace(1:sampleWindow);

    %compute spectrogram
    [s(:,k),f(:,k)] = spectrogram(trace,sampleWindow,sampleOverlap,NDft,Fs);
    
end
s = mean(abs(s),2);

%fit a mixture of gaussians to the pulse spectrum
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
songFit = fit(f(:,1),s, ft, opts );
songSpectrum = songFit(f(:,1));

end


