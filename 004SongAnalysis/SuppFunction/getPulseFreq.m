%% functions related to calculating carrier frequencies
% functions from Clemens, Coen, Roemscheid, et al. Current Biology. 2018
% (doi: 10.1016/j.cub.2018.06.011)

function [Fpeak, amp, F] = getPulseFreq(pulses,Fs,threshold)
% [Fpeak, amp, F] = getPulseFreq(pulses, threshold=1/e)
if ~exist('threshold','var')
   threshold = 1/exp(1);
end
fftLen = Fs/5; 
tmp = abs(fft(pulses, fftLen))';
amp = tmp(:,1:fftLen/2);
F = (1:ceil(fftLen/2)) * Fs/fftLen;
Fpeak = centerOfMass(F, amp, threshold)';
end


function com = centerOfMass(x,y, threshold)
% computes center of mass
% USAGE: centerOfMass(x,y, threshold)

if nargin==2
   threshold =  0.5;
end
y = normalizeMax(y);
y = limit(y-threshold,0,inf);
com = x*normalizeSum(y)';
end

function x = limit(x, varargin)
%x = limit(x, [low=0], [up=1])
if nargin==1
   up = 1;
   low = 0;
else
   low = varargin{1};
   up = varargin{2};
end
x = max(x,low); % Bound elements from below, x >= lowerBound
x = min(x,up); % Bound elements from above, x <= upperBound 
end

function data = normalizeSum(data)
%normalize data, such sum(data)=1
%   data = normalizeSum(data)
%
%args
%   data
%returns
%   normalized data
%
% created 20130209 jan
if (size(data,1)==1 || size(data,2)==1)
   data = data/nansum(data);
else
   for cat = 1:size(data,1)
      data(cat,:) = data(cat,:)/(nansum(data(cat,:)));
   end
end
end

function data = normalizeMax(data)
%NORMALIZE DATA, such that max=1, but preserves offset (min)
%   data = normalizeMax(data)
%
%ARGS
%   data
%RETURNS
%   normalized data
%
% created 06/09/19 Jan
if (size(data,1)==1 | size(data,2)==1)
   data = data/(nanmax(data)+eps);
else
   for cat = 1:size(data,1)
      data(cat,:) = data(cat,:)/(nanmax(data(cat,:))+eps);
   end
end
end