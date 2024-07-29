function [spike_idx, spike_height] = SpikeDetection2(respmatrix, MinPeakProminence, MaxPeakWidth, MinPeakHeight)
% [spike_timevec, spike_idxvec, spike_heightvec] = SpikeDetection2(MinPeakProminence, MaxPeakWidth, binsize)
% This function detects the spikes for one deoplarisation in one trials of 
% an experiment.
% Spikes are detected by 
% Spike times are defined as the times when the membrane potential reaches
% the maximum value of the action potential.
% The function relies on the Matlab function findpeaks with the parameter
% MinPeakprominence. 
% The function testspikefinding is used to control if the spike detection
% worked properly.
% -----
% Input:
%   respvec         vector of recorded membrane potential of one trial [mV]
%   timewindow      time where the stimulus was applied 
%                   [starttime endtime]
%   MinPeakProminence parameter for the findpeaks function, indicating the 
%                   prominence of a spike beneath the maximal value [mV] 
%   (Peak Width)
%   stop_criteria   criterion where the spike counting should stop. Height
%                   where the stimulus gets back to resting potential [mV]
%   binsize         duration of one recording bin [ms]
% Output:
%   spike_timevec   nxm cell array, spike times for n trials of m
%                   experiments [s]
%   spike_count     Number of spikes that occured during the stimulus
%   spike_latency_time   Time difference between the stimulus onset and the
%                   first spike [s]
% -----
% Example function call:
% [spike_timevec, spike_count, spike_latency_time] = SpikeDetectionLatency2(respvec, timewindow, 20, stop_criteria, binsize)
% -----
% Author:   Maren Duken
% Date:     20.06.2024
% modified from spikedetection written by Jutta Kretzberg (21.05.2024)
% -----

%% Convert the matrix with the recordings into a cell array 
% with one cell per trial
recData = mat2cell(respmatrix,ones(1,size(respmatrix,1)));
 
%% apply findpeaks to all cells
[spike_idx, spike_height] = cellfun(@(x)helpfindpeaks(x,MinPeakProminence, MaxPeakWidth,MinPeakHeight),... 
    recData,'UniformOutput',false);

% end of function
end

function [spike_idxvec, spike_heightvec]= helpfindpeaks(respvec, MinPeakProminence, MaxPeakWidth, MinPeakHeight)
% private subfunction to make a version of findpeaks that gives back only
% the spike times in seconds (not as value and index pair) for one trial
% See function testspikefinding for details.

%% find spikes in the vector
% ignore the values, we only need the times
[spike_heightvec,index] = findpeaks(respvec,"MinPeakProminence",MinPeakProminence,"MaxPeakWidth",MaxPeakWidth, "MinPeakHeight",MinPeakHeight);

%% convert spike indices to spike times in seconds
spike_idxvec = index;

% end of subfunction helpfindpeaks
end