function [spike_count, spike_count_time] = SpikeCount(data, spike_idx, spike_height, num_trials, time_window_count, binsize)
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
%   data            vector of recorded membrane potential of all trials [mV]
%   spike_idx       indices of all detected spikes
%   spike_height    height of all detected spike amplitudes [mV]
%   num_trials      Number of conducted trials
%   time_window_count      time where the stimulus was applied 
%                          [starttime : endtime]
% Output:
%   spike_count_mean    Mean of spike number during stimulus applied in
%                       time window
%   spike_count_sd      Standard deviation of spike number during stimulus 
%                       applied in time window
%   spike_count_median  Median of spike number during stimulus applied in
%                       time window
% -----
% Example function call:
% [spike_timevec, spike_idxvec, spike_heightvec] = SpikeDetection2(MinPeakProminence, MaxPeakWidth, binsize)
% -----
% Author:   Maren Duken
% Date:     20.06.2024
% -----

WindowNumSpikeCount = height(time_window_count);
% get spike count in specific window
for i = 1:WindowNumSpikeCount
    WindowStart = time_window_count(i,1);
    WindowEnd = time_window_count(i,end);
    for ii = 1 : num_trials
        TrialSpikes = spike_idx{ii, 1};
        TrialHeights = spike_height{ii, 1};
        window_idx = find(TrialSpikes >= WindowStart & TrialSpikes <= WindowEnd);
        spike_count_idx{ii, i} = TrialSpikes(window_idx);
        spike_count_height{ii, i} = TrialHeights(window_idx);
        spike_count_time{ii, i} = spike_count_idx{ii, i}*binsize/1000;
    end
end

% visualize results to control spike detection
figure('WindowState','maximized')
f3 = tiledlayout(num_trials, WindowNumSpikeCount);
f3.Title.String = 'Spike count control';
f3.Padding = 'none';
f3.TileSpacing = 'none';
for ii = 1:num_trials
    respvec = data(ii,:);
    for i = 1:WindowNumSpikeCount
        WindowStart = time_window_count(i,1);
        WindowEnd = time_window_count(i,end);
        PlotWindow = (time_window_count(i,1) - 1000):(time_window_count(i,end) + 1000);
        nexttile
        plot(PlotWindow, respvec(1,PlotWindow ));
        hold on
        plot(spike_count_idx{ii,i}, spike_count_height{ii,i}, '*', 'Color', 'r')
        xline(WindowStart, 'Color','r')
        xline(WindowEnd, 'Color','r')
        xlim([PlotWindow(1,1) PlotWindow(1,end)])
    end
end

% calculate spike count
spike_count = zeros(num_trials,WindowNumSpikeCount);
for ii = 1:num_trials
    for i = 1:WindowNumSpikeCount
        spike_count(ii,i) = length(spike_count_idx{ii,i});
    end
end
