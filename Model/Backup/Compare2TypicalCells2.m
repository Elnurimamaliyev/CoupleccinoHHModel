function [Compare,PreCellFeatures, PostCellFeatures] = Compare2TypicalCells2(PreCellFeatures, PostCellFeatures)

%% Calculate logistic regression of threshold
stimulus_vals = (0.1:0.1:1.2);
probability_factor = 0.5;
PreCellFeatures.intensity_threshold = LogregSpikeThreshold2_Typical(PreCellFeatures.spike_threshold, stimulus_vals, probability_factor);
PostCellFeatures.intensity_threshold = LogregSpikeThreshold2_Typical(PostCellFeatures.spike_threshold, stimulus_vals, probability_factor);

%% Calculate Median, mean, standarddeviation
% Amplitude
PreCellFeatures.amplitude.mean = mean(PreCellFeatures.amplitude.amplitude);
PreCellFeatures.amplitude.std = std(PreCellFeatures.amplitude.amplitude);
PreCellFeatures.amplitude.median = median(PreCellFeatures.amplitude.amplitude);
PostCellFeatures.amplitude.mean = mean(PostCellFeatures.amplitude.amplitude);
PostCellFeatures.amplitude.std = std(PostCellFeatures.amplitude.amplitude);
PostCellFeatures.amplitude.median = median(PostCellFeatures.amplitude.amplitude);

% Spike count
PreCellFeatures.spike_count.mean = mean(PreCellFeatures.spike_count.spike_count);
PreCellFeatures.spike_count.std = std(PreCellFeatures.spike_count.spike_count);
PreCellFeatures.spike_count.median = median(PreCellFeatures.spike_count.spike_count);
PostCellFeatures.spike_count.mean = mean(PostCellFeatures.spike_count.spike_count);
PostCellFeatures.spike_count.std = std(PostCellFeatures.spike_count.spike_count);
PostCellFeatures.spike_count.median = median(PostCellFeatures.spike_count.spike_count);

% Spike latency
PreCellFeatures.spike_latency.mean = mean(PreCellFeatures.spike_latency.spike_latency);
PreCellFeatures.spike_latency.std = std(PreCellFeatures.spike_latency.spike_latency);
PreCellFeatures.spike_latency.median = median(PreCellFeatures.spike_latency.spike_latency);
PostCellFeatures.spike_latency.mean = mean(PostCellFeatures.spike_latency.spike_latency);
PostCellFeatures.spike_latency.std = std(PostCellFeatures.spike_latency.spike_latency);
PostCellFeatures.spike_latency.median = median(PostCellFeatures.spike_latency.spike_latency);

% Inter-Spike Interval
PreCellFeatures.inter_spike_interval.mean = mean(PreCellFeatures.inter_spike_interval.inter_spike_interval);
PreCellFeatures.inter_spike_interval.std = std(PreCellFeatures.inter_spike_interval.inter_spike_interval);
PreCellFeatures.inter_spike_interval.median = median(PreCellFeatures.inter_spike_interval.inter_spike_interval);
PostCellFeatures.inter_spike_interval.mean = mean(PostCellFeatures.inter_spike_interval.inter_spike_interval);
PostCellFeatures.inter_spike_interval.std = std(PostCellFeatures.inter_spike_interval.inter_spike_interval);
PostCellFeatures.inter_spike_interval.median = median(PostCellFeatures.inter_spike_interval.inter_spike_interval);

%% Compare both cells mathematically (ratio/ difference)
% each trial 
Compare.trials.amplitude = PostCellFeatures.amplitude.amplitude ./ PreCellFeatures.amplitude.amplitude;
Compare.trials.spike_count = PostCellFeatures.spike_count.spike_count ./ PreCellFeatures.spike_count.spike_count;
Compare.trials.spike_latency = PostCellFeatures.spike_latency.spike_latency - PreCellFeatures.spike_latency.spike_latency;
Compare.trials.inter_spike_interval = PostCellFeatures.inter_spike_interval.inter_spike_interval ./ PreCellFeatures.inter_spike_interval.inter_spike_interval;
Compare.trials.spike_threshold = PostCellFeatures.spike_threshold ./ PreCellFeatures.spike_threshold;

% means
Compare.mean.amplitude = PostCellFeatures.amplitude.mean ./ PreCellFeatures.amplitude.mean;
Compare.mean.spike_count = PostCellFeatures.spike_count.mean ./ PreCellFeatures.spike_count.mean;
Compare.mean.spike_latency = PostCellFeatures.spike_latency.mean - PreCellFeatures.spike_latency.mean;
Compare.mean.inter_spike_interval = PostCellFeatures.inter_spike_interval.mean ./ PreCellFeatures.inter_spike_interval.mean;
Compare.mean.intensity_threshold = PostCellFeatures.intensity_threshold ./ PreCellFeatures.intensity_threshold;

%% End
end