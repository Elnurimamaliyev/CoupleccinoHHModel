%% Compare with and without pharma

%% Preparation
% start clean
close all
clear
clc
folder = "/Users/maren/Library/CloudStorage/OneDrive-CarlvonOssietzkyUniversitätOldenburg/Uni/Master/2. FS/Team Project/ClassicalStatistics/Conditions/Results";

%% load data
% load('240619-B-M-1.mat')
% name = "240619-B-M-1";

load('20240619_First_Complete_Trial_-_Blind_Pharma');
name = "20240619-S-AP-1";

% load("20240624-S-A-1.mat")
% name = "20240624-S-A-1";



%% settings
binsize = 0.1;

% save datasets
Pre_data_control_all = pre_cell_10;
Post_data_control_all = post_cell_10;

Pre_data_pharma_all = [pre_cell_pharma_10_2; pre_cell_pharma_10_3; pre_cell_pharma_10_4; pre_cell_pharma_10_5];
Post_data_pharma_all = [post_cell_pharma_10_2; post_cell_pharma_10_3; post_cell_pharma_10_4; post_cell_pharma_10_5];
%Pre_data_pharma_all = pre_cell_pharma_60;
%Post_data_pharma_all = post_cell_pharma_60;
%Pre_data_pharma_all = pre_cell_pharma_25;
%Post_data_pharma_all = post_cell_pharma_25;

%% Only Good trials
trials_includedP = 1:40;
trials_includedC = 1:10;
%trials_includedP = trials_included_postpharm;
%trials_includedC = trials_included_prepharm;

Pre_data_control_all = Pre_data_control_all(trials_includedC,:);
Post_data_control_all = Post_data_control_all(trials_includedC,:);
Pre_data_pharma_all = Pre_data_pharma_all(trials_includedP,:);
Post_data_pharma_all = Post_data_pharma_all(trials_includedP,:);

%% Decide which trials to take
% Maybe adjust these values
MinPeakProminence = 12;
MaxPeakWidth = 10000;
MinPeakHeight = -60;

% Good trials control 
[Pre_data_control, Post_data_control, trialsC] = TrialDecision3(Pre_data_control_all, Post_data_control_all, trials_includedC, binsize, MinPeakProminence, MaxPeakWidth, MinPeakHeight, name, "Control", folder);

% Good trials with pharma
[Pre_data_pharma, Post_data_pharma, trialsP] = TrialDecision3(Pre_data_pharma_all, Post_data_pharma_all, trials_includedP, binsize, MinPeakProminence, MaxPeakWidth, MinPeakHeight, name, "Pharma", folder);

close all
%% Get feature comparison of both cells for each condition
% Get feature comparison for trials without pharmacology 
[ControlDataFeatures, Pre_data_control, Post_data_control] = Compare2Cells2(Pre_data_control, Post_data_control, "Control", folder, name);

% Get feature comparison for trials with pharmacology 
[PharmaDataFeatures, Pre_data_pharma, Post_data_pharma] = Compare2Cells2(Pre_data_pharma, Post_data_pharma, "Pharma", folder, name);

% Plot the different features for the Pre and Post cell
PlotPrePost(Pre_data_control, Post_data_control, Pre_data_pharma, Post_data_pharma, name, folder)

%% Compare both conditions visually over trials

figure("WindowState","maximized")
Control = sprintf("Control (trials: %1.0f to %1.0f)", trialsC(1), trialsC(end));
Pharma = sprintf("Pharma (trials: %1.0f to %1.0f)", trialsP(1), trialsP(end));

% amplitude
x1 = 1:10;
x2 = 1:10;
num_trials = height(ControlDataFeatures.trials.amplitude);
subplot(5,1,1)
plot(x1,ControlDataFeatures.trials.amplitude, 'b+','MarkerSize',10,'LineWidth', 1.5)
hold on
plot(x2, PharmaDataFeatures.trials.amplitude, 'rx','MarkerSize',10,'LineWidth', 1.5)
miny = [min(ControlDataFeatures.trials.amplitude(:)), min(PharmaDataFeatures.trials.amplitude(:))];
maxy = [max(ControlDataFeatures.trials.amplitude(:)), max(PharmaDataFeatures.trials.amplitude(:))];
xlim([0.5 num_trials+0.5])
ylim([min(miny)-0.1 max(maxy)+0.1])
title('Hyperpolarization Amplitude ratio (Post/Pre) over trials')
xlabel('trials')
ylabel('Amplitude ratio')
legend(Control, Pharma, 'Location','eastoutside')

% spike count
subplot(5,1,2)
plot(x1, ControlDataFeatures.trials.spike_count(:,1), 'bo', ...
     x1, ControlDataFeatures.trials.spike_count(:,2), 'b+', ...
     x1, ControlDataFeatures.trials.spike_count(:,3), 'bx', ...
     'MarkerSize',10,'LineWidth', 2)
hold on
plot(x2, PharmaDataFeatures.trials.spike_count(:,1), 'ro', ...
     x2, PharmaDataFeatures.trials.spike_count(:,2), 'r+', ...
     x2, PharmaDataFeatures.trials.spike_count(:,3), 'rx', ...
     'MarkerSize',10,'LineWidth', 1)
maxy = [max(ControlDataFeatures.trials.spike_count(:)), max(PharmaDataFeatures.trials.spike_count(:))];
ylim([0 max(maxy)+1])
xlim([0.5 num_trials+0.5])
title('Spike count ratio (Post/Pre) over trials')
xlabel('trials')
ylabel('Spike count ratio')
legend(append(Control, ": Spontanious"), append(Control, ": stimulus intensity 1 [nA]"), append(Control, ": stimulus intensity 1.5 [nA]"), ...
    append(Pharma, ": Spontanious"), append(Pharma, ": stimulus intensity 1 [nA]"), ...
    append(Pharma, ": stimulus intensity 1.5 [nA]"),'Location','eastoutside')

% Spike latency
subplot(5,1,3)
plot(x1, ControlDataFeatures.trials.spike_latency(:,1), 'b+', ...
     x1, ControlDataFeatures.trials.spike_latency(:,2), 'bx', ...
     'MarkerSize',10,'LineWidth', 2)
hold on
plot(x2, PharmaDataFeatures.trials.spike_latency(:,1), 'r+', ...
     x2, PharmaDataFeatures.trials.spike_latency(:,2), 'rx', ...
     'MarkerSize',10,'LineWidth', 1)
maxy = [max(ControlDataFeatures.trials.spike_latency(:)), max(PharmaDataFeatures.trials.spike_latency(:))];
miny = [min(ControlDataFeatures.trials.spike_latency(:)), min(PharmaDataFeatures.trials.spike_latency(:))];
if min(miny) < 0
    ymin = min(miny);
else
    ymin = 0;
end
ylim([ymin max(maxy)+10])
xlim([0.5 num_trials+0.5])
title('Spike latency difference (Post - Pre) over trials')
xlabel('trials')
ylabel('Spike latency difference')
legend(append(Control, ": stimulus intensity 1 [nA]"), append(Control, ": stimulus intensity 1.5 [nA]"), ...
    append(Pharma, ": stimulus intensity 1 [nA]"), ...
    append(Pharma, ": stimulus intensity 1.5 [nA]"),'Location','eastoutside')

% Inter-Spike-Interval
subplot(5,1,4)
plot(x1, ControlDataFeatures.trials.inter_spike_interval(:,1), 'b+', ...
     x1, ControlDataFeatures.trials.inter_spike_interval(:,2), 'bx', ...
     'MarkerSize',10,'LineWidth', 2)
hold on
plot(x2, PharmaDataFeatures.trials.inter_spike_interval(:,1), 'r+', ...
     x2, PharmaDataFeatures.trials.inter_spike_interval(:,2), 'rx', ...
     'MarkerSize',10,'LineWidth', 1)
maxy = [max(ControlDataFeatures.trials.inter_spike_interval(:)), max(PharmaDataFeatures.trials.inter_spike_interval(:))];
miny = [min(ControlDataFeatures.trials.inter_spike_interval(:)), min(PharmaDataFeatures.trials.inter_spike_interval(:))];
ylim([min(miny)-10 max(maxy)+10])
xlim([0.5 num_trials+0.5])
title('Inter-Spike-Interval ratio (Post/Pre) over trials')
xlabel('trials')
ylabel('Inter-Spike-Interval ratio')
legend(append(Control, ": stimulus intensity 1 [nA]"), append(Control, ": stimulus intensity 1.5 [nA]"), ...
    append(Pharma, ": stimulus intensity 1 [nA]"), ...
    append(Pharma, ": stimulus intensity 1.5 [nA]"),'Location','eastoutside')

% Spike threshold
subplot(5,1,5)
plot(x1, ControlDataFeatures.trials.spike_threshold, 'b+','MarkerSize',10,'LineWidth', 1.5)
hold on
plot(x2, PharmaDataFeatures.trials.spike_threshold, 'rx','MarkerSize',10,'LineWidth', 1.5)
maxy = [max(ControlDataFeatures.trials.spike_threshold(:)), max(PharmaDataFeatures.trials.spike_threshold(:))];
miny = [min(ControlDataFeatures.trials.spike_threshold(:)), min(PharmaDataFeatures.trials.spike_threshold(:))];
ylim([min(miny)-1 max(maxy)+1])
xlim([0.5 num_trials+0.5])
title('Spike Threshold ratio (Post/Pre) over trials')
xlabel('trials')
ylabel('Spike Threshold ratio')
legend(Control, Pharma, 'Location','eastoutside')

saveas(gcf, fullfile(folder, append(name, "_ConditionComparison")))
saveas(gcf, fullfile(folder, append(name, "_ConditionComparison.svg")))

%% Save results
filename = append(name, "_Results");
save(fullfile(folder, filename), 'ControlDataFeatures', 'trialsC', 'PharmaDataFeatures', 'trialsP', ...
    'Pre_data_control', 'Post_data_control', 'Pre_data_pharma', 'Post_data_pharma')

%% end
