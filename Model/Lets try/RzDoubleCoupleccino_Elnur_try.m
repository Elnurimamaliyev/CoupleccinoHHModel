%% Two-compartment Leech Retzius cell Model (preliminary version)
% - MSc Course Neuroscience Team Project in Summer Semester 2024 
% - Modified by Elnur Imamaliyev from RzCell_prelim written by Go Ashida
% 
% 
% ------------------------------------------------------------------
% Author: Elnur Imamaliyev
% Date: 04.07.2024
%% Clean and clear workspace
clear
clc
%dataset to chose
addpath 'C:\Users\icbmadmin\Desktop\Team Project\Real Cell\Data'
addpath ('C:\Users\icbmadmin\Desktop\Team Project\Model\Functions')
load 20240628-B-M-3.mat


% simulation parameters
IampsModel = [-1.0, +1.5];  % [nA] step input amplitudes 
IampsReal =  [-1.0, +1.5];  % [nA] step input amplitudes 
% IampsModel = [-1.0, +1.0];  % [nA] step input amplitudes 
% IampsReal =  [-1.0, +1.5];  % [nA] step input amplitudes 

dt    = 0.1;                % [ms] time step (Bin size)
Tinit = 1000;               % [ms] initial silent period
Tisi  = 1000;               % [ms] inter-stimulus-interval 
Tinp  = 1000;               % [ms] duration of each step current
Ninit = round(Tinit/dt);
Nisi  = round(Tisi/dt); 
Ninp  = round(Tinp/dt); 
Nall  = Ninit + (Ninp+Nisi)*length(IampsModel); % total number of time steps
Tall  = Nall*dt/1000;       % [sec] total simulation time
t     = (0:Nall-1)*dt/1000; % [sec] time vector

% make current input
IinjModel = zeros(1,Nall);           % input vector
IinjReal = zeros(1,Nall);           % input vector

for i = 1:length(IampsModel)         
    istart = Ninit+(Ninp+Nisi)*(i-1); % stating index for each input

    IinjModel(istart+(1:Ninp)) = 1000*IampsModel(i); % [nA]->[pA]
    IinjReal(istart+(1:Ninp)) = 1000*IampsReal(i); % [nA]->[pA]

end


% Selecting real cell data
pro_cell = pre_cell_15(7,:);
post_cell = post_cell_15(7,:);



% Setting parameters of Double Rz cell model
% Ion Channel Conductance Parameters
gKFactor = 85.5;     % default: 57
gAFactor = 249;    % default: 166
gNFactor = 800;    % default: 650
gLFactor = 18;     % default: 12
gHFactor = 1.5;      % default: 1

% Compartment Conductance Parameters
gCompFactor = 82.5;  % default: 55
gGCoupFactor = 33; % default: [21:22]

% Reversal Potential Parameters
EL = -60;          % default:-40
EH = -30;          % default:-20
EN = +60;          % default:+43
EK = -70;          % default:-50

% Capacitance Parameters
c1_factor = 178.5; % Soma Capacitance
c2_factor = 15;  % Gap Junction Capacitance
c3_factor = 30;  % Spike Initation Zone Capacitance

% calling Rz cell model
gGCoupFactor_vals = gGCoupFactor;
% 20.2:0.2:22.2;


% Run the Model and save responses
for Gi = 1:length(gGCoupFactor_vals)
    gGCoupFactor = gGCoupFactor_vals(Gi);

    [v1_Pre_vals,v2_Pre_vals,v3_Pre_vals,v1_Post_vals,v2_Post_vals,v3_Post_vals] = DoubleRz_MultiComp(IinjModel,dt, gCompFactor, gGCoupFactor, gNFactor, gKFactor, gAFactor, gLFactor, gHFactor, EL, EH, EN, EK, c1_factor, c2_factor, c3_factor); % call the model
    
    % Saving responses
    Responses.v1_Pre(Gi, :) = v1_Pre_vals;
    Responses.v2_Pre(Gi, :) = v2_Pre_vals;
    Responses.v3_Pre(Gi, :) = v3_Pre_vals;
    Responses.v1_Post(Gi, :) = v1_Post_vals;
    Responses.v2_Post(Gi, :) = v2_Post_vals;
    Responses.v3_Post(Gi, :) = v3_Post_vals;
end

% Analyzing

    

for Gi = 1:length(gGCoupFactor_vals)
    gGCoupFactor = gGCoupFactor_vals(Gi);

    % hyperpolarization amplitude (Pre & Post)
    % time window parameters for hyperpolarization amplitude
        stim_onset = 10001;
        stim_offset = 20000;
        window_length = 5000;
        
        v1_Pre = Responses.v1_Pre(Gi, :);
        v1_Post = Responses.v1_Post(Gi, :);

    % Resting Membrane potential (Vm)
        Vm_Pre = mean(v1_Pre(stim_onset-window_length:stim_onset));
        Vm_Post = mean(v1_Post(stim_onset-window_length:stim_onset));

    % Calculate hyperpolarization amplitude of Pre and Post cell
        hyperpolarized_membrane_potential_Pre = mean(v1_Pre(stim_offset-window_length:stim_offset));
        hyperpolarized_membrane_potential_Post = mean(v1_Post(stim_offset-window_length:stim_offset));
        
        amplitude_Pre = Vm_Pre - hyperpolarized_membrane_potential_Pre;
        amplitude_Post = Vm_Post - hyperpolarized_membrane_potential_Post;

    % Spike Detection and Count
        MinPeakProminence = 3;
        MaxPeakWidth = 10000;
        MinPeakHeight = -50;

        [spike_idx_Pre, spike_height_Pre] = SpikeDetection2(v1_Pre, MinPeakProminence, MaxPeakWidth, MinPeakHeight); 
        [spike_idx_Post, spike_height_Post] = SpikeDetection2(v1_Post, MinPeakProminence, MaxPeakWidth, MinPeakHeight); 

    % spike count
        time_window_count = 30001:40000; 
        [spike_count_Pre, spike_count_time_Pre] = SpikeCount(v1_Pre, spike_idx_Pre, spike_height_Pre, 1, time_window_count, dt);
        [spike_count_Post, spike_count_time_Post] = SpikeCount(v1_Post, spike_idx_Post, spike_height_Post, 1, time_window_count, dt);

        close
        close

    % Spike latency
        if spike_count_Pre > 0
            % Calculation
            spike_latency_Pre  = spike_count_time_Pre{1}(1,1)  - time_window_count(1)*dt/1000;
            % Convertation 
            spike_latency_Pre = spike_latency_Pre * 1000; % convert to ms
        else
            spike_latency_Pre = 0;
        end

        if spike_count_Post > 0
            spike_latency_Post = spike_count_time_Post{1}(1,1) - time_window_count(1)*dt/1000;   
            spike_latency_Post = spike_latency_Post * 1000; % convert to ms
        
        else
            spike_latency_Post = 0;            
        end



    % Inter-Spike Interval        

        if spike_count_Pre > 0
            spike_count_Pre  = length(spike_count_time_Pre{1});
            total_time_Pre   = spike_count_time_Pre{1}(end) - spike_count_time_Pre{1}(1);
            isi_Pre = total_time_Pre/(spike_count_Pre-1);
            % Convertation
            inter_spike_interval_Pre = isi_Pre   * 1000; % convert to ms
        else
            inter_spike_interval_Pre = 0;
        end


        if spike_count_Post > 0
           spike_count_Post = length(spike_count_time_Post{1});
           total_time_Post  = spike_count_time_Post{1}(end) - spike_count_time_Post{1}(1);
           isi_Post = total_time_Post/(spike_count_Post-1);
           % Convertation
           inter_spike_interval_Post = isi_Post * 1000; % convert to ms
        else
           inter_spike_interval_Post = 0;
        end
        
    
    % plotting
    % color settings
    c0 = [0,0.6,0]; 
    c1 = [0,0,0.6];
    c2 = [0.6,0,0];

    % open fig 
    figure; clf; set(gcf,'Position',[50,50,1000,750]);

    % Presynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,1); cla; hold on; 
    % plot(t,v1_Pre,'-',color=c1); text(0.5,25,  sprintf("Spike count: %.0f", spike_count_Pre ), color=c1,fontsize=14);
    plot(t,Responses.v1_Pre,'-',color=c1); 

    txt_pre = "Spike count: " + num2str(spike_count_Pre);
    text(0.5,40,  txt_pre, color=c1, fontsize=10);

    txt_lat_pre = "Spike Latency: " + num2str(spike_latency_Pre);
    text(0.5,31,  txt_lat_pre, color=c1, fontsize=10);

    txt_hyper_pre = "Hyper. Amp.: " + num2str(amplitude_Pre);
    text(0.5,22,  txt_hyper_pre, color=c1, fontsize=10);

    txt_isi_pre = "ISI: " + num2str(inter_spike_interval_Pre);
    text(0.5,13,  txt_isi_pre, color=c1, fontsize=10);

    txt_gA = "gA: " + num2str(gAFactor);
    txt_gK = " gK: " + num2str(gKFactor);

    text(0.5,4,  txt_gA+txt_gK, color=c1, fontsize=10);

    title('Membrane Potential of Presynaptic cell Soma')

    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Presynaptic compartment 2 (Gap)
    subplot(3,2,3); cla; hold on; 
    title('Membrane Potential of Pressynaptic cell Gap')
    plot(t,Responses.v2_Pre,'-',color=c1); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 



    % Presynaptic compartment 3 (spike-initiation zone)
    subplot(3,2,5); cla; hold on; 
    plot(t,Responses.v3_Pre,'-',color=c1);
    title('Membrane Potential of Presynaptic cell SIZ')
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Postsynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,2); cla; hold on; 
    plot(t,Responses.v1_Post,'-',color=c2); 
    txt_post = "Spike count: " + num2str(spike_count_Post);
    text(0.5,40,  txt_post, color=c2, fontsize=10);

    txt_lat_post = "Spike Latency: " + num2str(spike_latency_Post);
    text(0.5,31,  txt_lat_post, color=c2, fontsize=10);

    txt_hyper_post = "Hyper. Amp.: " + num2str(amplitude_Post);
    text(0.5,22,  txt_hyper_post, color=c2, fontsize=10);

    txt_isi_post = "ISI: " + num2str(inter_spike_interval_Post);
    text(0.5,13,  txt_isi_post, color=c2, fontsize=10);

    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    title('Membrane Potential of Postsynaptic cell Soma')

    subplot(3,2,4); cla; hold on; 
    title('Membrane Potential of Postsynaptic cell Gap')
    plot(t,Responses.v2_Post,'-',color=c2); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Postsynaptic compartment 2 (Gap)
    subplot(3,2,6); cla; hold on; 
    title('Membrane Potential of Postsynaptic cell SIZ')
    plot(t,Responses.v3_Post,'-',color=c2); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Save 
    % Gi Values

    result.GCoup_values(Gi) = gGCoupFactor_vals(Gi);
    
    % Resting Membrane potential
    result.Vm_Pre(Gi) = Vm_Pre;
    result.Vm_Post(Gi) = Vm_Post;

    % Amplitude
    result.amplitude_Pre(Gi) = amplitude_Pre;
    result.amplitude_Post(Gi) = amplitude_Post;
    
    % Spike count
    result.spike_count_Pre(Gi) = spike_count_Pre;
    result.spike_count_Post(Gi) = spike_count_Post;
    
    % Spike latency
    result.spike_latency_Pre(Gi) = spike_latency_Pre;
    result.spike_latency_Post(Gi) = spike_latency_Post;
    
    % Inter-Spike Interval
    result.inter_spike_interval_Pre(Gi) = inter_spike_interval_Pre;
    result.inter_spike_interval_Post(Gi) = inter_spike_interval_Post;

end
% end % end of RzCell_prelim()

%
% color settings
ck = "k";
c0 = [0,0.6,0]; 
c1 = [0,0,0.6];
c2 = [0.6,0,0];
c3 = [0.6,0,0.6]; % Combination of real and model

% Plot Negative Responses
NegativeWindowPlot = 5001:25000;

for Gi = 1
    %Plotting Only Soma compartment of 2 Cells

    % open fig 
    figure; clf; set(gcf,'Position',[50,50,1000,750]);

    % Negative Input Current
    subplot(3,2,1); cla; hold on; 
    plot(t(NegativeWindowPlot),IinjReal(NegativeWindowPlot)/1000,'-',color=c3);
    xlim([5001/10000,25000/10000]); set(gca,'xtick',0:0.2:25000/10000);
    ylim([-1.5,2.0]); set(gca,'ytick',-2:0.5:2); 
    xlabel('time [s]'); ylabel('current [nA]'); 
    title('Stimulus (Presynaptic-Soma) ')
    % legend('Current (Real & Model)','Location', 'best')

    % Presynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,3); cla; hold on; 
    plot(t(NegativeWindowPlot),v1_Pre(NegativeWindowPlot),'-',color=c1);
    txt_res_pre = "Rm: " + num2str(result.Vm_Pre(Gi));
    text(0.6,25,  txt_res_pre, color=c1, fontsize=10);
    txt_hyper_pre = "HA: " + num2str(result.amplitude_Pre(Gi));
    text(0.6,05,  txt_hyper_pre, color=c1, fontsize=10);
    xlim([5001/10000,25000/10000]); set(gca,'xtick',0:0.2:25000/10000);
    set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]'); 
    title('Vm Presynaptic Soma (M)')

    % Postsynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,4); cla; hold on; 
    plot(t(NegativeWindowPlot),v1_Post(NegativeWindowPlot),'-',color=c1); 
    txt_res_post = "Rm: " + num2str(result.Vm_Post(Gi));
    text(0.6,25,  txt_res_post, color=c1, fontsize=10);
    txt_hyper_post = "HA: " + num2str(result.amplitude_Post(Gi));
    text(0.6,05,  txt_hyper_post, color=c1, fontsize=10);
    xlim([5001/10000,25000/10000]); set(gca,'xtick',0:0.2:25000/10000);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]');
    title('Vm Postsynaptic Soma (M)')

    % Concatenate the data if you want to plot it as one continuous plot
    VrealPreNegwindow = [pro_cell(179001:199000)];
    tRealNeg = (179001:199000)/10000;

    % Presynaptic cell Soma (soma - Actual Cell)
    subplot(3,2,5); cla; hold on; 
    plot(tRealNeg,VrealPreNegwindow,'-',color=c2);
    xlim([179001,199000]/10000); set(gca,'xtick',0:0.2:199000);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]');
    title('Vm Presynaptic Soma (R)')

    % Concatenate the data if you want to plot it as one continuous plot
    VrealPostNegwindow = [post_cell(179001:199000)];

    % Presynaptic cell Soma (soma - Actual Cell)    
    subplot(3,2,6); cla; hold on; 
    plot(tRealNeg,VrealPostNegwindow,'-',color=c2);
    xlim([179001,199000]/10000); set(gca,'xtick',0:0.2:199000);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]');
    title('Vm Postsynaptic Soma (R)')


end




%% Plot Positive responses
PositiveWindowPlot = 25001:45000; 
ck = "k";
c0 = [0,0.6,0]; 
c1 = [0,0,0.6];
c2 = [0.6,0,0];

for Gi = 10
% % Plotting Only Soma compartment of 2 Cells
    % color settings


    % open fig 
    figure; clf; set(gcf,'Position',[50,50,1000,750]);

    % Input Current
    subplot(3,2,1); cla; hold on; 
    plot(t(PositiveWindowPlot),IinjReal(PositiveWindowPlot)/1000,'-',color=c2);
    plot(t(PositiveWindowPlot),IinjModel(PositiveWindowPlot)/1000,'-',color=c1);
    xlim([25001/10000,45000/10000]); set(gca,'xtick',0:0.2:25000/10000);
    ylim([-1.5,2.0]); set(gca,'ytick',-2:0.5:2); 
    xlabel('time [s]'); ylabel('current [nA]'); 
    title('Stimulus (Presynaptic-Soma) ')
    legend('Real Current','Model Current','Location', 'south')

    % Presynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,3); cla; hold on; 
    plot(t(PositiveWindowPlot),v1_Pre(PositiveWindowPlot),'-',color=c1);
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]');
    title('Membrane Potential of Presynaptic cell Soma (Model)')
    % Text
    txt_pre = "Spike count: " + num2str(spike_count_Pre);
    text(0.5,35,  txt_pre, color=c1, fontsize=10);
    txt_lat_pre = "Spike Latency: " + num2str(spike_latency_Pre);
    text(0.5,20,  txt_lat_pre, color=c1, fontsize=10);
    txt_isi_pre = "ISI: " + num2str(inter_spike_interval_Pre);
    text(0.5,-10,  txt_isi_pre, color=c1, fontsize=10);
    % txt_gA = "gA: " + num2str(gAFactor);
    % txt_gK = " gK: " + num2str(gKFactor);
    % text(0.5,4,  txt_gA+txt_gK, color=c1, fontsize=10);


    % Postsynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,4); cla; hold on; 
    plot(t(PositiveWindowPlot),v1_Post(PositiveWindowPlot),'-',color=c1); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]');
    title('Membrane Potential of Postsynaptic cell Soma (Model)')
    % Text
    txt_post = "Spike count: " + num2str(spike_count_Post);
    text(0.5,35,  txt_post, color=c1, fontsize=10);
    txt_lat_post = "Spike Latency: " + num2str(spike_latency_Post);
    text(0.5,20,  txt_lat_post, color=c1, fontsize=10);
    txt_isi_post = "ISI: " + num2str(inter_spike_interval_Post);
    text(0.5,-10,  txt_isi_post, color=c1, fontsize=10);

    % Concatenate the data if you want to plot it as one continuous plot
    VrealPrePoswindow = [pro_cell(239001:264000)];
    tRealPos = (239001:264000)/10000;

    % Presynaptic cell Soma (soma - Actual Cell)
    subplot(3,2,5); cla; hold on; 
    plot(tRealPos,VrealPrePoswindow,'-',color=c2);
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]'); xline(2.5, '--');
    title('Membrane Potential of Presynaptic cell Soma (Real)')

    % Concatenate the data if you want to plot it as one continuous plot
    VrealPostPoswindow = [post_cell(239001:264000)];

    % Presynaptic cell Soma (soma - Actual Cell)    
    subplot(3,2,6); cla; hold on; 
    plot(tRealPos,VrealPostPoswindow,'-',color=c2);
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall); 
    ylim([-105,45]); set(gca,'ytick',-100:20:40); 
    xlabel('time [s]'); ylabel('potential [mV]'); xline(2.5, '--');
    title('Membrane Potential of Postsynaptic cell Soma (Real)')


end


%%
% figure; 
% x_ticks = 20:0.2:22.2;
% 
% subplot(5,1,1);
% 
% plot(result.GCoup_values,result.Vm_Pre, '+', LineWidth=1.2); 
% hold on; plot(result.GCoup_values,result.Vm_Post, 'x', LineWidth=1.2)
% xlim([20.0,22.2]); set(gca,'xtick',x_ticks); 
% ylabel('Voltage [mV]'); ylim([min(min(result.Vm_Pre), min(result.Vm_Post)) - 0.3, max(max(result.Vm_Pre), max(result.Vm_Post)) + 1]);
% legend('Rm-pre','Rm-post', Location='best')
% 
% % title('Membrane Potential');
% 
% 
% subplot(5,1,2);
% plot(result.GCoup_values,result.amplitude_Pre, '+', LineWidth=1.2); 
% hold on; plot(result.GCoup_values,result.amplitude_Post, 'x', LineWidth=1.2)
% xlim([20.0,22.2]); set(gca,'xtick',x_ticks); 
% ylabel('Amplitude [mV]'); ylim([min(min(result.amplitude_Pre), min(result.amplitude_Post)) - 4, max(max(result.amplitude_Pre), max(result.amplitude_Post)) + 4]);
% legend('HA-pre','HA-post', Location='best')
% % title('Hyperpolarization Amplidute');
% 
% 
% subplot(5,1,3);
% plot(result.GCoup_values,result.spike_count_Pre, '+', LineWidth=1.2); 
% hold on; plot(result.GCoup_values,result.spike_count_Post, 'x', LineWidth=1.2)
% xlim([20.0,22.2]); set(gca,'xtick',x_ticks); 
% ylabel('Spike Count [mV]'); ylim([min(min(result.spike_count_Pre), min(result.spike_count_Post)) - 5, max(max(result.spike_count_Pre), max(result.spike_count_Post)) + 5]);
% legend('SC-pre','SC-post', Location='best')
% % title('Spike Count'); 
% 
% 
% subplot(5,1,4);
% plot(result.GCoup_values,result.spike_latency_Pre, '+', LineWidth=1.2); 
% hold on; plot(result.GCoup_values,result.spike_latency_Post, 'x', LineWidth=1.2)
% ylabel('Latency [ms]'); ylim([min(min(result.spike_latency_Pre), min(result.spike_latency_Post)) - 25, max(max(result.spike_latency_Pre), max(result.spike_latency_Post)) + 25]);
% xlim([20.0,22.2]); set(gca,'xtick',x_ticks); 
% legend('Latency-pre','Latency-post', Location='best')
% % title('Spike Latency')
% 
% 
% subplot(5,1,5);
% plot(result.GCoup_values,result.inter_spike_interval_Pre, '+', LineWidth=1.2); 
% hold on; plot(result.GCoup_values,result.inter_spike_interval_Post, 'x', LineWidth=1.2)
% xlim([20.0,22.2]); set(gca,'xtick',x_ticks); 
% ylabel('ISI [ms]'); ylim([min(min(result.inter_spike_interval_Pre), min(result.inter_spike_interval_Post)) - 25, max(max(result.inter_spike_interval_Pre), max(result.inter_spike_interval_Post)) + 25]);
% legend('ISI-pre','ISI-post', Location='best')
% % title('Inter-Spike Interval')


%% Feature Comparison Plot
% Create a tiled layout with 5 rows and 1 column
figure;
tiledlayout(5, 1, "TileSpacing", "none", "Padding", "none");

x_ticks = 20:0.2:22.2;

% First subplot
nexttile;
plot(result.GCoup_values, result.Vm_Pre, '+', 'LineWidth', 1.2); 
hold on; 
plot(result.GCoup_values, result.Vm_Post, 'x', 'LineWidth', 1.2);
xlim([20.0, 22.2]); 
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('Voltage [mV]'); 
ylim([min(min(result.Vm_Pre), min(result.Vm_Post)) - 0.3, max(max(result.Vm_Pre), max(result.Vm_Post)) + 1]);
legend('Rm-pre', 'Rm-post', 'Location', 'best');

% Second subplot
nexttile;
plot(result.GCoup_values, result.amplitude_Pre, '+', 'LineWidth', 1.2); 
hold on; 
plot(result.GCoup_values, result.amplitude_Post, 'x', 'LineWidth', 1.2);
xlim([20.0, 22.2]); 
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('Amplitude [mV]'); 
ylim([min(min(result.amplitude_Pre), min(result.amplitude_Post)) - 4, max(max(result.amplitude_Pre), max(result.amplitude_Post)) + 4]);
legend('HA-pre', 'HA-post', 'Location', 'best');

% Third subplot
nexttile;
plot(result.GCoup_values, result.spike_count_Pre, '+', 'LineWidth', 1.2); 
hold on; 
plot(result.GCoup_values, result.spike_count_Post, 'x', 'LineWidth', 1.2);
xlim([20.0, 22.2]); 
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('Spike Count'); 
ylim([min(min(result.spike_count_Pre), min(result.spike_count_Post)) - 5, max(max(result.spike_count_Pre), max(result.spike_count_Post)) + 5]);
legend('SC-pre', 'SC-post', 'Location', 'best');

% Fourth subplot
nexttile;
plot(result.GCoup_values, result.inter_spike_interval_Pre, '+', 'LineWidth', 1.2); 
hold on; 
plot(result.GCoup_values, result.inter_spike_interval_Post, 'x', 'LineWidth', 1.2);
xlim([20.0, 22.2]); 
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('ISI [ms]'); ylim([min(min(result.inter_spike_interval_Pre), min(result.inter_spike_interval_Post)) - 25, max(max(result.inter_spike_interval_Pre), max(result.inter_spike_interval_Post)) + 25]);
legend('ISI-pre', 'ISI-post', 'Location', 'best');

% Fifth subplot
nexttile;
plot(result.GCoup_values, result.spike_latency_Pre, '+', 'LineWidth', 1.2); 
hold on; 
plot(result.GCoup_values, result.spike_latency_Post, 'x', 'LineWidth', 1.2);
xlim([20.0, 22.2]); 
set(gca, 'xtick', x_ticks);  % Keep x-ticks and labels for the bottom subplot
ylabel('Latency [ms]');  ylim([min(min(result.spike_latency_Pre), min(result.spike_latency_Post)) - 25, max(max(result.spike_latency_Pre), max(result.spike_latency_Post)) + 25]);
xlabel('Coupling Strength - gGCoup [mS]');
legend('Latency-pre', 'Latency-post', 'Location', 'best');

% Add a shared title
sgtitle('Feature Comparison');


%% Mathematical comparisson of pre and post cells (ratio and difference)
% Pre and post
    CompareModelCells.GCoup_values = result.GCoup_values;
    CompareModelCells.Vm = result.Vm_Post ./ result.Vm_Pre;
    CompareModelCells.amplitude = result.amplitude_Post ./ result.amplitude_Pre;
    CompareModelCells.spike_count = result.spike_count_Post ./ result.spike_count_Pre;
    CompareModelCells.inter_spike_interval = result.inter_spike_interval_Post ./ result.inter_spike_interval_Pre;
    CompareModelCells.spike_latency = result.spike_latency_Post - result.spike_latency_Pre;

%% Matematical Comparison Plot

figure;
% Create a tiled layout with 5 rows and 1 column
tiledlayout(5, 1, "TileSpacing", "none", "Padding", "none");

x_ticks = 20:0.2:22.2;

% First subplot
nexttile;
plot(CompareModelCells.GCoup_values, CompareModelCells.Vm, 'o', 'LineWidth', 1.2);
xlim([20.0, 22.2]);
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('Voltage (Ratio)');
ylim([0, 2]); set(gca,'ytick',0.0:0.25:1.75)
legend('Vm (Post/Pre)', 'Location', 'best');

% Second subplot
nexttile;
plot(CompareModelCells.GCoup_values, CompareModelCells.amplitude, 'o', 'LineWidth', 1.2);
xlim([20.0, 22.2]);
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('Amplitude (Ratio)'); ylim([0, 2]); set(gca,'ytick',0.0:0.25:1.75)
legend('Amplitude (Post/Pre)', 'Location', 'best');

% Third subplot
nexttile;
plot(CompareModelCells.GCoup_values, CompareModelCells.spike_count, 'o', 'LineWidth', 1.2);
xlim([20.0, 22.2]);
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylabel('Spike Count (Ratio)');
ylim([0, 2]); set(gca,'ytick',0.0:0.25:1.75)
legend('Spike Count (Post/Pre)', 'Location', 'best');

% Fourth subplot
nexttile;
plot(CompareModelCells.GCoup_values, CompareModelCells.inter_spike_interval, 'o', 'LineWidth', 1.2);
xlim([20.0, 22.2]);
ylabel('ISI (Ratio)');
set(gca, 'xtick', []);  % Remove x-ticks and labels
ylim([0, 2]); set(gca,'ytick',0.0:0.25:1.75);
legend('ISI (Post/Pre)', 'Location', 'best');

% Fifth subplot
nexttile;
plot(CompareModelCells.GCoup_values, CompareModelCells.spike_latency, 'o', 'LineWidth', 1.2);
xlim([20.0, 22.2]);
set(gca, 'xtick', x_ticks);  % Keep x-ticks and labels for the bottom subplot
ylabel('Latency Difference [ms]');
ylim([min(CompareModelCells.spike_latency) - 25, max(CompareModelCells.spike_latency) + 25]);
legend('Spike Latency (Post-Pre)', 'Location', 'best');
xlabel('Coupling Strength - gGCoup [mS]');

% Add a shared title
sgtitle('Ratios and Differences');



%%



%% function [V1_Pre,V2_Pre,V3_Pre, V1_Post,V2_Post,V3_Post] = RzDouble_Multi(Iinj,dt)
function [V1_Pre,V2_Pre,V3_Pre, V1_Post,V2_Post,V3_Post] = DoubleRz_MultiComp(Iinj,dt, gCompFactor, gGCoupFactor, gNFactor, gKFactor, gAFactor, gLFactor, gHFactor, EL,EH, EN, EK, c1_factor, c2_factor, c3_factor)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Two-compartment leech Retzius cell model 
    %
    % - Input:
    %  Iinj       vector for injected current input [pA]
    %  dt         time step [ms]
    %
    % - Output:
    %  V1_Pre     membrane potential of the Presynaptic soma [mV]
    %  V2_Pre     membrane potential of the Presynaptic Gap Junction (GJ) [mV]
    %  V3_Pre     membrane potential of the Presynaptic spike-initiation zone (SIZ) [mV]
    %  V1_Post    membrane potential of the Presynaptic soma [mV]
    %  V2_Post    membrane potential of the Presynaptic Gap Junction (GJ) [mV]
    %  V3_Post    membrane potential of the Presynaptic spike-initiation zone (SIZ) [mV]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% assigning membrane parameters
    
    % Pre capacitances
    c1_Pre = c1_factor * 1e-6;  % [uF] membrane capacitance of Presynaptic soma: default 119 [pF]
    c2_Pre = c2_factor * 1e-6;  % [uF] membrane capacitance of Presynaptic Gap Compartment: default 10 [pF]
    c3_Pre = c3_factor * 1e-6;  % [uF] membrane capacitance of Presynaptic SIZ: default 20 [pF]

    % Post capacitances
    c1_Post = c1_factor * 1e-6; % [uF] membrane capacitance of Postsynaptic soma: default 119 [pF]
    c2_Post = c2_factor * 1e-6; % [uF] membrane capacitance of Postsynaptic Gap Compartment: default 10 [pF]
    c3_Post = c3_factor * 1e-6; % [uF] membrane capacitance of Postsynaptic SIZ: default 20 [pF]

    % Soma Ion Channels conductance 
    gL =  gLFactor * 1e-6; % [mS] leak conductance of soma: default 8 [nS]
    gH =  gHFactor * 1e-6; % [mS] Ih conductance of soma:   default 4 [nS]   % To simulate a neuron without a sag, you can set gH=0 and increase gL to match the measured input resistance. 

    % SIZ Ion Channels conductance 
    gN = gNFactor * 1e-6;  % [mS] Na conductance of SIZ: default 1800 [nS]   % (!) Increase it for having spike with lower intensities. Caution, it's cursed.
    gK = gKFactor * 1e-6;  % [mS] K conductance of SIZ:   default 180 [nS]  
    gA = gAFactor * 1e-6;  % [mS] KA conductance of SIZ:  default 180 [nS]

    % Compartment conductance 
    gComp =  gCompFactor * 1e-6; % [mS] conductance between the two compartments: default 80 [nS]
    
    % Gap Junction conductance
    gGCoup = gGCoupFactor * 1e-6; % [mS] conductance of Gap junction between two Cells: default 80 [nS]

    % % reversal potentials
    % EL = -25; % [mV] leak reversal potential (soma only)  Default: -25
    % EH = -25; % [mV] Ih reversal potential (soma only)    Default: -25
    % EN = +50; % [mV] Na reversal potential (SIZ only)     Default: +60
    % EK = -55; % [mV] K reversal potential (SIZ only)      Default: -60

    %% data vectors and initial values
    Ntotal = length(Iinj); % total number of data points

    % Initiate Pre cell membrane potential matrix
    V1_Pre = zeros(1,Ntotal); % [mV] membrane potential of somatic compartment of Presynaptic cell
    V2_Pre = zeros(1,Ntotal); % [mV] membrane potential of Gap Junction compartment of Presynaptic cell
    V3_Pre = zeros(1,Ntotal); % [mV] membrane potential of spike-initiation zone of Presynaptic cell
    
    % % Initiate Post cell membrane potential matrix
    V1_Post = zeros(1,Ntotal); % [mV] membrane potential of somatic compartment of Postsynaptic cell
    V2_Post = zeros(1,Ntotal); % [mV] membrane potential of Gap Junction compartment of Postsynaptic cell
    V3_Post = zeros(1,Ntotal); % [mV] membrane potential of spike-initiation zone of Postsynaptic cell

    % initial values of Pre and Post Cell 
    Vinit1 = -50; % [mV]
    Vinit2 = -52; % [mV]
    Vinit3 = -52; % [mV]

    % First Membrane potential values (Pre)
    V1_Pre(1) = Vinit1; 
    V2_Pre(1) = Vinit2; 
    V3_Pre(1) = Vinit3; 
    % First Activation/inactivation variables (Pre)
    z_Pre = infZ(Vinit1); % Ih activation variable (soma only)
    m_Pre = infM(Vinit3); % Na activation variable (SIZ only)
    h_Pre = infH(Vinit3); % Na inactivation variable (SIZ only)
    n_Pre = infN(Vinit3); % K activation variable (SIZ only)
    b_Pre = infB(Vinit3); % KA inactivation variable (SIZ only)

    % First Membrane potential values (Post)
    V1_Post(1) = Vinit1; 
    V2_Post(1) = Vinit2; 
    V3_Post(1) = Vinit3; 
    % First Activation/inactivation variables (Post)
    z_Post = infZ(Vinit1); % Ih activation variable (soma only)
    m_Post = infM(Vinit3); % Na activation variable (SIZ only)
    h_Post = infH(Vinit3); % Na inactivation variable (SIZ only)
    n_Post = infN(Vinit3); % K activation variable (SIZ only)
    b_Post = infB(Vinit3); % KA inactivation variable (SIZ only)

    %% calculate membrane response step-by-step 
    for j = 1:Ntotal-1

        % ionic currents (Gap Junctions): gGCoup[mS] * (V1[mV] - V2[mV]) = IGJ[uA]
        IGCoup = gGCoup     * ((V2_Pre(j)-V2_Post(j))) ; % current from Gap Junction

        %% Pre
        % ionic currents (soma): g[mS] * V[mV] = I[uA]
        IL_Pre = gL         * (EL-V1_Pre(j));        % leak current
        IH_Pre = gH * z_Pre * (EH-V1_Pre(j));        % Ih current

        % ionic currents (Compartments)
        IC2_Pre = gComp     * (V1_Pre(j)-V2_Pre(j)); % current from compartment 2 to compartment 1
        IC3_Pre = gComp     * (V2_Pre(j)-V3_Pre(j)); % current from compartment 3 to compartment 2

        % ionic currents (spike-initiation zone)
        IN_Pre = gN * m_Pre^4 * h_Pre * (EN-V3_Pre(j)); % Na current
        IK_Pre = gK * n_Pre^2     * (EK-V3_Pre(j));     % K current
        IA_Pre = gA * b_Pre       * (EK-V3_Pre(j));     % KA current

        % derivatives: I[uA] / C[uF] * dt[ms] = dv[mV]
        dv1_dt_Pre = ( IL_Pre + IH_Pre - IC2_Pre + Iinj(j)*1e-6) / c1_Pre; 
        dv2_dt_Pre = ( IC2_Pre - IC3_Pre - IGCoup) / c2_Pre; % -IGap
        dv3_dt_Pre = ( IN_Pre + IK_Pre + IA_Pre + IC3_Pre) / c3_Pre; 

        dz_dt_Pre = ( infZ(V1_Pre(j)) - z_Pre ) / tauZ(V1_Pre(j));
        dh_dt_Pre = ( infH(V3_Pre(j)) - h_Pre ) / tauH(V2_Pre(j));
        dn_dt_Pre = ( infN(V3_Pre(j)) - n_Pre ) / tauN(V2_Pre(j));
        db_dt_Pre = ( infB(V3_Pre(j)) - b_Pre ) / tauB(V2_Pre(j));

        %% Post
        % ionic currents (soma): g[mS] * V[mV] = I[uA]
        IL_Post = gL     * (EL-V1_Post(j));         % leak current
        IH_Post = gH * z_Post * (EH-V1_Post(j));    % Ih current

        % ionic currents (Compartments)
        IC2_Post = gComp     * (V1_Post(j)-V2_Post(j)); % current from compartment 2 to compartment 1
        IC3_Post = gComp     * (V2_Post(j)-V3_Post(j)); % current from compartment 3 to compartment 2

        % ionic currents (spike-initiation zone)
        IN_Post = gN * m_Post^4 * h_Post * (EN-V3_Post(j));  % Na current
        IK_Post = gK * n_Post^2          * (EK-V3_Post(j));  % K current
        IA_Post = gA * b_Post            * (EK-V3_Post(j));  % KA current

        % derivatives: I[uA] / C[uF] * dt[ms] = dv[mV]
        dv1_dt_Post = ( IL_Post + IH_Post - IC2_Post) / c1_Post;  
        dv2_dt_Post = ( IC2_Post - IC3_Post + IGCoup) / c2_Post; % +IGap
        dv3_dt_Post = ( IN_Post + IK_Post + IA_Post + IC3_Post) / c3_Post; 

        dz_dt_Post = ( infZ(V1_Post(j)) - z_Post ) / tauZ(V1_Post(j));
        dh_dt_Post = ( infH(V2_Post(j)) - h_Post ) / tauH(V2_Post(j));
        dn_dt_Post = ( infN(V2_Post(j)) - n_Post ) / tauN(V2_Post(j));
        db_dt_Post = ( infB(V2_Post(j)) - b_Post ) / tauB(V2_Post(j));


        %% Pre
        % Presynaptic Soma - calculate next step 
        V1_Pre(j+1) = V1_Pre(j) + dv1_dt_Pre * dt; 

        % Presynaptic Gap Junction - calculate next step 
        V2_Pre(j+1) = V2_Pre(j) + dv2_dt_Pre * dt; 

        % Presynaptic SIZ - calculate next step 
        V3_Pre(j+1) = V3_Pre(j) + dv3_dt_Pre * dt; 

        % Update the Parameters
        m_Pre = infM(V3_Pre(j)); % Na activation is assumed to be instantaneous
        z_Pre = z_Pre + dz_dt_Pre * dt; 
        h_Pre = h_Pre + dh_dt_Pre * dt; 
        n_Pre = n_Pre + dn_dt_Pre * dt; 
        b_Pre = b_Pre + db_dt_Pre * dt; 

        %% Post
        % Postsynaptic Soma - calculate next step 
        V1_Post(j+1) = V1_Post(j) + dv1_dt_Post * dt; 

        % Postsynaptic Gap Junction - calculate next step 
        V2_Post(j+1) = V2_Post(j) + dv2_dt_Post * dt; 

        % Postsynaptic SIZ - calculate next step 
        V3_Post(j+1) = V3_Post(j) + dv3_dt_Post * dt; 
        
        % Update the Parameters
        m_Post = infM(V2_Post(j)); % Na activation is assumed to be instantaneous
        z_Post = z_Post + dz_dt_Post * dt; 
        h_Post = h_Post + dh_dt_Post * dt; 
        n_Post = n_Post + dn_dt_Post * dt; 
        b_Post = b_Post + db_dt_Post * dt; 

    end

end % end of RzDouble()

%% activation/inactivation functions
function x = infM(v) % steady-state function for Na activation 
  x = 1 ./ ( 1 + exp(-(v+28)/8) ); 
end
function x = infH(v) % steady-state function for Na inactivation 
  x = 1 ./ ( 1 + exp(+(v+45)/10) ); 
end

function x = infN(v) % steady-state function for K activation 
  x = 1 ./ ( 1 + exp(-(v+35)/10) ); 
end
function x = infB(v) % steady-state function for KA inactivation 
  x = 1 ./ ( 1 + exp(+(v+55)/10) ); 
end
function x = infZ(v) % steady-state function for Ih activation 
  x = 1 ./ ( 1 + exp(+(v+55)/10) ); 
end

function x = tauH(v) % time scale [ms] for Na inactivation Default: 4.0     Changed to: 4.0
  x = 4.0 * ones(size(v)); 
end
function x = tauN(v) % time scale [ms] for K activation    Default: 20.0     Changed to: 19.55
  x = 19.55 * ones(size(v)); 
end
function x = tauB(v) % time scale [ms] for KA inactivation Default: 20.0    Changed to: 20.1
  x = 20.1 * ones(size(v)); 
end
function x = tauZ(v) % time scale [ms] for Ih activation   Default: 70.0
  x = 70 * ones(size(v)); 
end
