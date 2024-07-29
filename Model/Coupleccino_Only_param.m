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
% load 20240628-B-M-3.mat

%% simulation parameters
IampsModel = [-1.0, +1.5];  % [nA] step input amplitudes 

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
% IinjReal = zeros(1,Nall);           % input vector

for i = 1:length(IampsModel)         
    istart = Ninit+(Ninp+Nisi)*(i-1); % stating index for each input

    IinjModel(istart+(1:Ninp)) = 1000*IampsModel(i); % [nA]->[pA]
    % IinjReal(istart+(1:Ninp)) = 1000*IampsReal(i); % [nA]->[pA]

end


%%
% Setting parameters of Double Rz cell model
% Ion Channel Conductance Parameters
gKFactor = 58;     % default: 58
gAFactor = 162;    % default: 162
gNFactor = 645;    % default: 645
gLFactor = 10;     % default: 10
gHFactor = 1;      % default: 1

% Compartment Conductance Parameters
gCompFactor = 35;  % default: 35

% Reversal Potential Parameters
EL = -50;          % default:-50 leak
EH = -20;          % default:-20 Hyperpolarization
EN = +45;          % default:+45 Sodium
EK = -55;          % default:-55 Potassium 

% Capacitance Parameters
c1_factor = 130; % Soma Capacitance
c2_factor = 20;  % Gap Junction Capacitance
c3_factor = 30;  % Spike Initation Zone Capacitance

factorTauH = 10.0; 
factorTauN = 40.0; 
factorTauB = 30.0; 
factorTauZ = 100.0; 

% calling Rz cell model
gGCoupFactor_vals = [0 3 4 7 10 13 16 19 22];

% Run the Model and save responses
for Gi = 1:length(gGCoupFactor_vals)
    gGCoupFactor = gGCoupFactor_vals(Gi);
    [v1_Pre_vals,v2_Pre_vals,v3_Pre_vals,v1_Post_vals,v2_Post_vals,v3_Post_vals] = DoubleRz_MultiComp(IinjModel,dt, gCompFactor, gGCoupFactor, gNFactor, gKFactor, gAFactor, gLFactor, gHFactor, EL, EH, EN, EK, c1_factor, c2_factor, c3_factor, factorTauH, factorTauN, factorTauB, factorTauZ); % call the model

    stim_onset = 10001;
    stim_offset = 20000;
    window_length = 5000;

    v1_Pre = v1_Pre_vals;
    v1_Post = v2_Post_vals;
    v2_Pre = v2_Pre_vals;
    v2_Post = v2_Post_vals;
    v3_Pre = v3_Pre_vals;
    v3_Post = v3_Post_vals;

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
    c1 = [0,0,0.6];
    c2 = [0.6,0,0];

    % open fig 
    figure; clf; set(gcf,'Position',[50,50,1000,750]);

    % Presynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,1); cla; hold on; 
    % plot(t,v1_Pre,'-',color=c1); text(0.5,25,  sprintf("Spike count: %.0f", spike_count_Pre ), color=c1,fontsize=14);
    plot(t,v1_Pre,'-',color=c1); 

    txt_pre = "Spike count: " + num2str(spike_count_Pre);
    text(0.5,40,  txt_pre, color=c1, fontsize=10);

    txt_lat_pre = "Spike Latency: " + num2str(spike_latency_Pre);
    text(0.5,31,  txt_lat_pre, color=c1, fontsize=10);

    txt_hyper_pre = "Hyper. Amp.: " + num2str(amplitude_Pre);
    text(0.5,22,  txt_hyper_pre, color=c1, fontsize=10);

    txt_isi_pre = "ISI: " + num2str(inter_spike_interval_Pre);
    text(0.5,13,  txt_isi_pre, color=c1, fontsize=10);

    txt_gA = "gGC: " + num2str(gGCoupFactor);
    txt_gK = " EL: " + num2str(EL);

    text(0.5,4,  txt_gA+txt_gK, color=c1, fontsize=10);
    title('Membrane Potential of Presynaptic cell Soma')

    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Presynaptic compartment 2 (Gap)
    subplot(3,2,3); cla; hold on; 
    title('Membrane Potential of Pressynaptic cell Gap')
    plot(t,v2_Pre,'-',color=c1); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Presynaptic compartment 3 (spike-initiation zone)
    subplot(3,2,5); cla; hold on; 
    plot(t,v3_Pre,'-',color=c1);
    title('Membrane Potential of Presynaptic cell SIZ')
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Postsynaptic compartment 1 (soma = non-spiking)
    subplot(3,2,2); cla; hold on; 
    plot(t,v1_Post,'-',color=c2); 
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
    plot(t,v2_Post,'-',color=c2); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    % Postsynaptic compartment 2 (Gap)
    subplot(3,2,6); cla; hold on; 
    title('Membrane Potential of Postsynaptic cell SIZ')
    plot(t,v3_Post,'-',color=c2); 
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [s]'); ylabel('potential [mV]'); 

    %
    % ParameterFeatureData.Parameters.factorTauH;
    % ParameterFeatureData.Parameters.factorTauN;
    % ParameterFeatureData.Parameters.factorTauB;
    % ParameterFeatureData.Parameters.factorTauZ;
    % 
    % ParameterFeatureData.Parameters.EL;
    % ParameterFeatureData.Parameters.EH;
    % ParameterFeatureData.Parameters.EN;
    % ParameterFeatureData.Parameters.EK;
    % 
    % ParameterFeatureData.Parameters.c1_factor;
    % ParameterFeatureData.Parameters.c2_factor;
    % ParameterFeatureData.Parameters.c3_factor;
    % 
    % ParameterFeatureData.Parameters.gKFactor;
    % ParameterFeatureData.Parameters.gAFactor;
    % ParameterFeatureData.Parameters.gNFactor;
    % ParameterFeatureData.Parameters.gLFactor;
    % ParameterFeatureData.Parameters.gKFactor;
    % 
    % ParameterFeatureData.Parameters.gCompFactor;
    % ParameterFeatureData.Parameters.gGCoupFactor;
    % 
    % 
    % % 
    % ParameterFeatureData.Features.spike_count_Pre
    % ParameterFeatureData.Features.spike_latency_Pre
    % ParameterFeatureData.Features.amplitude_Pre
    % ParameterFeatureData.Features.inter_spike_interval_Pre
    % 
    % ParameterFeatureData.Features.spike_count_Post
    % ParameterFeatureData.Features.spike_latency_Post
    % ParameterFeatureData.Features.amplitude_Post
    % ParameterFeatureData.Features.inter_spike_interval_Post
    % 

end

%% function [V1_Pre,V2_Pre,V3_Pre, V1_Post,V2_Post,V3_Post] = RzDouble_Multi(Iinj,dt)
function [V1_Pre,V2_Pre,V3_Pre, V1_Post,V2_Post,V3_Post] = DoubleRz_MultiComp(Iinj,dt, gCompFactor, gGCoupFactor, gNFactor, gKFactor, gAFactor, gLFactor, gHFactor, EL,EH, EN, EK, c1_factor, c2_factor, c3_factor, factorTauH, factorTauN, factorTauB, factorTauZ)

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
    Vinit1 = -55; % [mV]
    Vinit2 = -55; % [mV]
    Vinit3 = -55; % [mV]

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

        dz_dt_Pre = ( infZ(V1_Pre(j)) - z_Pre ) / tauZ(V1_Pre(j), factorTauZ);
        dh_dt_Pre = ( infH(V3_Pre(j)) - h_Pre ) / tauH(V2_Pre(j), factorTauH);
        dn_dt_Pre = ( infN(V3_Pre(j)) - n_Pre ) / tauN(V2_Pre(j), factorTauN);
        db_dt_Pre = ( infB(V3_Pre(j)) - b_Pre ) / tauB(V2_Pre(j), factorTauB);

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

        dz_dt_Post = ( infZ(V1_Post(j)) - z_Post ) / tauZ(V1_Post(j), factorTauZ);
        dh_dt_Post = ( infH(V2_Post(j)) - h_Post ) / tauH(V2_Post(j), factorTauH);
        dn_dt_Post = ( infN(V2_Post(j)) - n_Post ) / tauN(V2_Post(j), factorTauN);
        db_dt_Post = ( infB(V2_Post(j)) - b_Post ) / tauB(V2_Post(j), factorTauB);


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


function x = tauH(v, factorTauH) % time scale [ms] for Na inactivation Default: 4.0     Changed to: 4.0  
    x = factorTauH * ones(size(v)); 
    % x = 10.0 * ones(size(v)); 
end
function x = tauN(v, factorTauN) % time scale [ms] for K activation    Default: 20.0     Changed to: 19.55
x = factorTauN * ones(size(v)); 
% x = 50 * ones(size(v)); 
end

function x = tauB(v, factorTauB) % time scale [ms] for KA inactivation Default: 20.0    Changed to: 20.1
  x = factorTauB * ones(size(v)); 
  % x = 30 * ones(size(v)); 
end
function x = tauZ(v, factorTauZ) % time scale [ms] for Ih activation   Default: 70.0
  x = factorTauZ * ones(size(v)); 
    % x = 70 * ones(size(v)); 
end

%% technical notes (by GA, June 2024)
%
%  + The model contains the hyperpolarization-activated current (Ih) 
%    in the soma to replicate the so-called voltage sag, which is a slow 
%    increase of the membrane potentilal with negative current injection.
%    To simulate a neuron without a sag, you can set gH=0 and increase 
%    gL to match the measured input resistance. 
%    
%    Modified: Changed to 
%    
%
%  + The KA current (inactivating K channel) was introduced to replicate 
%    the long first spike timing latency of the Rz cell. By changing this
%    conductance gA, the first spike latency can be adjusted. 
%
%    Modified: 
% 
%
%  + To simplify the model structure and to reduce the number of 
%    parameters, the following assumptions have been made. Some of these 
%    assumptions might be considered unphysiological. 
%  -- Na activation is instantaneous. Namely, m = infM(V) for all time. 
%  -- The time constants for all activation/inactivation functions are 
%     voltage independent. (This is clearly unrealistic)
%  -- The voltage-gated activation of KA current is negected. Only its 
%     inactivation is considered. 
%  -- Assuming that the soma is a sphere with a diameter of D=60 um, its 
%     surface area corresponds to pi*D^2 = 11000 um^2 and its capacitance
%     will be 110 pF under the standard capacitance density of 1 uF/cm^2.
%     Thus the used value of 200 pF in the model can be considered as the 
%     total capacitance of the soma and neighboring membrane areas. 
%  -- The capacitance of the spike initiation zone (SIZ) was chosen 
%     totally arbitrarily. 
%  -- The connection parameter of the two compartments gC greatly affects  
%     the input rsistance and was (somewhat arbitrarily) chosen based 
%     on subthreshold response data.
%  -- There are no direct measurements of the reversal potenitals. The 
%     used values are based on an "educated guess". 
% 
%  + The subthreshold response of the model is relatively robust to
%    changes in the model parameters. Namely, a small change in the 
%    setting leads to a small change in the outcome. 
% 
%  + The suprathreshold spiking activity, however, is very sensitive to 
%    parameters including gNa/gK/gA and the voltage dependence of these 
%    currents. A small change in these parameters can lead to a large 
%    change in spike counts and timings. 
%
%  + Model fitting was done by hand to replicate some of known response 
%    properties of the Rz cell. The model is clearly far from perfect 
%    and its predictability and usefulness are not at all guaranteed. 
%
%  + The model produces apparently more spikes than the real Rz cell in 
%    response to large input currents. It is not clear (at least to me) if
%    this issue can be easily fixed only by tweaking the model parameters. 
%    A complete change of the model structure might be needed instead.
%
%
%   Notes for meeting with Go:
% - Changing parameters for decreasing the spike frequency to 10-15 spikes,
%   and probably parameter combination. Changing activation or 
%   inactivation parameters.
%   
%   
%
%
