function [v1,v2, t, Iinj, z, m, h, n, b] = RzCell_prelim

    %% Two-compartment Leech Retzius cell Model (preliminary version)
    % --- MSc Course Neuroscience Team Project in Summer Semester 2024 
    % --- written by Go Ashida
    %
    % - Output: 
    %  v1  simulated membrane potential of the soma [mV]
    %  v2  simulated membrane potential of the spike initiation zone (SIZ) [mV]
    %  t   time vector [ms]
    % 
    
    %% simulation parameters 
    Iamps = [-2.0,-1.0,+0.5,+1.5]; % [nA] step input amplitudes 
    dt = 0.01; % [ms] time step 
    Tinit = 1000; % [ms] initial silent period
    Tisi  = 1000; % [ms] inter-stimulus-interval 
    Tinp  = 500; % [ms] duration of each step current
    Ninit = round(Tinit/dt);
    Nisi  = round(Tisi/dt); 
    Ninp  = round(Tinp/dt); 
    Nall = Ninit + (Ninp+Nisi)*length(Iamps); % total number of time steps
    Tall = Nall*dt/1000; % [sec] total simulation time
    t = (0:Nall-1)*dt/1000; % [sec] time vector
    
    %% make current input
    Iinj = zeros(1,Nall); % input vector
    for i = 1:length(Iamps)
     istart = Ninit+(Ninp+Nisi)*(i-1); % stating index for each input
     Iinj(istart+(1:Ninp)) = 1000*Iamps(i); % [nA]->[pA]
    end
    
    %% calling Rz cell model
    [v1,v2, z, m, h, n, b] = RzDouble(Iinj,dt); % call the model
    
    %% plotting
    
    % downsample data vectors to reduce plotting time
    s = t(1:10:end);
    u0 = Iinj(1:10:end)/1000;
    u1 = v1(1:10:end);
    u2 = v2(1:10:end);
    
    % color settings
    c0 = [0,0.6,0]; 
    c1 = [0,0,0.6]; 
    c2 = [0.6,0,0];
    
    % open fig 
    figure(123); clf; set(gcf,'Position',[50,50,1000,750]);
    
    % input current
    subplot(3,1,1); cla; hold on; 
    plot(s,u0,'-',color=c0); text(0.5,1.8,'input',color=c0,fontsize=14);
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-2.5,2.5]); set(gca,'ytick',-3:1:3);
    xlabel('time [sec]'); ylabel('current [nA]'); 
    
    % compartment 1 (soma = non-spiking)
    subplot(3,1,2); cla; hold on; 
    plot(s,u1,'-',color=c1); text(0.5,25,'soma',color=c1,fontsize=14);
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [sec]'); ylabel('potential [mV]'); 
    
    % compartment 2 (spike-initiation zone)
    subplot(3,1,3); cla; hold on; 
    plot(s,u2,'-',color=c2); text(0.5,25,'spike-initiation zone',color=c2,fontsize=14);
    xlim([0,Tall]); set(gca,'xtick',0:0.5:Tall);
    ylim([-105,45]); set(gca,'ytick',-100:20:40);
    xlabel('time [sec]'); ylabel('potential [mV]'); 
    


end % end of RzCell_prelim()


function [V1,V2, z, m, h, n, b] = RzDouble(Iinj,dt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Two-compartment leech Retzius cell model 
    % - Input:
    %  Iinj  vector for injected current input [pA]
    %  dt    time step [ms]
    % - Output:
    %  V1    membrane potential of the soma [mV]
    %  V2    membrane potential of the spike-initiation zone (SIZ) [mV]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% assigning membrane parameters
    
    % capacitances
    c1 = 200 * 1e-6; % [uF] membrane capacitance of soma: default 200 [pF]
    c2 = 200 * 1e-6; % [uF] membrane capacitance of SIZ: default 200 [pF]
    
    % conductances 
    gL =  8 * 1e-6; % [mS] leak conductance of soma: default 8 [nS]
    gH =  4 * 1e-6; % [mS] Ih conductance of soma: default 4 [nS]      % To simulate a neuron without a sag, you can set gH=0 and increase gL to match the measured input resistance. 
    gN = 1800 * 1e-6; % [mS] Na conductance of SIZ: default 1800 [nS]
    gK = 180 * 1e-6; % [mS] K conductance of SIZ: default 180 [nS]
    gA = 180 * 1e-6; % [mS] KA conductance of SIZ: default 180 [nS]
    gC =  80 * 1e-6; % [mS] conductance between the two compartments: default 80 [nS]
    
    % reversal potentials
    EL = -25; % [mV] leak reversal potential (soma only)
    EH = -25; % [mV] Ih reversal potential (soma only)
    EN = +60; % [mV] Na reversal potential (SIZ only)
    EK = -60; % [mV] K reversal potential (SIZ only)
    
    %% data vectors and initial values
    Ntotal = length(Iinj); % total number of data points
    V1 = zeros(1,Ntotal); % [mV] membrane potential of somatic compartment
    V2 = zeros(1,Ntotal); % [mV] membrane potential of spike-initiation zone
    
    % initial values
    Vinit1 = -53; % [mV]
    Vinit2 = -57; % [mV]
    V1(1) = Vinit1; 
    V2(1) = Vinit2; 
    z = infZ(Vinit1); % Ih activation variable (soma only)
    m = infM(Vinit2); % Na activation variable (SIZ only)
    h = infH(Vinit2); % Na inactivation variable (SIZ only)
    n = infN(Vinit2); % K activation variable (SIZ only)
    b = infB(Vinit2); % KA inactivation variable (SIZ only)
    
    %% calculate membrane response step-by-step 
    for j = 1:Ntotal-1
    
        % ionic currents (soma): g[mS] * V[mV] = I[uA]
        IL = gL     * (EL-V1(j));   % leak current
        IH = gH * z * (EH-V1(j));   % Ih current
        IC = gC     * (V2(j)-V1(j)); % current from compartment 2 to compartment 1
    
        % ionic currents (spike-initiation zone)
        IN = gN * m^4 * h * (EN-V2(j)); % Na current
        IK = gK * n^2     * (EK-V2(j)); % K current
        IA = gA * b       * (EK-V2(j)); % KA current
    
        % derivatives: I[uA] / C[uF] * dt[ms] = dv[mV]
        dv1_dt = ( IL + IH + IC + Iinj(j)*1e-6 ) / c1; % -IGap
        dv2_dt = ( IN + IK + IA - IC  ) / c2; 
        dz_dt = ( infZ(V1(j)) - z ) / tauZ(V1(j));
        dh_dt = ( infH(V2(j)) - h ) / tauH(V2(j));
        dn_dt = ( infN(V2(j)) - n ) / tauN(V2(j));
        db_dt = ( infB(V2(j)) - b ) / tauB(V2(j));
    
        % Compartment 1 - calculate next step 
        V1(j+1) = V1(j) + dv1_dt * dt; 

        % Compartment 2 - calculate next step 
        V2(j+1) = V2(j) + dv2_dt * dt; 

        m = infM(V2(j)); % Na activation is assumed to be instantaneous
        z = z + dz_dt * dt; 
        h = h + dh_dt * dt; 
        n = n + dn_dt * dt; 
        b = b + db_dt * dt; 
    
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
function x = tauH(v) % time scale [ms] for Na inactivation 
  x = 4.0 * ones(size(v)); 
end
function x = tauN(v) % time scale [ms] for K activation 
  x = 6.0 * ones(size(v)); 
end
function x = tauB(v) % time scale [ms] for KA inactivation 
  x = 20.0 * ones(size(v)); 
end
function x = tauZ(v) % time scale [ms] for Ih activation 
  x = 200 * ones(size(v)); 
end

%% technical notes (by GA, June 2024)
%
%  + The model contains the hyperpolarization-activated current (Ih) 
%    in the soma to replicate the so-called voltage sag, which is a slow 
%    increase of the membrane potentilal with negative current injection.
%    To simulate a neuron without a sag, you can set gH=0 and increase 
%    gL to match the measured input resistance. 
%
%  + The KA current (inactivating K channel) was introduced to replicate 
%    the long first spike timing latency of the Rz cell. By changing this
%    conductance gA, the first spike latency can be adjusted. 
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


