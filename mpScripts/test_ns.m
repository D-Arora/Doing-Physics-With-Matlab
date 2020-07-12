% LIFspike_ELIF.m
%
% The code simulates the leaky integrate-and-fire (LIF) neuron then the
% exponential leaky integrate-and-fire (ELIF) neuron, to compare traces of
% the membrane potential produced by the two models. A spike is added
% artificially in the LIF model.
%
% This code is used to produce Figure 2.6 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the default plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
    figure(1)
    clf;

%% Simulation parameters
dt = 0.0001;                % time-step
t = 0:dt:0.5;               % vector of time-points
ton = 0.15;                 % time to begin applied current (onset)
toff = 0.35;                % time to end applied current (offset)
non = round(ton/dt);        % time-point index of current onset
noff = round(toff/dt);      % time-point index of current offset

%% Parameters for the model neurons
tau = 0.010;                % membrane time constant
E_L = -0.070;               % leak potential (also resting potential)
Vth = -0.050;               % threshold potential (to produce spike)
Vreset = -0.080;            % reset potential (post-spike)
Cm = 100e-12;               % total membrane capacitance
G_L = Cm/tau;               % total membrane conductance (leak conductance)

% The following are specific parameters for the ELIF model neuron
Vmax = 0.050;               % Clip spikes at this value of membrane potential
Delta_th = 0.005;           % Range of membrane potential for spike acceleration

Iapp = [210e-12 210e-12];   % Same applied current for the 2 models
Ntrials = length(Iapp);     % One trial for each model

for trial = 1:Ntrials;              % loop through different trials
    I = zeros(size(t));             % vector for current at each time-point
    I(non:noff) = Iapp(trial);      % add the applied current for the trial
    V = E_L*ones(size(t));          % initialize the membrane potential vector
    spikes = zeros(size(t));        % initialize a vector to record spikes
    
    
    if ( trial == 1 )               % Standard LIF model (see LIF_model.m)
        for i = 2:length(t);            % loop through all time points
            % next line: Forward Euler method to update membrane potential
            % see Eq. 1.9 in the book
            V(i) = V(i-1) + dt*(I(i) +G_L*(E_L-V(i-1)))/Cm;
            if ( V(i) > Vth)            % if potential is above threshold
                spikes(i) = 1;          % record the spike at that time-point
                V(i) = Vreset;          % reset the potential
                V(i-1) = Vmax;          % add a spike on prior time-point for visual purposes
            end;
        end;                            % end the loop & go to next time-point
    else;                               % Otherwise simulate the ELIF model  
        for i = 2:length(t)             % loop through all time points
            % next line: Forward Euler update of membrane potential 
            V(i) = V(i-1) + dt*(I(i) +G_L*(E_L-V(i-1) + ...
                Delta_th*exp((V(i-1)-Vth)/Delta_th) ))/Cm;
            if ( V(i) > Vmax )          % if potential is greater than Vmax
                spikes(i) = 1;          % record the spiketime
                V(i) = Vreset;          % Reset the membrane potential
                V(i-1) = Vmax;          % add a spike on prior time-point for visual purposes
            end
        end;                            % Loop to next time-point
        
    end
    
    % Set up the sub-panel for the figure
    subplot('position',[0.46*trial-0.3 0.23 0.34 0.64])
    plot(t,1000*V,'k');          % plot membrane potential versus time
    if ( trial == 1 )       % if LIF model used
        ylabel('V_m (mV)')
        title('LIF + spike')
    else;                   % if ELIF model used
        title('ELIF')
    end
    axis([0.2 0.3 1000*(Vreset-0.005) 1000*(Vmax+0.005)])
    xlabel('Time (sec)')
    
end

annotation('textbox',[0.00 0.97 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.53 0.97 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')