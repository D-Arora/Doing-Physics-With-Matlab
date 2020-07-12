% neurons_002c.m
% Finite difference method
% Leaky linear integrate and fire method
% Sub-threshold regime
% SI units unless specified

clear all
close all
clc

% INPUTS -----------------------------------------------------------------
   R = 1e4;                % membrane resistance  [ohms]
   C = 1e-8;               % membrane capacitance [F]
   I0 = 2e-3;              % external current amplitude  [mA]
   N = 5000;                % number of data points for calculations

   vTH = 15;               % threshold voltage
   vReset = -1;            % reset potential
   
% SETUP ------------------------------------------------------------------
   tau = R * C;             % time constant  [s]
   tmin = 0;                % time interval for modelling
   tmax = 10*tau;
   t = linspace(tmin, tmax, N);     % time steps [s]
   dt = t(2) - t(1);

   Iext = zeros(1,N);        % initialize membrane current  [mA]
   V = zeros(1,N);           % initialize membrane voltage  [mV]
   K = dt / tau;             % constant

   
% EXTERNAL STIMULUS  -----------------------------------------------------
  flagS = 2;        % enter value for type of external input

switch flagS
    
    case 1  % free solutions
        Iext = zeros(1,N);
        V(1) = -10;         % initial value of membrane potential [mV]
        
    case 2   % series of pulses - subthreshold repsonse
        num_start = 0.05 * N;
        num_end   = 0.05 * N;
        num_d     = 0.1 * N;
        num_width = round(0.05 * N);

       cn_max = round((N - num_end - num_start)/num_d);
       num_index = zeros(cn_max,1);

       for cn = 1 : cn_max
          num_index(cn) = round(num_start + cn * num_d);
          Iext(num_index(cn) : num_index(cn) + num_width) = I0;
       end  % if
       
    case 3   % step input current  
       num_start = 0.05 * N;
       Iext(num_start : end) = I0;
       
    
    case 4    % pulse input   
       num_start = 0.05 * N;
       num_end   = 0.2  * N;  
       num_width = round(num_end - num_start);
       Iext(num_start : 1*num_start + 1*num_width) = I0/1;
       
       
       
end   % switch
       
       
Vext = R .* Iext;           % external voltage stimulus


% Finite Difference method  ==============================================
for c = 1 : N - 1
   V(c+1) = V(c) - K * (V(c) - R * Iext(c));
   if V(c+1) > vTH; V(c) = vTH; V(c+1) = vReset; end;
       
end

% firing rates 

% GRAPHICS  ==============================================================
figure(1)
fs = 12;
set(gcf,'units','normalized');
set(gcf,'position',[0.02 0.40 0.2 0.25]);

subplot(2,1,1)
sfx = 1e3; sfy = 1; x = sfx .* t;
y = sfy .* Iext;
plot(x,y,'r','linewidth',2);
ylabel('Iext  [mA]','fontsize',fs);
axis([0 max(x) 0 2.2e-3])
grid on
set(gca,'fontsize',fs);


subplot(2,1,2)
sfx = 1e3; sfy = 1; x = sfx .* t;   y = sfy .* V;
plot(x,y,'linewidth',2);
xlabel('time  t  [ms]','fontsize',fs);
ylabel('V  [mV]','fontsize',fs);
grid on
set(gca,'fontsize',fs);
axis([0 max(x) -2 1.2*max(y)])


