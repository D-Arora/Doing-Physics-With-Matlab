% ns_LIF_002.m
% Finite difference method
% Leaky linear integrate and fire method
% Sub-threshold regime
% SI units unless specified

clear
close all
clc

% INPUTS -----------------------------------------------------------------
   R = 10e6;                % membrane resistance  [10e6 ohms]
   tau = 10;                % time constant  [20 s] 
   I0 = 1.2e-6;               % external current amplitude  [mA]
   N = 5000;                % number of data points for calculations [5000]

   vTH = -50;               % threshold voltage
   vReset = -80;            % reset potential
   vRest = -60;             % rest potential
   vSpike = 20;             % spike potential
   
   tARP = 5;                % absolute refractrory period
   
% SETUP ------------------------------------------------------------------
   
   tmin = 0;                % time interval for modelling
   tmax = 100; %50*tau;
   t = linspace(tmin, tmax, N);     % time steps [s]
   dt = t(2) - t(1);

   Iext = zeros(1,N);        % initialize membrane current  [mA]
   V = vRest .* ones(1,N);           % vM  initialize membrane voltage  [mV]
   K = dt / tau;             % constant

   nARP = round(tARP / dt);
   
% EXTERNAL STIMULUS  -----------------------------------------------------
 
flagS = 2;        % enter value for type of external input

switch flagS
    
    case 1  % free solutions
        Iext = zeros(1,N);
        V(1) = -80;
        
    case 2   % series of pulses - subthreshold repsonse
        num_start = 0.05 * N;
        num_end   = 0.05 * N;
        num_d     = 0.05 * N;
        num_width = round(0.025 * N);

       cn_max = round((N - num_end - num_start)/num_d);
       num_index = zeros(cn_max,1);

       for cn = 1 : cn_max
          num_index(cn) = round(num_start + cn * num_d);
          Iext(num_index(cn) : num_index(cn) + num_width) = I0;
       end  % if
       
    case 3   % step input current  
       num_start = 0.1 * N;
       Iext(num_start : end) = I0;
       
    
    case 4    % pulse input   
       num_start = 0.1 * N;
       num_end   = 0.3  * N;  
       num_width = round(num_end - num_start);
       Iext(num_start : num_start + num_width) = I0;
       %Iext(num_start : 0.5*num_start + 0.5*num_width) = 2*I0;
       
    case 5    % ramp input 
       Iext = 5e-8 .* t;
    
    case 6    % nosiy input   
       Iext = I0 .* rand(1,N);
       
end   % switch
   
       
       Vext = R .* Iext;           % external voltage stimulus


% Finite Difference method  ==============================================
d = -1;
for c = 1 : N - 1
   V(c+1) = V(c) + K * (-V(c) + vRest + Vext(c));
  
   if V(c+1) > vTH
      V(c) = vSpike;
      V(c+1) = vReset;
      d = nARP;
   end 
   
   if d > 0
      V(c+1) = vReset;
      d = d - 1;
   end
   
end


% Firing rates: flagF = 1 (yes) flagF = 0 (no) ===========================
   flagF = 0;
   if flagF == 1
      indexFire = find(V == vSpike);          % indeices for spikes
      indexFire2 = indexFire(2:end);
      indexFire1 = indexFire(1:end-1);
      ISI = dt .* (indexFire2 - indexFire1);  % interspike interval
      fireRate = 1000 ./ISI;                  % firning rate
      f = mean(fireRate);                     % frequency: mean firing rate
      
   disp('Interspike times   ISI   [ms]  ');
   fprintf(' %2.2f ',ISI);
   disp('   ');
   disp('Neuron firing rate  [Hz]   ');
   fprintf(' %2.2f ',fireRate);
   disp('  ');
   fprintf('mean firing rate  f =  %2.2f \n',f);
   disp('  ');
   
   if flagS == 5     % ramp input stimuls
      figure(2)
      fs = 12;
      set(gcf,'units','normalized');
      set(gcf,'position',[0.42 0.40 0.2 0.2]);
      xP = 1e6 .* Iext(indexFire2);  yP = fireRate;
      plot(xP,yP,'b','linewidth',2);
      ylabel('f  [Hz]','fontsize',fs);
      xlabel('I_{ext}  [nA]','fontsize',fs);
      axis([0 max(xP) 0 1.2*max(yP)])
      grid on
      set(gca,'fontsize',fs);
   end  % if
   
   end  % flagF

% GRAPHICS  ==============================================================
figure(1)
fs = 12;
set(gcf,'units','normalized');
%set(gcf,'position',[0.02 0.40 0.4 0.4]);
set(gcf,'position',[0.02 0.40 0.25 0.35]);

hsp = subplot(2,1,1);
set(hsp,'Position', [0.1600 0.65 0.77 0.25]);
sfx = 1; sfy = 1e6;
xP = sfx .* t;
yP = sfy .* Iext;
plot(xP,yP,'r','linewidth',2);
ylabel('I_{ext}  [ nA ]','fontsize',fs);
%set(gca,'yLim',[0 1.1]);
%axis([0 max(xP) 0 1.2*max(yP)])
grid on
set(gca,'fontsize',fs);
box on

hsp = subplot(2,1,2);
set(hsp,'Position', [0.1600 0.15 0.77 0.40]);
sfx = 1; sfy = 1;
xP = sfx .* t;   yP = sfy .* V;
plot(xP,yP,'linewidth',2);
xlabel('time  t  [ ms ]','fontsize',fs);
ylabel('v_M  [ mV ]','fontsize',fs);
grid on
set(gca,'fontsize',fs);
%set(gca,'yLim',[-60 -55]);
%axis([0 max(xP) 1.1*min(yP) 1.4*max(yP)])
box on

