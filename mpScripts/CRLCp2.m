% CRLCp1.m

% Modelling a resonance circuit
% Parallel LC voltage divider circuit

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 171220 


clear all 
close all
clc

% ========================================================================
%   INPUTS   default values [ ]
% =======================================================================
% series resistance Z1 [ 1e4 ohms]    
   RS = 1.00e4;
% OUTPUT (LOAD) resistance Z2 [1e6 ohms]
   ROUT = 1e6;
% inductance and inductor resistance Z4 [10.3e-3 H  0 ohms]
   L = 10.3e-3;
   RL = 20;
% capacitance Z2 [10.4e-9 F]
   C = 10.4e-9;
% input voltage emf [10 V]
   V_IN = 10;
% frequency range [1000 to 30e3 Hz   5000]   
   fMin = 1000; fMax = 50e3;
   N = 5000;

   
% Experimental Data  frequency [kHz}
fE = [5 10 12.5 15.5 16 16.5 17 18 20 25 30 35 40 45 15 15.25 15.75];
GE = [0.039 0.115 0.232 0.810 0.690 0.554 0.464 0.304 0.196 0.098 0.073 0.058 0.051 0.045 0.776 0.81 0.789];

   
% =======================================================================  
%   CALCULATIONS
% =======================================================================   
   f = linspace(fMin,fMax, N); 
   w = (2*pi).*f;
   
% impedances
   Z1 = RS;                     % series resistance
   Z2 = ROUT;                   % output or load resistance
   Z3 = -1i ./ (w .*C);         % capacitive impedance (reactance)
   Z4 = RL + 1i .* w .* L;      % inductive impedance (resistance + reactance)

   Z5 = 1./ (1./Z2 + 1./Z3 + 1./Z4);    % parallel combination
   Z6 = Z1 + Z5;                        % total circuit impedance

% currents [A] and voltages [V]
   I1 = V_IN ./ Z6;
   V1 = I1 .* Z1;
   V_OUT = V_IN - V1;
   
   I2 = V_OUT ./ Z2;
   I3 = V_OUT ./ Z3;
   I4 = V_OUT ./ Z4;
     
  
   
% phases
   phi_OUT = angle(V_OUT);
   phi_1   = angle(V1);
   
   theta_1 = angle(I1); 
   theta_2 = angle(I2);
   theta_3 = angle(I3);
   theta_4 = angle(I4);
   
% Resonance frequencies and Bandwidth calculations
   f0 = 1/(2*pi*sqrt(L*C));
   G_V = abs(V_OUT ./ V_IN);     % voltage gain   
   Vpeak = max(G_V);             % max voltage gain
   VG3dB = Vpeak/sqrt(2);        % 3 dB points
   k = find(G_V == Vpeak);       % index for peak voltage gain
   f_peak = f(k);                % frequency at peak
   kB = find(G_V > VG3dB);       % indices for 3dB peak
   k1 = min(kB); f1 = f(k1);
   k2 = max(kB); f2 = f(k2);
   df = f2-f1;                   % bandwidth 
   Q = f0 / df;                  % quality factor
   
   P_OUT = V_OUT .* I2;          % power delivered to load
   


  
% =======================================================================
%  OUTPUTS IN COMMAND WINDOW
% =======================================================================
  fprintf('theoretical resonance frequency f0 = %3.0f ohms \n',f0);
  fprintf('peak frequency f_peak = %3.0f ohms \n',f_peak);
  fprintf('half power frequencies f1 = %3.0f Hz  %3.0f  Hz \n',f1,f2);
  fprintf('bandwidth   df = %3.0f Hz \n',df);
  fprintf('quality factor  Q  =  %3.2f  \n',Q);
  


% =======================================================================
%   GRAPHICS
% =======================================================================
   fss = 12;
figure (1)   % 
   pos = [0.07 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = fE; yP = GE;        % experimental data
   plot(xP, yP,'ko','linewidth',2);   
   
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('voltage gain  G');
   %set(gca,'yLim',[0 2e4]);
   set(gca,'fontsize',fss)
   
   
   figure (2)   % 
   pos = [0.47 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = fE; yP = GE;        
   plot(xP, yP,'ko','linewidth',2);   
   hold on
   xP = f./1e3; yP = G_V;                     
   plot(xP, yP,'r','linewidth',2);   

   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('voltage gain  G');
    