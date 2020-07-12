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
% inductance Z1 [10e-3 H]
   L = 3.8e-3;
% capacitance Z2 [1.0e-8 F]
   C = 1.1e-8;  
% series resistance Z3 [ 1e2 ohms]    
   R = 1.1e3;
% OUTPUT (LOAD) resistance Z4 [1e6 ohms]
   ROUT = 1e6;


% input voltage emf [10 V]
   V_IN = 10;
% frequency range [2000 to 50e3 Hz   5000]   
   fMin =1e3; fMax = 50e3; N = 5000;
   

   
% Experimental Data  frequency [kHz}   current [mA]
fE = [10 15 20 25 30 35 40 50];     
IE = [5 7 8.5 9 9 8.2 7.5 5.8]; 

   
% =======================================================================  
%   CALCULATIONS
% =======================================================================   
   f = linspace(fMin,fMax, N); 
   w = (2*pi).*f;
   
% impedances
   Z1 = 1i .* w .* L;            % inductive impedance (reactance)
   Z2 = -1i ./ (w .*C);    % capacitive impedance (reactance)
   Z3 = R;                 % series resistance
   Z4 = ROUT;              % output or load resistance
  
   Z5 = 1./ (1./Z3 + 1./Z4);    % parallel combination
   Z6 = Z1 + Z2 + Z5;                % total circuit impedance

% currents [A] and voltages [V]
   I1 = V_IN ./ Z6;
   I2 = I1;
   V1 = I1 .* Z1;
   V2 = I2 .* Z2;
   V_OUT = V_IN - V1 - V2;
   V3 = V_OUT; V4 = V_OUT;
   I3 = V_OUT ./ Z3;
   I4 = V_OUT ./ Z4;
   
% phases
     phi_1   = angle(V1);
     phi_2   = angle(V2);
     phi_3   = angle(V3);
     
     theta_1 = angle(I1); 
     theta_2 = angle(I2);
     theta_3 = angle(I3);
     theta_4 = angle(I4);
    
% Resonance frequencies and Bandwidth calculations
   f0 = 1/(2*pi*sqrt(L*C));
   
   Ipeak = max(abs(I1));           % max input current
   k = find(abs(I1) == Ipeak);     % index for peak voltage gain
   f_peak = f(k);                  % frequency at peak 
 
   I3dB = Ipeak/sqrt(2);           % 3 dB points
   kB = find(abs(I1) > I3dB);       % indices for 3dB peak
   k1 = min(kB); f1 = f(k1);
   k2 = max(kB); f2 = f(k2);
   df = f2-f1;                   % bandwidth 
   Q = f0 / df;                  % quality factor
   
% P_OUT = V_OUT .* I2;          % power delivered to load
   P_OUT = V_OUT .* I4;  
   P_emf = V_IN .* I1;


  
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

   xP = fE; yP = IE;        % experimental data
   plot(xP, yP,'ko','linewidth',2);   
   
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('current  I_{IN}  [ mA]');
   set(gca,'yLim',[0 10]);
   set(gca,'fontsize',fss)
   
   
   figure (2)   % 
   pos = [0.47 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = fE; yP = IE;        
   plot(xP, yP,'ko','linewidth',2);   
   hold on
   xP = f./1e3; yP = 1e3.*abs(I1);                     
   plot(xP, yP,'r','linewidth',2);   
   ylabel('current  I_{IN}  [ mA ]');
   grid on
   xlabel('frequency  f  [ kHz ]')
   set(gca,'fontsize',fss)
    