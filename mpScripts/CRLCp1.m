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
   RS = 1e-1;
% OUTPUT (LOAD) resistance Z2 [1e6 ohms]
   ROUT = 1e6;
% inductance and inductor resistance Z4 [10.3e-3 H  0 ohms]
   L = 100e-3;
   RL = 100;
% capacitance Z2 [10.4e-9 F]
   C = 10e-6;
% input voltage emf [10 V]
   V_IN = 1;
% frequency range [1000 to 30e3 Hz   5000]   
   fMin = 100; fMax = 200;
   N = 5000;
   
   
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
     
   I_sum = abs(I1 - I2 - I3 - I4);
   
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
% TIME DOMAIN    select source frequency 
% =======================================================================
     fs = f_peak; kk = k;
   % fs = f1; kk = k1;
   % fs = f2; kk = k2;
   
   Ns = 500;
   ws = 2*pi*fs;
   Ts = 1/fs;
   tMin = 0;
   tMax = 3*Ts;
   t = linspace(tMin,tMax,Ns);
   
   emf = real(V_IN .* exp(1j*ws*t));
   v_OUT = real(abs(V_OUT(kk)) .* exp(1j*(ws*t + phi_OUT(kk))));
   v1 = real(abs(V1(kk)) .* exp(1j*(ws*t + phi_1(kk))));
   
   i1 = real(abs(I1(kk)) .* exp(1j*(ws*t + theta_1(kk))));
   i2 = real(abs(I2(kk)) .* exp(1j*(ws*t + theta_2(kk))));
   i3 = real(abs(I3(kk)) .* exp(1j*(ws*t + theta_3(kk))));
   i4 = real(abs(I4(kk)) .* exp(1j*(ws*t + theta_4(kk))));
   
   
% =======================================================================
%  OUTPUTS IN COMMAND WINDOW
% =======================================================================
  fprintf('theoretical resonance frequency f0 = %3.0f Hz \n',f0);
  fprintf('peak frequency f_peak = %3.0f Hz \n',f_peak);
  fprintf('half power frequencies f1 = %3.0f Hz  %3.0f  Hz \n',f1,f2);
  fprintf('bandwidth   df = %3.0f Hz \n',df);
  fprintf('quality factor  Q  =  %3.2f  \n',Q);
  fprintf('current at junction I_sum  =  %3.2f mA  \n',max(1e3*I_sum));


% =======================================================================
%   GRAPHICS
% =======================================================================
   fss = 12;
figure(1)   % impedances -----------------------------------------------
   pos = [0.07 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f./1e3; yP = abs(Z3);        % capacitve reactance
   plot(xP, yP,'k','linewidth',2);   
   hold on
   yP = abs(Z4);                     % inductive impedance
   plot(xP, yP,'m','linewidth',2);   
   yP = abs(Z5);                     % output impedance
   plot(xP, yP,'b','linewidth',2);
   legend(' Z_C',' Z_L',' Z_{OUT}');
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('impedances  | Z |  [ \Omega ]');
  % set(gca,'yLim',[0 2e4]);
   set(gca,'fontsize',fss)
   
figure(2)    % output voltage ------------------------------------------
   pos = [0.05 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   yyaxis left
   xP = f./1e3; yP = abs(G_V);
   plot(xP,yP,'b','lineWidth',2);
   %set(gca,'yLim',[0 4e4]);
    xlabel('frequency  f  [ kHz ]')
   ylabel('Gain  G_V =  |V_{OUT} /  V_{IN}|');
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',fss);
   grid on
   
   yyaxis right
   yP = phi_OUT./pi;
   plot(xP,yP,'r','lineWidth',2);
   %set(gca,'yLim',[0 4e4]);
   ax.YAxis(2).Color = 'r';
   set(gca,'yTick',-0.5:0.25:0.5);
   set(gca,'fontsize',fss);
   ylabel('phase   \phi_{OUT}    [ rad / \pi ]');
   
figure(3)   % current source -------------------------------------------
   pos = [0.38 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   yyaxis left
   xP = f./1e3; yP = 1e3.*abs(I1);
   plot(xP,yP,'b','lineWidth',2);
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',14);
   ylabel('source current  I_S  [mA]');
   
   xlabel('frequency  f  [ kHz ]')
   grid on
   
   yyaxis right
   yP = theta_1./pi;
   plot(xP,yP,'r','lineWidth',2);
   %set(gca,'yLim',[0 4e4]);
   ax.YAxis(2).Color = 'r';
   set(gca,'yTick',-0.5:0.25:0.5);
   set(gca,'fontsize',fss);
   ylabel('phase   \theta_{S}    [ rad / \pi ]');
  
figure(4)   % current phases -------------------------------------------
   pos = [0.38 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f./1e3; yP = theta_3 ./ pi;        
   plot(xP, yP,'k','linewidth',2);   
   hold on
   yP = theta_4 ./ pi;
   plot(xP, yP,'m','linewidth',2);  
   legend(' I_C',' I_L');
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('current phase    [ rad / \pi ]');
   set(gca,'fontsize',fss)   

figure(5)   % current parallel combination ------------------------------
   pos = [0.69 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f./1e3; yP = 1e3.*abs(I3);        
   plot(xP, yP,'k','linewidth',2);   
   hold on
   yP = 1e3.*abs(I4);
   plot(xP, yP,'m','linewidth',2);  
   legend(' I_C',' I_L');
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('current magnitude [ mA ]');
   set(gca,'fontsize',fss)
   
figure(6)   % power -----------------------------------------------------
   pos = [0.69 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f./1e3; yP = 1e3.*abs(P_OUT);        
   plot(xP, yP,'b','linewidth',2);   
      grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('P_{OUT}   [ mW ]');
   set(gca,'fontsize',fss)   
   
figure (7)   % voltages: time domain -------------------------------------
   pos = [0.2 0.5 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = t; yP = emf;        
   plot(xP, yP,'b','linewidth',2);   
   hold on
   yP = v_OUT;
   plot(xP, yP,'r','linewidth',2); 
   yP = v1;
   plot(xP, yP,'k','linewidth',2); 
   set(gca,'xLIm',[0 tMax]);
   
   legend(' emf',' v_{OUT}', ' v_1');
   grid on
   xlabel('time [ s ]')
   ylabel('voltages  V  [ V ]');
   
   tm1 = 'f  =  ';
   tm2 = num2str(fs, '%3.0f \n');
   tm3 = '  Hz';
   tm = [tm1 tm2 tm3];
   h_text = text(1e-5,8,tm);
   set(h_text,'fontsize',fss);
   
   set(gca,'fontsize',fss)      

   
   figure (8)   % currrents: time domain ---------------------------------
   pos = [0.6 0.5 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   i2Max = max(i2).*1e3;
   mt1 = 'i_{OUTmax}  =  ';
   mt2 = num2str(i2Max,'%3.3f \n');
   mt3 = '  mA';
   mt = [mt1 mt2 mt3];
   
   xP = t; yP = 1e3.*i1;        
   plot(xP, yP,'b','linewidth',2);   
   hold on
   yP = 1e3.*i2;
   plot(xP, yP,'r','linewidth',2); 
   yP = 1e3.*i3;
   plot(xP, yP,'k','linewidth',2);
   yP = 1e3.*i4;
   plot(xP, yP,'m','linewidth',2);
   set(gca,'xLIm',[0 tMax]);
   legend(' i_S',' i_{OUT}', ' i_C',' i_L');
   grid on
   xlabel('time [ s ]')
   ylabel('currents   I  [ mA ]');
   title(mt,'fontweight','normal');
   
   tm1 = 'f  =  ';
   tm2 = num2str(fs, '%3.0f \n');
   tm3 = '  Hz';
   tm = [tm1 tm2 tm3];
   h_text = text(1e-5,8,tm);
   set(h_text,'fontsize',fss);
   set(gca,'yLim',[-10 10]);
   set(gca,'fontsize',fss)  