% CRLCs1.m

% Modelling a resonance circuit
% RLC series resoance circuit

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% 180103 


clear all
close all
clc

% ========================================================================
%   INPUTS   default values [ ]
% =======================================================================

% inductance Z1 [10e-3 H]
   L = 100e-3;
% capacitance Z2 [1.0e-8 F]
   C = 10e-6;  
% series resistance Z3 [ 1e2 ohms]    
   R = 1000;
% OUTPUT (LOAD) resistance Z4 [1e6 ohms]
   ROUT = 1e6;
% input voltage emf [10 V]
   V_IN = 110;
% frequency range [2000 to 50e3 Hz   5000]   
   fMin = 60; fMax = 400; N = 5000;
   
   
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
% TIME DOMAIN    select source frequency 
% =======================================================================
   c = 1;
%  c = 1 fs = f_peak; 
%  c = 2 fs = f1;
%  c = 3 fs = f2
   if c == 1; kk = k; fs = f_peak; kk = k; end
   if c == 2; kk = k1; fs = f1; end
   if c == 3; kk = k2; fs = f2; end
fs = 60;
    Ns = 500;
    ws = 2*pi*fs;
    Ts = 1/fs;
    tMin = 0;
    tMax = 3*Ts;
    t = linspace(tMin,tMax,Ns);
    
    emf = real(V_IN .* exp(1j*ws*t));
    v1 = real(abs(V1(kk)) .* exp(1j*(ws*t + phi_1(kk))));
    v2 = real(abs(V2(kk)) .* exp(1j*(ws*t + phi_2(kk))));
    v3 = real(abs(V3(kk)) .* exp(1j*(ws*t + phi_3(kk))));
    
    i1 = real(abs(I1(kk)) .* exp(1j*(ws*t + theta_1(kk))));
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
%   fprintf('current at junction I_sum  =  %3.2f mA  \n',max(1e3*I_sum));


% =======================================================================
%   GRAPHICS
% =======================================================================
   fss = 12;
figure(1)   % impedances -----------------------------------------------
   pos = [0.07 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f./1e3; yP = abs(Z2)./1e3;        % capacitve reactance
   plot(xP, yP,'k','linewidth',2);   
   hold on
   yP = abs(Z1)./1e3;                     % inductive impedance
   plot(xP, yP,'m','linewidth',2);   
   yP = abs(Z6)./1e3;                     % total impedance
   plot(xP, yP,'b','linewidth',2);
   legend(' Z_C',' Z_L',' Z_{TOTAL}');
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('impedances  | Z |  [ k\Omega ]');
   %set(gca,'yLim',[0 5]);
   set(gca,'fontsize',fss)
   
figure(2)    % current   ------------------------------------------
   pos = [0.05 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   yyaxis left
   xP = f./1e3; yP = abs(I1.*1e3);
   plot(xP,yP,'b','lineWidth',2);
   %set(gca,'yLim',[0 100]);
   xlabel('frequency  f  [ kHz ]')
   ylabel('input current  I_{IN}  [ mA ]');
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',fss);
   grid on
   
   yyaxis right
   yP = theta_1./pi;
   plot(xP,yP,'r','lineWidth',2);
   %set(gca,'yLim',[0 4e4]);
   ylabel('phase   I_{IN}    [ rad / \pi ]');
   ax.YAxis(2).Color = 'r';
   set(gca,'yLim',[-0.52 0.52]);
   set(gca,'yTick',-0.5:0.25:0.5);
   set(gca,'fontsize',fss);
   
   
figure(3)   % voltages -------------------------------------------
   pos = [0.38 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   
   xP = f./1e3; yP = V_IN .* ones(N,1);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   yP = abs(V1);
   plot(xP,yP,'m','lineWidth',2);
   yP = abs(V2);
   plot(xP,yP,'k','lineWidth',2);
   yP = abs(V3);
   plot(xP,yP,'r','lineWidth',2);
   yP = abs(V1+V2);
   plot(xP,yP,'g-','lineWidth',2);
   
   set(gca,'fontsize',14);
   %set(gca,'yLim',[-0.20 V_IN + 0.2]);
   ylabel('voltages  [ V ]');
   xlabel('frequency  f  [ kHz ]')
   legend('emf','V_L','V_C','V_{OUT}','V_L + V_C');
   grid on
   %set(gca,'yLim',[0 100]);

%   
figure(4)   % power -------------------------------------------
   pos = [0.38 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   yyaxis left
   xP = f./1e3; yP = 1e3.*abs(P_OUT);        
   plot(xP, yP,'b','linewidth',2);   
   xlabel('frequency  f  [ kHz ]')
   ylabel('P_{OUT}   [ mW ]');
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',fss);
   grid on
   
   
   yyaxis right
   yP = abs(P_emf);        
   plot(xP, yP,'r','linewidth',2);  
   grid on
   ylabel('P_{IN}   [ W ]');
   set(gca,'fontsize',fss) 
   set(gca,'yLim',[0 1]);
   

figure(5)   % voltage for L and C ------------------------------
   pos = [0.69 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f./1e3; yP = phi_1 ./ pi;        
   plot(xP, yP,'m','linewidth',2);   
   hold on
   yP = phi_2 ./ pi;
   plot(xP, yP,'k','linewidth',2);  
   legend(' V_L',' V_C');
   grid on
   xlabel('frequency  f  [ kHz ]')
   ylabel('voltage phase    [ rad / \pi ]');
   set(gca,'fontsize',fss)   
   
figure(6)   % time domain: voltages  -------------------------------------
   pos = [0.69 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = t; yP = emf;        
   plot(xP, yP,'b','linewidth',3);   
   hold on
   yP = v1;
   plot(xP, yP,'m','linewidth',2); 
   yP = v2;
   plot(xP, yP,'k','linewidth',2); 
   set(gca,'xLIm',[0 tMax]);
   yP = v3;
   plot(xP, yP,'r','linewidth',1); 
   set(gca,'xLIm',[0 tMax]);
      
   legend(' emf',' v_L', ' v_C',' v_{OUT}');
   grid on
   xlabel('time t  [ s ]')
   ylabel('voltages  V  [ V ]');
   
   tm1 = 'f  =  ';
   tm2 = num2str(fs, '%3.0f \n');
   tm3 = '  Hz';
   tm = [tm1 tm2 tm3];
   h_tm = title(tm,'fontWeight','normal');
   set(h_tm,'fontsize',fss);
   
   set(gca,'fontsize',fss)      

   
   
figure(7)   % currrents: time domain ---------------------------------
   pos = [0.6 0.5 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = t; yP = 1e3.*i1;        
   plot(xP, yP,'m','linewidth',3);   
   hold on
   yP = 1e3.*i3;
   plot(xP, yP,'b','linewidth',2);
   yP = 1e3.*i4;
   plot(xP, yP,'r','linewidth',2);
   yP = max(1e3.*i1) .* emf./V_IN;
   plot(xP, yP,'g','linewidth',1)
   set(gca,'xLIm',[0 tMax]);
   legend(' i_{IN}',' i_{R}', ' i_{Load}', 'emf [a.u.]');
   grid on
   xlabel('time [ s ]')
   ylabel('currents   I  [ mA ]');
   
   tm1 = 'f  =  ';
   tm2 = num2str(fs, '%3.0f \n');
   tm3 = '  Hz';
   tm = [tm1 tm2 tm3];
   h_tm = title(tm,'fontWeight','normal');
   set(h_tm,'fontsize',fss);
   %set(gca,'yLim',[-100 100])
   set(gca,'fontsize',fss)  
  
  
   
   