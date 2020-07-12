% CacFilters1.m

% Modelling FILTER circuits
% 

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180103 


clear all
close all
clc

% ========================================================================
%   INPUTS   default values [ ]
% =======================================================================

% series resistance Z1 [ 1e3 ohms] 
   R1 = 1e3;
% capacitance Z2 [1.0e-6 F]
   C2 = 1e-6;
% capacitance Z3 [1.0e-6 F]
   C3 = 1e-6; 
% resistance Z4 [ 1e3 ohms] 
   R4 = 1e3;   
% output resistance Z5 [ 1e3 ohms]    
   R5 = 1e3;  R_OUT = R5;
% input peak voltage for apllied emf [10 V]   
    V_IN = 10;
% frequency range {N = 1000 fMIn = 10 Hz fMAx = 5e3 Hz]    
    N = 5000;
    fMin = 10;
    fMax =10e3;
% input signal frequency fs for time domain graph [319.6897] Hz
    fs = 319.6897;

% =======================================================================  
%   CALCULATIONS
% =======================================================================   
   f = linspace(fMin,fMax, N); 
   w = (2*pi).*f;
   
% impedances
   Z1 = R1;
   Z2 = -1i ./ (w .*C2);    
   Z3 = -1i ./ (w .*C3);
   Z4 = R4;
   Z5 = R5;
   Z6 = 1./(1./Z4 + 1./Z5);
   Z7 = Z3 + Z6;
   Z8 = 1./(1./Z2 + 1./Z7);
   Z9 = Z1 + Z8;
   
% currents [A] and voltages [V] 
   i1 = V_IN ./ Z9; i_IN = i1;
   v1 = i1 .* Z1;
   v2 = V_IN - v1;
   i3 = v2./Z3;
   i7 = i1 - i3;
   v6 = i7 .* Z6;
   i5 = v6./Z5;
   i_OUT = i5;
   I_OUT = abs(i_OUT);
   v_OUT = v6;
   V_OUT = abs(v_OUT);
   P_OUT = V_OUT .* I_OUT;
   
  
%    v3 = v2; v_OUT = v2;
%    i2 = v2 ./ Z2; 
%    i3 = v3 ./ Z3; i_OUT = i3; 
%    
% % Magnitudes: currents [mA] and voltages [V]
%    V1 = abs(v1);
%    V2 = abs(v2);
%    V3 = V2;
%    V_OUT = V2;
%    I1 = 1e3.*abs(i1); I_IN = I1;
%    I2 = 1e3.*abs(i2);
%    I3 = 1e3.*abs(i3);
%    I_OUT = I3;
%    
% % Phases   [rad /pi]
%      phi_1   = angle(v1)./pi;
%      phi_2   = angle(v2)./pi;
%      phi_3   = angle(v3)./pi;
%      phi_OUT = phi_3;
%      
%      theta_1 = angle(i1)./pi; 
%      theta_2 = angle(i2)./pi;
%      theta_3 = angle(i3)./pi;
%  
% % Powers  [mW]
%      P1 = V1.*I1;
%      P2 = V2.*I2;
%      P3 = V3.*I3; P_OUT = P3;
%      P_IN = V_IN .* I_IN;
% 
% % Voltage gain Av and -3dB point
%     Av = V_OUT ./ V_IN;
%     AvdB = real(20*log10(Av));
%     
%     AvdB_max = max(AvdB);
%     AvdB_3 = AvdB_max-3;
%     K = find(AvdB < AvdB_3,1);
%     if flagF == 2
%       K = find(AvdB > AvdB_3,1);
%     end
%     fC = f(K);
%     AvdB_C = AvdB(K);
%     
% % =======================================================================   
% % TIME DOMAIN    select source frequency fs
% % =======================================================================
%     
%     KK = find(f > fs,1);
%     Ns = 1000;
%     ws = 2*pi*fs;
%     Ts = 1/fs;
%     tMin = 0;
%     tMax = 3*Ts;
%     t = linspace(tMin,tMax,Ns);
%  % input signal emf   
%     emf = real(V_IN .* exp(1j*ws*t));
%  % output or load voltage
%     vLoad = real(abs(V3(KK)) .* exp(1j*(ws*t + pi*phi_3(KK))));
%    
%     
% =======================================================================
%   GRAPHICS
% =======================================================================
 
 figure(1)   % Av vs f -----------------------------------------------
    pos = [0.07 0.05 0.25 0.28];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    set(gcf,'color','w');
    fss = 12;
    xP = f/1e3; yP = P_OUT.*1e3;        
    plot(xP, yP,'b','linewidth',2); 
%    hold on
%    xP = 1e-3.*[fC fC]; yP = [0 max(Av)];
%    plot(xP, yP,'r','linewidth',1);
    grid on
    xlabel('frequency  f  [ kHz ]')
    ylabel('p_{OUT} [ mW ]');
%    tm1 = 'f_C  =  ';
%    tm2 = num2str(fC, '%4.1f Hz \n');
%    tm = [tm1 tm2];
%    title(tm,'fontweight','normal');
    set(gca,'fontsize',fss)
%    

figure(2)   % Av vs f -----------------------------------------------
    pos = [0.33 0.05 0.25 0.28];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    set(gcf,'color','w');
    fss = 12;
    xP = f; yP = P_OUT.*1e3;        
    semilogx(xP, yP,'b','linewidth',2); 
%    hold on
%    xP = 1e-3.*[fC fC]; yP = [0 max(Av)];
%    plot(xP, yP,'r','linewidth',1);
    grid on
    xlabel('frequency  f  [ kHz ]')
    ylabel('p_{OUT} [ mW ]');
%    tm1 = 'f_C  =  ';
%    tm2 = num2str(fC, '%4.1f Hz \n');
%    tm = [tm1 tm2];
%    title(tm,'fontweight','normal');
    set(gca,'fontsize',fss)
%    

% figure(2)   % Av vs log(f) ----------------------------------------------
%    pos = [0.33 0.05 0.25 0.28];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    fss = 12;
%    xP = f; yP = AvdB;        
%    semilogx(xP, yP,'b','linewidth',2); 
%    hold on
%    xP = fC.*ones(length(yticks),1);
%    yP = yticks;
%    semilogx(xP, yP,'r','linewidth',1);
%    yP = AvdB_C.*ones(length(xticks),1);
%    xP = xticks;
%    semilogx(xP, yP,'r','linewidth',1);
%    
%    grid on
%    xlabel('f_{IN}  [ Hz ]');
%    ylabel('voltage gain  A_v  [dB]');
%    tm1 = 'f_C  =  ';
%    tm2 = num2str(fC, '%4.1f Hz \n');
%    tm = [tm1 tm2];
%    title(tm,'fontweight','normal');
%    set(gca,'fontsize',fss) 
%    
%    
% figure(3)   % output phase vs log10(f) ----------------------------------
%    pos = [0.59 0.05 0.25 0.28];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    fss = 12;
%    xP = f; yP = phi_OUT;        
%    semilogx(xP, yP,'b','linewidth',2);
%    hold on
%    xP = [fC fC];
%    yP = [-0.5 0];
%    if flagF == 2
%      yP = [0 0.5];
%    end
%    semilogx(xP, yP,'r','linewidth',1);
%    grid on
%    xlabel('f_{IN}  [ Hz ]');
%    ylabel('\phi_{OUT} / \pi');
%    tm1 = 'f_C  =  ';
%    tm2 = num2str(fC, '%4.1f Hz \n');
%    tm = [tm1 tm2];
%    title(tm,'fontweight','normal');
%    set(gca,'fontsize',fss)  
%    set(gca,'yTIck',-0.5:0.125:0);
%    if flagF == 2
%    set(gca,'yTIck',0:0.125:0.5);
%    end
%    
%  figure(4)   % output power  [mW] ----------------------------------
%    pos = [0.07 0.42 0.25 0.28];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    fss = 12;
%    xP = f; yP = P_OUT;        % capacitve reactance
%    semilogx(xP, yP,'b','linewidth',2);
%    hold on
%    xP = [fC fC]; yP = [0 max(P_OUT)];
%    semilogx(xP, yP,'r','linewidth',1);
%    
%    yP = max(P_OUT)/2.*ones(length(xticks),1);
%    xP = xticks;
%    semilogx(xP, yP,'r','linewidth',1);
%    
%    grid on
%    xlabel('f_{IN}  [ Hz ]');
%    ylabel('P_{OUT}  [ mW]');
%    tm1 = 'f_C  =  ';
%    tm2 = num2str(fC, '%4.1f Hz \n');
%    tm = [tm1 tm2];
%    title(tm,'fontweight','normal');
%    set(gca,'fontsize',fss)  
%    
%   figure(5)   % time domain: voltages  -------------------------------------
%    pos = [0.33 0.42 0.25 0.28];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    
%    xP = 1e3.*t; yP = emf;        
%    plot(xP, yP,'b','linewidth',2);   
%    hold on
%    yP = vLoad;
%    plot(xP, yP,'r','linewidth',2); 
%    legend(' emf',' v_{OUT}');
%    grid on
%    xlabel('time t  [ ms ]')
%    ylabel('voltages  V  [ V ]');
%    
%     tm1 = 'f_{IN}  =  ';
%     tm2 = num2str(fs, '%3.0f \n');
%     tm3 = ' Hz';
%     tm4 = '   \phi_{OUT} / \pi =  ';
%     tm5 = num2str(phi_OUT(KK), '%3.2f \n');
%     tm6 = ' ';
%     tm = [tm1 tm2 tm3 tm4 tm5 tm6];
%     h_tm = title(tm,'fontWeight','normal');
%     set(h_tm,'fontsize',fss);
%     set(gca,'fontsize',fss)      
%  
%    
%    