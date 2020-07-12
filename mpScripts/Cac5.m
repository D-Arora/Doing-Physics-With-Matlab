% Cac51.m

% RLC series circuits

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180103 


clear all
close all
clc

% ========================================================================
%   INPUTS   SI units default values [ ]
% =======================================================================
% source emf (peak or amplidue value) [10 10e3]
   VS = 10;        
   f = 15.91549e3;
% resistance Z1    [1e3]
   R = 1e3;
% capacitance Z2 [1.0e-8 F]
   C = 1.0e-8;   
% inductance Z1 [10e-3 H]
   L = 10e-3;
 
% =======================================================================  
%   CALCULATIONS   SI units
% =======================================================================   
   w = 2*pi*f;                   % angular frequency
   T = 1/f;                      % period 
   t = linspace(0,3*T,5001);     % time  N must be an odd number for Simpson's Rule
   vS = VS .* exp(1i*w*t);       % emf as a function of time
   f0 = 1/(2*pi*sqrt(L*C));       % Resonance frequency
% Impedances
   Z1 = R;                       % resistance
   X2 = 1/(w*C);                 % capacitive reactance
   Z2 = -1i * X2;                % capacitive impedance 
   X3 = w*L;                     % inductive reactance
   Z3 = 1i * X3;                 % inductive impedance 
   Z4 = Z1 + Z2 + Z3;            % total circuit impedance
  
% Currents and Voltages 
   iS = vS ./ Z4;               % source current
   i1 = iS; i2 = iS; i3 = iS;    % element currents
   theta = angle(iS);            % current phase
   IS = max(abs(iS));            % peak current
   
   v1 = i1 .* Z1;                % element voltages
   v2 = i2 .* Z2;
   v3 = i3 .* Z3;
   
   V1 = max(abs(v1));            % peak voltages
   V2 = max(abs(v2));
   V3 = max(abs(v3));
   
   phi1 = angle(v1);             % voltage phases 
   phi2 = angle(v2);
   phi3 = angle(v3);
   phiS = 0;

   pS = real(vS) .* real(iS);   % powers
   p1 = real(v1) .* real(i1);
   p2 = real(v2) .* real(i2);
   p3 = real(v3) .* real(i3);
   
   Prms = 0.5*real((vS(1) .* conj(iS(1))));
   Prms_N = simpson1d(pS,0,t(end)/t(end));
   P1 = simpson1d(p1,0,t(end)/t(end));
   P2 = simpson1d(p2,0,t(end)/t(end));
   P3 = simpson1d(p3,0,t(end)/t(end));

% =======================================================================
% DISPLAY RESULTS  actual (real) values
% =======================================================================
   disp('Inputs   ');
     fprintf(' Source peak voltage VS =  %3.2f  V \n',VS);
     fprintf(' Source frequency f  =  %3.2e Hz \n',f);
     fprintf('  R =  %3.2f ohms \n',R);
     fprintf('  C =  %3.2e F \n',C);
     fprintf('  L =  %3.2e H \n',L);
   disp('  ') 
   disp('Outputs   ');
     fprintf(' Resonance freq  f0 =  %3.2f  Hz  \n',f0);
     fprintf('  XC =  %3.2f  ohms  \n',X2);
     fprintf('  XL =  %3.2f  ohms  \n',X3);
     fprintf('  peak current  IS =  %3.2f  \n',1e3*IS);
   disp('Peak Values')
     fprintf('  IS =  %3.2f mA \n',1e3*IS);
     fprintf('  emf  VS =  %3.2f V \n',VS);
     fprintf('  VR =  %3.2f V \n',V1);
     fprintf('  VC =  %3.2f V \n',V2);
     fprintf('  VL =  %3.2f V \n',V3);
   disp('Phases')
     fprintf('  phi_S =  %3.2f pi rad\n',phiS(1)/pi);
     fprintf('  phi_R =  %3.2f pi rad \n',phi1(1)/pi);
     fprintf('  phi_C =  %3.2f pi rad \n',phi2(1)/pi);
     fprintf('  phi_L =  %3.2f pi rad \n',phi3(1)/pi);
   disp('Power  rms values')
     fprintf('  Prms =  %3.2f pi mW  \n',1e3*Prms);
     fprintf('  Prms (Simpsons Rule) Prms =  %3.2f mW \n',1e3*Prms_N);   
     fprintf('  Simpsons Rule PR =  %3.2f mW \n',1e3*P1);
     fprintf('  Simpsons Rule PC =  %3.2f mW \n',1e3*P2);
     fprintf('  Simpsons Rule PL =  %3.2f mW \n',1e3*P3);
     
     
% =======================================================================
%   GRAPHICS
% =======================================================================
   
figure(1)  % source voltage and current 
   FS = 14;
   pos = [0.07 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   yyaxis left
   xP = t; yP = real(vS);
   plot(xP,yP,'b','linewidth',2');
   set(gca,'xLim',[0 max(t)]);
   grid on
   xlabel('time  t  [ms]');
   ylabel('emf  [ V ]');
   tm  = sprintf('phase difference / \\pi    \\Delta\\phi = %3.2f \n',theta(1)/pi);
   title(tm,'fontweight','normal');
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',FS);
   
   
   yyaxis right
   xP = t; yP = real(iS).*1e3;
   plot(xP,yP,'r','linewidth',2');
   ylabel('i_S  [ mA ]');
   ax = gca;
   ax.YAxis(2).Color = 'r';
   
figure(2)   % voltages --------------------------------------------------  
   FS = 14;
   pos = [0.07 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   xP = t; yP = real(vS);
   plot(xP,yP,'b','linewidth',2');
   hold on
   yP = real(v1);
   plot(xP,yP,'r','linewidth',2');
   yP = real(v2);
   plot(xP,yP,'k','linewidth',2');
   yP = real(v3);
   plot(xP,yP,'m','linewidth',2');
   set(gca,'xLim',[0 max(t)]);
   legend(' emf',' v_R',' v_C', ' v_L');
   grid on
   xlabel('time  t  [ms]');
   ylabel('voltage  [ V ]');
   set(gca,'fontsize',FS);
   
figure(3)   % powers --------------------------------------------------  
   FS = 14;
   pos = [0.38 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   xP = t; yP = 1e3.*pS;
   plot(xP,yP,'b','linewidth',2');
   hold on
   yP = 1e3.*p1;
   plot(xP,yP,'r','linewidth',2');
   yP = 1e3.*p2;
   plot(xP,yP,'k','linewidth',2');
   yP = 1e3.*p3;
   plot(xP,yP,'m','linewidth',2');
   set(gca,'xLim',[0 max(t)]);
   legend(' P_S',' P_R',' P_C', ' P_L');
   grid on
   xlabel('time  t  [ms]');
   ylabel('power  [mW ]');
   set(gca,'fontsize',FS);  

   figure(4)   % phasors --------------------------------------------------  
   FS = 14;
   pos = [0.38 0.48 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   xP = [0 cos(phiS(1))]; yP = [0 sin(phiS(1))];
   plot(xP,yP,'b','linewidth',2');
   hold on
   xP = (V1/VS)*[0 cos(phi1(1))]; yP = [0 (V1/VS)*sin(phi1(1))];
   plot(xP,yP,'r','linewidth',2')
   xP = (V2/VS)*[0 cos(phi2(1))]; yP = (V2/VS)*[0 sin(phi2(1))];
   plot(xP,yP,'k','linewidth',2')
   xP = (V3/VS)*[0 cos(phi3(1))]; yP = (V3/VS)*[0 sin(phi3(1))];
   plot(xP,yP,'m','linewidth',2')
   legend(' emf',' v_R',' v_C', ' v_L','Location','EastOutside');
   grid on
   xlabel('Re');
   ylabel('Im');
   set(gca,'xlim',[-1.1 1.1]);
   set(gca,'ylim',[-1.1 1.1]);
   set(gca,'xTick',-1:0.5:1);
   set(gca,'xTick',-1:0.5:1);
   axis square
   tm = 'Voltage Phasors';
   title(tm,'fontweight','normal');
   set(gca,'fontsize',FS);  