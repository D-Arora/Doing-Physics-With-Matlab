%CN01.m

% V/I connections for a resistor, capacitor and inductor
% Finite Difference Method for circuit anlysis

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% 180121

clear all
close all
clc

% INPUTS: S.I. UNITS ====================================================
  R = 1.0e3;
  C = 1e-6;
  L = 10e-3;
  IR = 0.10;
  IC = 0.10;
  VL = 10;
  f = 1e3;
  N = 200;

% CALCULATIONS ==========================================================  
  T = 1/f;
  w = 2*pi*f;
  tMax = 2*T;
  t = linspace(0,tMax,N);
  dt = t(2)-t(1);
  
  iR = IR.*cos(w*t);
  vR = R .* iR;
  
  kC = 2*dt/C;
  iC = IC.* cos(w*t);
  vC = zeros(1,N);
  vC(2) = vC(1)+(dt/C)*iC(1);
  
  kL = 2*dt/L;
  vL = VL .* sin(w*t);
  iL = zeros(1,N);
  iL(2) = iL(1)+(dt/L)*vL(1);

for c = 3 : N-1
 vC(c) = vC(c-2) + kC*iC(c-1);
 iL(c) = iL(c-2) + kL*vL(c-1);
 %iL(c) = iL(c-1)+(dt/L)*vL(c-1);
end
 
 iL = iL - max(iL)/2;   % not sure why need to shift the current
 
% Peak values /Impedances / reactances
  VR = max(vR);
  IR = max(iR);
  VC = max(vC);
  IC = max(iC);
  VL = max(vL);
  IL = max(iL);

  ZR = VR/IR;        % resistance
  ZC = VC/IC;        % capacitive reactance
  ZL = VL/IL;        % inductive reactance

% Theorteical values for impedance
  ZR_T = VR/IR;
  ZC_T = 1/(w*C);
  ZL_T = w*L;


% =====================================================================
%   DISPLAY RESULTS
% =====================================================================
     fprintf(' input frequency    f =  %3.2f   Hz \n',f);
     disp('  ') 
disp(' RESISTOR  ');
     fprintf(' resistance   R =  %3.2f   ohms \n',R);
     fprintf(' peak voltage VR =  %3.2f  V  \n',VR);
     fprintf(' peak current IR =  %3.2f  mA  \n',1e3*IR);
     fprintf(' impedance    ZR =  %3.2f  ohms  \n',ZR);
     fprintf(' theoretical impedance ZR_T =  %3.2f ohms \n',ZR_T);
     disp('  ') 
disp('  CAPACITOR  ');
     fprintf(' capacitance  C =  %3.2e  F \n',C);
     fprintf(' peak voltage VC =  %3.2f  V  \n',VC);
     fprintf(' peak current IC =  %3.2f  mA  \n',1e3*IC);
     fprintf(' impedance    ZC =  %3.2f  ohms  \n',ZC);
     fprintf(' theoretical impedance ZC_T =  %3.2f ohms \n',ZC_T);
     disp('  ') 
disp('  INDUCTOR  ');
     fprintf(' inductance   L =  %3.2e  H \n',L);
     fprintf(' peak voltage VL =  %3.2f  V  \n',VL);
     fprintf(' peak current IL =  %3.2f  mA  \n',1e3*IL);
     fprintf(' impedance    ZL =  %3.2f  ohms  \n',ZL);
     fprintf(' theoretical impedance ZL_T =  %3.2f ohms \n',ZL_T);
     
% ======================================================================
% GRAPHICS
% ======================================================================

figure(1)    % RESISTOR ------------------------------------------------
   FS = 14;
   pos = [0.07 0.05 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   yyaxis left 
   xP = 1e3.*t; yP = vR;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('time  t  [ms]');
   ylabel('v_R  [ V ]');
   
   
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',FS);
   set(gca,'xLim',[0 1e3*max(t)]);
   
   yyaxis right
   xP = 1e3.*t; yP = iR.*1e3;
   plot(xP,yP,'r--','linewidth',2');
   ylabel('i_R  [ mA ]');
   ax = gca;
   ax.YAxis(2).Color = 'r';
   
   tm  = ' RESISTOR: voltage and current in phase'; 
   title(tm,'fontweight','normal','fontsize',12);


figure(2);  % CAPACITOR -------------------------------------------------
   FS = 14;
   pos = [0.37 0.05 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   yyaxis left 
   xP = 1e3.*t; yP = vC;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('time  t  [ms]');
   ylabel('v_C  [ V ]');
    
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',FS);
   set(gca,'xLim',[0 1e3*max(t)]);
   
   yyaxis right
   xP = 1e3.*t; yP = 1e3.*iC;
   plot(xP,yP,'r','linewidth',2');
   ylabel('i_C  [ mA ]');
   ax = gca;
   ax.YAxis(2).Color = 'r';
   
   tm  = ' CAPACITOR: voltage phase lags current phase by \phi = -\pi / 2 rad'; 
   title(tm,'fontweight','normal','fontsize',12);

figure(3);   % INDUCTOR --------------------------------------------------
   FS = 14;
   pos = [0.67 0.05 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   yyaxis left 
   xP = 1e3.*t; yP = vL;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('time  t  [ms]');
   ylabel('v_L  [ V ]');
    
   ax = gca;
   ax.YAxis(1).Color = 'b';
   set(gca,'fontsize',FS);
   set(gca,'xLim',[0 1e3*max(t)]);
   
   yyaxis right
   xP = 1e3.*t; yP = 1e3.*iL;
   plot(xP,yP,'r','linewidth',2');
   ylabel('i_L  [ mA ]');
   ax = gca;
   ax.YAxis(2).Color = 'r';
   
   tm  = ' INDUCTOR: voltage phase leads current phase by \phi = +\pi / 2 rad'; 
   title(tm,'fontweight','normal','fontsize',12);
   
figure(4)   % RESISTOR I/V curve ----------------------------------------
   FS = 14;
   pos = [0.07 0.45 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
        
   xP = vR; yP = 1e3.* iR;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('v_R [ V ]');
   ylabel('I_R  [ mA ]');
   set(gca,'fontsize',FS);
   %set(gca,'xLim',[0 1e3*max(t)]);
   tm = 'RESISTOR';
   title(tm,'fontweight','normal','fontsize',12);
   
 figure(5)   % CAPACITOR I/V curve ----------------------------------------
   FS = 14;
   pos = [0.37 0.45 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
        
   xP = vC; yP = 1e3.* iC;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('v_C [ V ]');
   ylabel('I_C  [ mA ]');
   set(gca,'fontsize',FS);
   %set(gca,'xLim',[0 1e3*max(t)]);
   tm = 'CAPACITOR';
   title(tm,'fontweight','normal','fontsize',12);  
   
figure(6)   % INDUCTOR I/V curve ----------------------------------------
   FS = 14;
   pos = [0.67 0.45 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
        
   xP = vL; yP = 1e3.* iL;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('v_L [ V ]');
   ylabel('I_L  [ mA ]');
   set(gca,'fontsize',FS);
   %set(gca,'xLim',[0 1e3*max(t)]);
   tm = 'INDUCTOR';
   title(tm,'fontweight','normal','fontsize',12);
   
   