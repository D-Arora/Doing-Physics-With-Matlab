%CNsRCL.m

% SERIES RCL CIRCUIT
% Finite Difference Method for serires RCL circuit analysis
% Voltage inputs: step function OFF/ON or sinusoidal

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% 180131

clear all
close all
clc

% =====================================================================
% INPUTS: S.I. UNITS [default values] 
% ======================================================================
%  series resistance [40] 
    R  = 40;
%  capacitance  [1e-6]    
    C  = 1e-6;
%  inductance  [0.025330295910584]
    L  = 0.025330295910584 ;
%  peak value for source (input) emf   [10]  
    VS = 10;       % [10]
 
% Select voltage source
% flagV = 1: step function / flagV = 2: complex sinusoidal function 
% Change parameters for source within CASE script
   flagV = 2; 
  
% Resonance Frequency f0, period T and time step dt
   f0 = 1/(2*pi*sqrt(L*C));
   T0 = 1/f0;
   dt = T0 /1000;
   
switch flagV
  case 1  % Step function  OFF to ON
    % max simulation time [5*T] 
       %tMax = 5*T0;  
       tMax = 4e-3;
    % percentage off time  [10]
       pOFF = 10;         
                
         t = 0:dt:tMax;              % time grid               
         N = length(t);              % number of time steps
         nOFF = round(N*pOFF/100);   % number of grid point OFF
         vS = zeros(1,N);            % source emf
         vS(nOFF:end) = VS;
     
  case 2     % Complex sinusoidal input
    % source frequency  [1000]          
       fS = 2000;
    % max simulation time [5*T] 
     %  tMax = 10*T0;  
       tMax = 10e-3;
       
       w = 2*pi*fS;                      % angular frequency
       t = 0:dt:tMax;                    % time grid
       N = length(t);                    % number of time steps
       vS = VS .* exp(1j*(w*t - pi/2));  % source emf
     
    % Impedance calculations  
       XC = 1/(w*C); ZC = -1j*XC;
       XL = w*L;     ZL = 1j*XL;
                     ZR = R;
       Z = ZR + ZC + ZL;
  
end   % end switch
   
% =======================================================================  
% CALCULATIONS: Finite Difference Approximations
%               Initially capacitor is uncharged
% =======================================================================
% Initialise arrays: voltages and current iS = iR = iC = iL  
    vR = zeros(1,N);
    vC = zeros(1,N);
    vL = zeros(1,N);
    iS = zeros(1,N); 
  
%  Time Step #1
   vR(1)  = vS(1);
   vC(1)  = 0;            
   vL(1) = vS(1) - vR(1) - vC(1);
   iS(1) = vR(1) / R;
   iS(1) = iS(1) + vL(1) * dt / L;
   vC(1) = vC(1) + iS(1) * dt / C;
   vR(1) = iS(1) * R;
   vL(1) = vS(1) - vR(1) - vC(1);
  
% Time Steps #2 to #N
  for c = 2 : N
    iS(c) = iS(c-1) + vL(c-1) * dt/L;
    vR(c) = iS(c) * R;
    vC(c) = vC(c-1) + 0.5*(iS(c)+iS(c-1)) * dt / C;
    vL(c) = vS(c) - vR(c) - vC(c);
  end
 
% Peak values 
   VSp  = max(vS);
   VRp  = max(vR);
   VCp  = max(vC);
   VLp  = max(vL);
   ISp  = max(iS);

% Powers and energy
   pS = real(vS) .* real(iS);
   pR = real(vR) .* real(iS);
   pC = real(vC) .* real(iS);
   pL = real(vL) .* real(iS);
     
   uS = zeros(1,N); uR = zeros(1,N);
   uC = zeros(1,N); uL = zeros(1,N);
   for c = 2 : N
       uS(c) = uS(c-1) + pS(c)*dt;
       uR(c) = uR(c-1) + pR(c)*dt;
       uC(c) = uC(c-1) + pC(c)*dt;
       uL(c) = uL(c-1) + pL(c)*dt;
   end

% Phases voltage: phi and current theta  at time step nP  [degrees]
    nP = N-500;
    
    phiS   = rad2deg(angle(vS(nP)));
    phiR   = rad2deg(angle(vR(nP)));
    phiC   = rad2deg(angle(vC(nP)));
    phiL   = rad2deg(angle(vL(nP)));
    thetaS = rad2deg(angle(iS(nP)));
    phiSR = phiS - phiR;
    

% Find peaks in vS and corresponding times
%   Calculate period of oscillations
%   May need to change number of peaks for estimate of period: nPeaks
    [iS_Peaks, t_Peaks] = findpeaks(real(iS));
    nPeaks = 3;
    T_Peaks = (t(t_Peaks(end)) - t(t_Peaks(end-nPeaks)))/(nPeaks);
    f_Peaks = 1 / T_Peaks;
%   f_Peaks  = 0;
%   T_Peaks  = 0;
    
% =====================================================================
%   DISPLAY RESULTS
% =====================================================================
   disp('Circuit Elements  ') 
      fprintf(' Resistance   R  =  %3.2e  ohms  \n',R);
      fprintf(' Capacitance  C  =  %3.3e  F  \n',C);
      fprintf(' Inductance   L  =  %3.3e  H  \n',L);
   disp('  ')     
      if flagV == 1
        disp(' Source emf: step Function OFF/ON');
      else
        disp(' Source emf: Sinusoidal Function ');  
      end
      
      fprintf(' Peak emf    VS  =  %3.2e  V  \n',VS);
   disp('  ')   
      fprintf(' Resonance Frequency f0   =  %3.3e  Hz \n',f0);
      fprintf(' Resonance Period     T0  =  %3.3e  s  \n',T0);
      fprintf(' time increment      dt   = %3.3e   s  \n',dt);
      fprintf('                 dt / T   =  %3.2e    \n',dt/T0);
  
   if flagV == 2
      disp('  ');
      disp('Sinusoidal Source emf');
      fprintf('  Source frequency  fS  =  %3.3e  Hz  \n',fS);
      fprintf('  Impedance: resistance   ZR  =  %3.3e  ohms  \n',ZR);
      fprintf('  Impedance: capacitance  ZC  =  %3.3e %3.3e ohms  \n',real(ZC),imag(ZC));
      fprintf('    capacitance phase angle    phiC  =  %3.1f deg  \n',rad2deg(angle(ZC)));
      fprintf('  Impedance: inductance   ZL  =  %3.3e %3.3e ohms  \n',real(ZL),imag(ZL));
      fprintf('    inductance phase angle    phiL   =  %3.1f deg  \n',rad2deg(angle(ZL)));
      fprintf('  Impedance: total        Z   =  %3.3e %3.3e ohms  \n',real(Z),imag(Z));
      fprintf('  Impedance: magnitude    |Z| =  %3.2e   ohms  \n',abs(Z));
      fprintf('    impedance phase angle phiZ       =  %3.1f deg  \n',rad2deg(angle(Z)));
      disp('  ');
      fprintf(' Frequency of peaks f_Peaks  =  %3.3e  Hz \n',f_Peaks);
      fprintf(' Period of peaks    T_Peaks  =  %3.3e  s  \n',T_Peaks);
      disp('  ')
      fprintf(' Phases [degrees] at time t = %3.2f ms  \n',1e3*t(nP));
      fprintf(' phiS   = %3.0f deg  \n',phiS);
      fprintf(' phiR   = %3.0f deg  \n',phiR);  
      fprintf(' phiC   = %3.0f deg  \n',phiC);
      fprintf(' phiL   = %3.0f deg  \n',phiL);
      fprintf(' thetaS = %3.0f deg  \n',thetaS);
      fprintf(' Phase difference between source emf & current thetaSR = %3.0f deg  \n',phiSR);
   end 
   
% ======================================================================
% GRAPHICS
% ======================================================================

figure(1)  % voltages ------------------------------------------------ 
   FS = 14;
   pos = [0.07 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   subplot(3,1,1)
     xP = 1e3.*t; yP = real(vS);
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = real(vR);
     plot(xP,yP,'r','linewidth',2');
     ylabel('v_R  [ V ]');
     set(gca,'fontsize',FS);
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.0f');
     tm3 = ' \Omega   C = ';
     tm4 = num2str(C, '%3.1e');
     tm5 = ' F   L = ';
     tm6 = num2str(L, '%3.3e');
     tm7 = '  H';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',12);    
     tm1 = 'f_0 =  ';
     tm2 = num2str(f0, '%5.0f');
     tm3 = '  Hz';
     tm4 = '   f_{Peaks} = ';
     tm5 = num2str(f_Peaks, '%5.0f');
     tm = [tm1 tm2 tm3 tm4 tm5 tm3];
     xlabel(tm,'fontsize',12);
     
   subplot(3,1,2)
     xP = 1e3.*t; yP = real(vS);
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = real(vC);
     plot(xP,yP,'k','linewidth',2');
     ylabel('v_C  [ V ]');
     set(gca,'fontsize',FS);
     tm1 = 'T_0 =  ';
     tm2 = num2str(1e3*T0, '%3.2f');
     tm3 = '  ms';
     tm4 = '   T_{Peaks} = ';
     tm5 = num2str(1e3*T_Peaks, '%3.2f');
     tm = [tm1 tm2 tm3 tm4 tm5 tm3];
     xlabel(tm,'fontsize',12);
   subplot(3,1,3)
     xP = 1e3.*t; yP = real(vS);
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = real(vL);
     plot(xP,yP,'m','linewidth',2');
     xlabel('t  [ ms ]');
     ylabel('v_L  [ V ]');
     set(gca,'fontsize',FS);


figure(2) % powers --------------------------------------------------  
   FS = 14;
   pos = [0.37 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   subplot(3,1,1)
     xP = 1e3.*t; yP = pS;
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = pR;
     plot(xP,yP,'r','linewidth',2');
     if flagV == 1
       set(gca,'yLim',[-1 1]);
       set(gca,'yTick',-1:0.5:1);
     end
     ylabel('p_R  [ W ]');
     set(gca,'fontsize',FS);
     
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.0f');
     tm3 = ' \Omega   C = ';
     tm4 = num2str(C, '%3.1e');
     tm5 = ' F   L = ';
     tm6 = num2str(L, '%3.3e');
     tm7 = '  H';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',12);    
     
   subplot(3,1,2)
     xP = 1e3.*t; yP = pS;
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = pC;
     plot(xP,yP,'k','linewidth',2');
     if flagV == 1
       set(gca,'yLim',[-1 1]);
       set(gca,'yTick',-1:0.5:1);
     end
     ylabel('p_C  [ W ]');
     set(gca,'fontsize',FS);
     
   subplot(3,1,3)
     xP = 1e3.*t; yP = pS;
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = pL;
     plot(xP,yP,'m','linewidth',2')
     if flagV == 1
       set(gca,'yLim',[-1 1]);
       set(gca,'yTick',-1:0.5:1);
     end
     xlabel('t  [ ms ]');
     ylabel('p_L  [ W ]');
     set(gca,'fontsize',FS);
     

 figure(3)  % energies ------------------------------------------------ 
   FS = 14;
   pos = [0.67 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   subplot(3,1,1)
     xP = 1e3.*t; yP = 1e3.*uS;
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = 1e3.*uR;
     plot(xP,yP,'r','linewidth',2');
     ylabel('u_R  [ mJ ]');
     set(gca,'fontsize',FS);
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.0f');
     tm3 = ' \Omega   C = ';
     tm4 = num2str(C, '%3.1e');
     tm5 = ' F   L = ';
     tm6 = num2str(L, '%3.3e');
     tm7 = '  H';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',12);    
     
   subplot(3,1,2)
     xP = 1e3.*t; yP = 1e3.*uS;
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = 1e3.*uC;
     plot(xP,yP,'k','linewidth',2');
     ylabel('u_C  [ mJ ]');
     set(gca,'fontsize',FS);
     
   subplot(3,1,3)
     xP = 1e3.*t; yP = 1e3.*uS;
     plot(xP,yP,'b','linewidth',1');
     grid on
     hold on
     yP = 1e3.*uL;
     plot(xP,yP,'m','linewidth',2');
     xlabel('t  [ ms ]');
     ylabel('u_L  [ mJ ]');
     set(gca,'fontsize',FS);
     
 figure(4)   % current --------------------------------------------------
   LimY = 100;
   TickY = -LimY:20:LimY;
   FS = 14;
   pos = [0.40 0.15 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   yyaxis left 
   xP = 1e3.*t; yP = 1e3.*real(iS);
   plot(xP,yP,'r','linewidth',2');
   grid on
   xlabel('t  [ms]');
   ylabel('i_S  [ mA ]');
   ax = gca;
   ax.YAxis(1).Color = 'r';
   set(gca,'fontsize',FS);
   set(gca,'yLim',[-LimY, LimY]);
   set(gca,'yTick',TickY); 
   
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.0f');
     tm3 = ' \Omega   C = ';
     tm4 = num2str(C, '%3.1e');
     tm5 = ' F   L = ';
     tm6 = num2str(L, '%3.3e');
     tm7 = '  H';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',12);      
    
  yyaxis right
     xP = 1e3.*t; yP = real(vS);
     plot(xP,yP,'b','linewidth',2');
     ylabel('v_S  [ V ]');
     ax.YAxis(2).Color = 'b';
     set(gca,'yLim',[-1, 11]);
   set(gca,'yTick',0:2:10);
     
    
 if flagV == 2        
  figure(5)  % V/I characteristics   ------------------------------------
   FS = 14;
   pos = [0.10 0.15 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
  subplot(3,1,1)
     xP = real(vR); yP = 1e3.*real(iS);
     plot(xP,yP,'b','linewidth',1');
     grid on
     xlabel('v_R  [ V ]');
     ylabel('i_R  [ mA ]');
     
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.0f');
     tm3 = ' \Omega   C = ';
     tm4 = num2str(C, '%3.1e');
     tm5 = ' F   L = ';
     tm6 = num2str(L, '%3.3e');
     tm7 = '  H';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',12);    
     set(gca,'fontsize',FS);
    % set(gca,'xLim',[0 1e3*max(t)]);
    % set(gca,'yLim',[-VS-0.2 VS+0.2]);
    % set(gca,'yTick',-10:5:10);

  subplot(3,1,2)
     xP = real(vC); 
     plot(xP,yP,'b','linewidth',1');
     grid on
     xlabel('v_C  [ V ]');
     ylabel('i_C  [ mA ]');
     set(gca,'fontsize',FS);
     
     subplot(3,1,3)
     xP = real(vL); 
     plot(xP,yP,'b','linewidth',1');
     grid on
     xlabel('v_L  [ V ]');
     ylabel('i_L  [ mA ]');
     set(gca,'fontsize',FS);
             
%%
figure(6)   % phasors ------------------------------------------------
   xyLim = 20;
   xyTick = -xyLim:5:xyLim;
   
   FS = 14;
   pos = [0.50 0.15 0.28 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = [0 real(vR(nP))]; yP = [0 imag(vR(nP))];
   plot(xP,yP,'r','linewidth',2');
   hold on
   xP = [0 real(vC(nP))]; yP = [0 imag(vC(nP))];
   plot(xP,yP,'k','linewidth',2'); 
   xP = [0 real(vL(nP))]; yP = [0 imag(vL(nP))];
   plot(xP,yP,'m','linewidth',2');  
   xP = [0 real(vS(nP))]; yP = [0 imag(vS(nP))];
   plot(xP,yP,'b','linewidth',2');
  % vA = vR+vC+vL;
  % xP = [0 real(vA(nP))]; yP = [0 imag(vA(nP))];
  % plot(xP,yP,'g','linewidth',1');
   xlabel('Re(V)  [ V ]');
   ylabel('Im(V)  [ V ]');
   set(gca,'fontsize',FS);
   set(gca,'xLim',[-xyLim, xyLim]);
   set(gca,'yLim',[-xyLim, xyLim]);
   set(gca,'xTick',xyTick);
   set(gca,'yTick',xyTick);
   grid on
   axis square
  % axis equal
  
     tm1 = 'Phasors at time t  =   ';
     tm2 = num2str(1e3*t(nP), '%3.2f');
     tm3 = ' ms';
     tm = [tm1 tm2 tm3];
     title(tm,'fontweight','normal','fontsize',12);  
end
