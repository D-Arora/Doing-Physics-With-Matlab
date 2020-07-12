%CNpRCL.m

% PARALLEL RCL CIRCUIT
% Finite Difference Method for PARALLEL RCL circuit analysis
% Voltage inputs: step function OFF/ON or sinusoidal
% Calls the function arrow.m to draw arrows 


% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180207

clear all
close all
clc

% =====================================================================
% INPUTS: S.I. UNITS [default values] 
% ======================================================================
%  series resistance [80] 
    R  = 80;
%  capacitance  [1e-6]    
    C  = 1e-6;
%  inductance  [3.8e-3]
    L  = 0.5e-3;
%  Inductor resistance  [1]    min value 0.0001
    RL = 30;
%  Output (load) resi stance [1e4]
    RO = 80;
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
   
   if RL < 0.0001; RL = 0.0001; end
   
switch flagV
  case 1  % Step function  OFF to ON
    % INPUT max simulation time [5*T] 
       %tMax = 5*T0;  
       tMax = 1e-3;
    % INPUT percentage off time  [10]
       pOFF = 10;         
                
         t = 0:dt:tMax;              % time grid               
         N = length(t);              % number of time steps
         nOFF = round(N*pOFF/100);   % number of grid point OFF
         vS = zeros(1,N);            % source emf
         vS(nOFF:end) = VS;
     
  case 2     % Complex sinusoidal input
    % INPUT source frequency  [f0]          
       fS = 2000;
    % INPUT max simulation time [5*T] 
    %  tMax = 10*T0;  
       tMax = 1e-3;
        
       w = 2*pi*fS;                      % angular frequency
       t = 0:dt:tMax;                    % time grid
       N = length(t);                    % number of time steps
       vS = VS .* exp(1j*(w*t - pi/2));  % source emf
     
    % Impedance calculations  
       XC = 1/(w*C); ZC = -1j*XC;
       XL = w*L;     ZL = 1j*XL;
                     ZR = R;
       
  
end   % end switch
   

% =======================================================================  
% CALCULATIONS: Finite Difference Approximations
%               Initially capacitor is uncharged
% =======================================================================
% Initialise arrays: voltages and current
%  vS = vR + vO    v0 = vC = vL + vRL
%  iS = iR = iC + iL + iRL    iO = iRL
    vR = zeros(1,N);
    vC = zeros(1,N);
    vL = zeros(1,N);
    vRL= zeros(1,N);
    iS = zeros(1,N); 
    iC = zeros(1,N);
    iL = zeros(1,N);
    iO = zeros(1,N);
    
% Time Step #1
   vR(1)  = vS(1);
   iS(1) = vR(1) / R;
   iC(1) = iS(1);
  
% Time Steps #2 to #N
  for c = 2 : N
     vC(c)  = vC(c-1) + iC(c-1) * dt/C;
     vR(c)  = vS(c) - vC(c);
     iS(c)  = vR(c) / R;
     iO(c)  = vC(c) / RO;
     vRL(c) = (iL(c-1) + vL(c-1)*dt/L)*RL;
     vL(c)  = vC(c) - vRL(c);
     iL(c)  = vRL(c)/RL;
     iC(c)  = iS(c) - iO(c) - iL(c);
  end
 
  vO = vC;

% Peak values 
    VSp = max(vS);
    VRp = max(vR);
    VCp = max(vC);
    VLp = max(vL);
    ISp = max(iS);
    ICp = max(iC);
    IOp = max(iO);
% 
% Powers 
    pS = real(vS) .* real(iS);
    pO = real(vO) .* real(iO);

% Phases voltage: phi and current theta  at time step nP  [degrees]
% INPUT value for nP
     nP = N-500;
     
     phi_vS = rad2deg(angle(vS(nP)));
     phi_vO = rad2deg(angle(vO(nP)));
     theta_iS = rad2deg(angle(iS(nP)));
     theta_iO = rad2deg(angle(iO(nP)));
     theta_iC = rad2deg(angle(iC(nP)));
     theta_iL = rad2deg(angle(iL(nP)));
     
% Total circuit impedance 
     Z = (vS(nP)/iS(nP));
     
      
% Find peaks in vS and corresponding times
%   Calculate period of oscillations
%   May need to change number of peaks for estimate of period: nPeaks
     if flagV == 2
       [iS_Peaks, t_Peaks] = findpeaks(real(iS));
       nPeaks = length(t_Peaks)-1;
       T_Peaks = (t(t_Peaks(end)) - t(t_Peaks(end-nPeaks)))/(nPeaks);
       f_Peaks = 1 / T_Peaks;
     end
    
     
% ======================================================================
% GRAPHICS
% ======================================================================
 
 figure(1)  % voltages ------------------------------------------------ 
   FS = 14;
   pos = [0.07 0.05 0.28 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   xP = 1e3.*t; yP = real(vS);
   plot(xP,yP,'b','linewidth',1');
   grid on
   hold on
   yP = real(vO);
   plot(xP,yP,'r','linewidth',2');
   set(gca,'fontsize',FS);
   xlabel('t  [ms]','fontsize',12);
   ylabel('v  [ V ]','fontsize',12);
   legend('v_S','v_O');
   tm1 = 'f_S = ';
   
   tm3 = '  kHz';
   tm4 = '    f_0 = ';
   tm5 = num2str(f0/1e3, '%3.3f');
   tm = [tm4 tm5 tm3];
   if flagV == 2
     tm2 = num2str(fS/1e3, '%3.3f');  
     tm = [tm1 tm2 tm3 tm4 tm5 tm3];
   end
   title(tm,'fontweight','normal','fontsize',12); 
   set(gca,'fontsize',FS);


 figure(2)  % voltages ------------------------------------------------ 
   FS = 14;
   pos = [0.37 0.05 0.28 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   xP = 1e3.*t; yP = real(vS);
   plot(xP,yP,'b','linewidth',1');
   grid on
   hold on
   yP = real(vO);
   plot(xP,yP,'k','linewidth',2');
   yP = real(vL);
   plot(xP,yP,'m','linewidth',2');
   xlabel('t  [ms]','fontsize',12);
   ylabel('v  [ V ]','fontsize',12);
   legend('v_S','v_C','v_L');
   tm1 = 'f_S = ';
   tm3 = '  kHz';
   tm4 = '    f_0 = ';
   tm5 = num2str(f0/1e3, '%3.3f');
   tm = [tm4 tm5 tm3];
   if flagV == 2
     tm2 = num2str(fS/1e3, '%3.3f');
     tm = [tm1 tm2 tm3 tm4 tm5 tm3];
   end
   title(tm,'fontweight','normal','fontsize',12); 
   set(gca,'fontsize',FS);

   
 
   figure(3)  % TEXT------------------------------------------------ 
   %FS = 14;
   pos = [0.27 0.05 0.35 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   tm1 = 'R  = ';
   tm2 = num2str(R,'%3.1e');
   tm3 = '  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(0,100,tm);
   set(h_text,'fontsize',14);
   
   tm1 = 'R_L = ';
   tm2 = num2str(RL,'%3.1f');
   tm3 = '  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(0,85,tm);
   set(h_text,'fontsize',14);
   
   tm1 = 'R_O = '; 
   tm2 = num2str(RO,'%3.1e');
   tm3 = '  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(0,70,tm);
   set(h_text,'fontsize',14); 
   
   tm1 = 'C  = ';
   tm2 = num2str(C,'%3.2e');
   tm3 = '  F';
   tm = [tm1 tm2 tm3];
   h_text = text(0,55,tm);
   set(h_text,'fontsize',14); 
   
   tm1 = 'L  = ';
   tm2 = num2str(L,'%3.2e');
   tm3 = '  H';
   tm = [tm1 tm2 tm3];
   h_text = text(0,40,tm);
   set(h_text,'fontsize',14); 
   
   tm1 = 'V_S = ';
   tm2 = num2str(VS,'%3.2f');
   tm3 = '  V';
   tm = [tm1 tm2 tm3];
   h_text = text(0,25,tm);
   set(h_text,'fontsize',14); 
   
   tm1 = 'f_0 = ';
   tm2 = num2str(f0/1e3,'%3.3e');
   tm3 = '  kHz';
   tm = [tm1 tm2 tm3];
   h_text = text(0,10,tm);
   set(h_text,'fontsize',14); 
   
 if flagV == 2
   x = 36; y = 100; dy = -15;   
   tm1 = 'f_S = ';
   tm2 = num2str(fS/1e3,'%3.3e');
   tm3 = '  kHz';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14); 
   
   y = y+dy;
   tm1 = 'f_{Peaks} = ';
   tm2 = num2str(f_Peaks,'%3.3e');
   tm3 = '  kHz';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14); 
   
   y = y+dy;
   AV = max(real(vO))/max(real(abs(vS)));
   tm1 = 'A_V = ';
   tm2 = num2str(AV,'%3.2f');
   tm3 = '';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
  
   y = y+dy;
   tm1 = 'Z_C = ';
   tm2 = num2str(imag(ZC),'%3.2f');
   tm3 = ' j   \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
     
   y = y+dy;
   tm1 = 'Z_L = ( ';
   tm2 = num2str(RL,'%3.2f');
   tm3 = ' + ';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = ' ';
   tm2 = num2str(imag(ZL),'%3.2f');
   tm3 = ' j )  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x+10,y+10,tm);
   set(h_text,'fontsize',14);
   
   %y = y+dy;
   tm1 = 'Z = ( ';
   tm2 = num2str(real(Z),'%3.2f');
   tm3 = ' + ';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = ' ';
   tm2 = num2str(imag(Z),'%3.2f');
   tm3 = ' j )  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x+10,y+10,tm);
   set(h_text,'fontsize',14); 
   
   y = y+dy+10;
   tm1 = '|Z| = ';
   tm2 = num2str(abs(Z),'%3.2f');
   tm3 = ' \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   
   x = 78; y = 100;
   tm1 = 'At time t = ';
   tm2 = num2str(t(nP)*1e3,'%3.2f');
   tm3 = '  ms';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = '\phi_S = ';
   tm2 = num2str(phi_vS,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14); 
   
   y = y + dy;
   tm1 = '\phi_O = ';
   tm2 = num2str(phi_vO,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y + dy;
   tm1 = '\theta_S = ';
   tm2 = num2str(theta_iS,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y + dy;
   tm1 = '\theta_O = ';
   tm2 = num2str(theta_iO,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y + dy;
   tm1 = '\theta_C = ';
   tm2 = num2str(theta_iC,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y + dy;
   tm1 = '\theta_L = ';
   tm2 = num2str(theta_iL,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
 end
 
   axis off
   set(gca,'xLim',[0 100]);
   set(gca,'yLim',[0 100]);
   set(gca,'Position', [0.05 0.1 0.8 0.85]);
   
   figure(4)  % powers ------------------------------------------------ 
   FS = 14;
   pos = [0.07 0.45 0.28 0.38];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   subplot(2,1,1)
   xP = 1e3.*t; yP = pS;
   plot(xP,yP,'b','linewidth',1');
   grid on
   xlabel('t  [ms]','fontsize',12);
   ylabel('p_S  [ W ]','fontsize',12);
   tm1 = 'f_S = ';
   tm3 = '  kHz';
   tm4 = '    f_0 = ';
   tm5 = num2str(f0/1e3, '%3.3f');
   tm = [tm4 tm5 tm3];
   if flagV == 2
     tm2 = num2str(fS/1e3, '%3.3f');
     tm = [tm1 tm2 tm3 tm4 tm5 tm3];
   end
   title(tm,'fontweight','normal','fontsize',12); 
   set(gca,'fontsize',FS);

   subplot(2,1,2)
   yP = pO;
   plot(xP,yP,'r','linewidth',2');
   xlabel('t  [ms]','fontsize',12);
   ylabel('p_O  [ W ]','fontsize',12);
   grid on
   set(gca,'fontsize',FS);
   

   figure(5)  % currents ------------------------------------------------ 
   FS = 14;
   pos = [0.47 0.15 0.28 0.68];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   subplot(4,1,1)
   xP = 1e3.*t; yP = real(iS).*1e3;
   plot(xP,yP,'b','linewidth',1');
   grid on
   %xlabel('t  [ms]','fontsize',12);
   ylabel('i_S  [ mA ]','fontsize',12);
   tm1 = 'f_S = ';
   tm3 = '  kHz';
   tm4 = '    f_0 = ';
   tm5 = num2str(f0/1e3, '%3.3f');
   tm = [tm4 tm5 tm3];
   if flagV == 2
     tm2 = num2str(fS/1e3, '%3.3f');
     tm = [tm1 tm2 tm3 tm4 tm5 tm3];
   end
   title(tm,'fontweight','normal','fontsize',12); 
   set(gca,'fontsize',FS);

   subplot(4,1,2)
   yP = real(iO).*1e3;
   plot(xP,yP,'r','linewidth',2');
   xlabel('t  [ms]','fontsize',12);
   ylabel('i_O  [ mA ]','fontsize',12);
   grid on
   set(gca,'fontsize',FS);
   
   subplot(4,1,3)
   yP = real(iC).*1e3;
   plot(xP,yP,'k','linewidth',2');
   %xlabel('t  [ms]','fontsize',12);
   ylabel('iC  [ mA ]','fontsize',12);
   grid on
   set(gca,'fontsize',FS);
   
   subplot(4,1,4)
   yP = real(iL).*1e3;
   plot(xP,yP,'m','linewidth',2');
   xlabel('t  [ms]','fontsize',12);
   ylabel('i_L  [ mA ]','fontsize',12);
   grid on
   set(gca,'fontsize',FS);
   
   
%%  
figure(6)
   FS = 14;
   pos = [0.67 0.45 0.25 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   xP = [0 0]; yP = [-10,37];
   plot(xP,yP,'k','linewidth',1);
   hold on
   
   a = 8; xP = 0;yP = 45;
   th = linspace(0,2*pi,200);
   xunit = a * cos(th) + xP;
   yunit = a * sin(th) + yP;
   h = plot(xunit, yunit,'k','linewidth',2);
   
   xP = [0 0]; yP = [53,90];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [0 110]; yP = [-10,-10];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [0 30 30]; yP = [90 90 80];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [25 35]; yP = [80 80];
   plot(xP,yP,'k','linewidth',2);
   
   xP = [25 35]; yP = [65 65];
   plot(xP,yP,'k','linewidth',2);
   xP = [25 25]; yP = [65 80];
   plot(xP,yP,'k','linewidth',2);
   xP = [35 35]; yP = [65 80];
   plot(xP,yP,'k','linewidth',2);
   xP = [30 30]; yP = [29 65];
   plot(xP,yP,'k','linewidth',1);
   
   % C
   xP = [22 38]; yP = [29 29];
   plot(xP,yP,'k','linewidth',2);
   xP = [22 38]; yP = [25 25];
   plot(xP,yP,'k','linewidth',2);
   
   xP = [30 30]; yP = [-10 25];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [30 110]; yP = [55 55];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [85 85]; yP = [35 55];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [85 85]; yP = [-10 20];
   plot(xP,yP,'k','linewidth',1);
   
   % RO
   xP = [80 90]; yP = [35 35];
   plot(xP,yP,'k','linewidth',2);
   xP = [80 90]; yP = [20 20];
   plot(xP,yP,'k','linewidth',2);
   xP = [80 80]; yP = [20 35];
   plot(xP,yP,'k','linewidth',2);
   xP = [90 90]; yP = [20 35];
   plot(xP,yP,'k','linewidth',2);
   
   % Inductor branch
   xP = [55 55]; yP = [45 55];
   plot(xP,yP,'k','linewidth',1);
   
   xP = [50 60]; yP = [45 45];
   plot(xP,yP,'k','linewidth',2);
   xP = [50 60]; yP = [30 30];
   plot(xP,yP,'k','linewidth',2);
   xP = [50 50]; yP = [30 45];
   plot(xP,yP,'k','linewidth',2);
   xP = [60 60]; yP = [30 45];
   plot(xP,yP,'k','linewidth',2);
   
   xP = [55 55]; yP = [25 30];
   plot(xP,yP,'k','linewidth',1);
   
   % inductor
   th = linspace(-pi/2,pi/2,200);
   a = 5; 
   aX = a * cos(th);
   aY = a * sin(th);
   xP = 55+aX; yP = 20+aY;
   plot(xP,yP,'k','linewidth',2);
   xP = 55+aX; yP = 10+aY;
   plot(xP,yP,'k','linewidth',2);
   xP = 55+aX; yP = 0+aY;
   plot(xP,yP,'k','linewidth',2);
   
   xP = [55 55]; yP = [-10 -5];
   plot(xP,yP,'k','linewidth',1);
   
   % TEXT
   text(-3,45,'v_S','fontsize',14);
   text(13,27,'C','fontsize',14);
   text(13,70,'R','fontsize',14);
   text(48,10,'L','fontsize',14);
   text(38,37,'R_L','fontsize',14);
   text(67,27,'R_O','fontsize',14);
   text(108,27,'v_O','fontsize',14);
   arrow([105,-10,0],[105,55,0], 'EdgeColor','r','FaceColor','r','linewidth',1);
   text(13,103,'i_S','fontsize',14);
   arrow([5,95,0],[20,95,0], 'EdgeColor','b','FaceColor','b','linewidth',2);
   text(13,10,'i_C','fontsize',14);
   arrow([25,15,0],[25,-2,0], 'EdgeColor','k','FaceColor','k','linewidth',2);
   text(70,10,'i_L','fontsize',14);
   arrow([65,15,0],[65,-2,0], 'EdgeColor','m','FaceColor','m','linewidth',2);
   text(93,10,'i_O','fontsize',14);
   arrow([90,15,0],[90,-2,0], 'EdgeColor','r','FaceColor','r','linewidth',2);
   
   
   axis off
   grid off
   axis square
   set(gca,'xLim',[-10 110]);
   set(gca,'yLim',[-10 110]);
  