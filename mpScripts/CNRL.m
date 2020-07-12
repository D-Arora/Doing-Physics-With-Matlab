%CNRL.m

% RL CIRCUIT
% Finite Difference Method for RL circuit analysis

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
% 180212

clear all
close all
clc

% =======================================================================
% INPUTS: S.I. UNITS  default values [  ]
% =======================================================================
% Resistance [1e3]
    R  = 1.0e3;    
% Inductance [1e-3];      
    L  = 10e-3;   
% Peak source emf [10] 
  VS = 10;         
% Source emf: flagV = 1  step function off/on
%             flagV = 2  step function on/off
%             flagV = 3  pulses (square wave)
%             flagV = 4  sinusoidal emf
%             flagV = 5  rectified sinusoidal emf
%             flagV = 6  superposition of sinusoidal emfs
  flagV = 6;
   
   
%  ======================================================================   
%  Define emf
% =======================================================================
% Time constant and time step
   tau = L/R;
   dt = tau/1000;
   
   switch flagV
    case 1  % step function  off to on
       % max time interval  [10*tau]
         tMax = 10*tau;
       %  tMax = 500e-6;        
       % percentage off time  [10]
         pOFF = 10;         
                
         t = 0:dt:tMax;
         N = length(t);
         nOFF = round(N*pOFF/100);
         vS = zeros(1,N);
         vS(nOFF:end) = VS;
     
    case 2  % step function  on to off
        % max time interval  [10e-3] 
         tMax = 10e-6;        
       % percentage on time  [50]
         pON = 10;         
       
         t = 0:dt:tMax;
         N = length(t);
         nON = round(N*pON/100);
         vS = zeros(1,N);
         vS(1:nON) = VS;
         
    case 3    % pulses ON / OFF
       % max time interval  [10e-3] 
         tMax = 100e-6;        
       % percentage off time  [10]
         pOFF = 50;         
       % Number of periods     [3]
         nT = 4;
       
       % Calculations  
         t = 0:dt:tMax;
         N = length(t);  
                  numT = round(N/nT-0.5);
         nON = round(numT*pOFF/100);
         nOFF = numT-nON;
         vS = zeros(1,N);
         for c = 1 : nT
           vS(nOFF + numT*(c-1) : nOFF + numT*(c-1)+nON) = VS;
         end

    case 4     % sinusoidal input    [fS = 1e4]
       fS = 10e3;
       w = 2*pi*fS;
       T = 1/fS;
       nT = 8; 
       tMax = nT*T;
       t = 0:dt:tMax;
       N = length(t);
       vS = VS .* exp(1j*(w*t - pi/2));
    
       
     case 5     % sinusoidal input - rectified
       f = 10e3;
       w = 2*pi*f;
       T = 1/f;
       nT = 8;
       tMax = nT*T;
       t = 0:dt:tMax;
       N = length(t);
       vS = abs(VS .* sin(w*t));
            
     case 6     % Input: superposition of sinusoidal signals
       f = 10e3;
       w = 2*pi*f;
       T = 1/f;
       nT = 8;
       tMax = nT*T;
       t = 0:dt:tMax;
       N = length(t);
       vS = VS .* sin(w*t) + (0.5*VS) .* sin(10*w*t);

   end
  
   
% ======================================================================  
% CALCULATIONS 
% ======================================================================
   k = dt/L;
%  Initialise values 
  vR = zeros(1,N);
  vL = zeros(1,N);
  iS = zeros(1,N);
    
%  Time Step #1
  if flagV == 2
     iS(1) = vS(1)/R;
     vR(1) = vS(1);
  end
  
% Time Steps #2 to #N
  for c = 2 : N
    iS(c) = iS(c-1) + k*vL(c-1);
    vR(c) = iS(c) * R;
    vL(c) = vS(c) - vR(c);
  end
 

% Powers and energy
     pS = real(vS) .* real(iS);
     pR = real(vR) .* real(iS);
     pL = real(vL) .* real(iS);
     
      uS = zeros(1,N); uR = zeros(1,N); uL = zeros(1,N);
     if flagV == 2
         uS(1) = pS(1)*dt;
         uR(1) = uS(1);
     end
     
    
     for c = 2 : N
         uS(c) = uS(c-1) + pS(c)*dt;
         uR(c) = uR(c-1) + pR(c)*dt;
         uL(c) = uL(c-1) + pL(c)*dt;
     end

     
% Phases  [degrees]
% INPUT time to calculate phasors
     tP = 0.75*tMax;
     
     nP = find(t > tP,1);
     phi_vS   = rad2deg(angle(vS(nP)));
     phi_vR   = rad2deg(angle(vR(nP)));
     phi_vL   = rad2deg(angle(vL(nP)));
     theta_iS = rad2deg(angle(iS(nP)));
     
     if flagV == 4
        [iS_Peaks, t_Peaks] = findpeaks(real(iS));
       nPeaks = length(t_Peaks)-1;
       T_Peaks = (t(t_Peaks(end)) - t(t_Peaks(end-nPeaks)))/(nPeaks);
       f_Peaks = 1 / T_Peaks;
       
       ZR = R;
       ZL = 1j * w*L;
       Z = ZR + ZL;
       Zmag = abs(Z);
       phi_ZR = rad2deg(angle(ZR));
       phi_ZL = rad2deg(angle(ZL));
       phi_Z = rad2deg(angle(Z));
     end


% ======================================================================
% GRAPHICS
% ======================================================================

figure(1)    
   FS = 14;
   pos = [0.07 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
subplot(3,1,1)   % -----------------------------------------------------
     xP = 1e3.*t; yP = real(vS);
     plot(xP,yP,'b','linewidth',2');
     grid on
     set(gca,'xLim',[0 1e3*max(t)]);
  %  set(gca,'yTick',-10:5:10);
     set(gca,'yLim',[-VS-0.2 VS+0.2]);
     tm1 = '   R = ';
     tm2 = num2str(R/1e3, '%3.2f');
     tm3 = ' k\Omega   L = ';
     tm4 = num2str(L*1e3, '%3.2f');
     tm5 = ' mH   \tau = ';
     tm6 = num2str(tau*1e3, '%3.4f');
     tm7 = ' ms';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',10);
     ylabel('v_S  [ V ]');
     set(gca,'fontsize',FS);
     
subplot(3,1,2) % ---------------------------------------------------------
     xP = 1e3.*t; yP = real(vR);
     plot(xP,yP,'r','linewidth',2');
     grid on
     set(gca,'xLim',[0 1e3*max(t)]); 
     set(gca,'yLim',[-VS-0.2 VS+0.2]);
%    set(gca,'yTick',-10:5:10);
     ylabel('v_R   [ V ]');
     set(gca,'fontsize',FS);
     
subplot(3,1,3) % ---------------------------------------------------------
     xP = 1e3.*t; yP = real(vL);
     plot(xP,yP,'m','linewidth',2');
     grid on
     set(gca,'yLim',[-VS-0.2 VS+0.2]);
     set(gca,'xLim',[0 1e3*max(t)]); 
     xlabel('t  [ ms ]');
     ylabel('v_L  [ V ]');
     set(gca,'fontsize',FS);

figure(2) % =============================================================    
   FS = 14;
   pos = [0.37 0.05 0.28 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   xP = 1e3.*t; yP = real(iS.*1e3);
   plot(xP,yP,'r','linewidth',2');
   grid on
  
   tm1 = '   R = ';
   tm2 = num2str(R/1e3, '%3.2f');
   tm3 = ' k\Omega   L = ';
   tm4 = num2str(L*1e3, '%3.2f');
   tm5 = ' mH   \tau = ';
   tm6 = num2str(tau*1e3, '%3.4f');
   tm7 = ' ms';
   tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
   title(tm,'fontweight','normal','fontsize',10);
   xlabel('t  [ ms ]');
   ylabel('i_S = i_R = i_L [ mA ]');
   set(gca,'fontsize',FS)

figure(3)   % ==========================================================   
   FS = 12;
   pos = [0.67 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
 subplot(3,1,1) % ------------------------------------------------------
     xP = 1e3.*t; yP = pS;
     plot(xP,yP,'b','linewidth',2');
     grid on
     ylabel('p_S  [ W ]');
     tm1 = '   R = ';
     tm2 = num2str(R/1e3, '%3.2f');
     tm3 = ' k\Omega   L = ';
     tm4 = num2str(L*1e3, '%3.2f');
     tm5 = ' mH   \tau = ';
     tm6 = num2str(tau*1e3, '%3.4f');
     tm7 = ' ms';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',10);
     set(gca,'fontsize',FS);
     
subplot(3,1,2) % ------------------------------------------------------
     xP = 1e3.*t; yP = pR;
     plot(xP,yP,'r','linewidth',2');
     grid on
     ylabel('p_R  [ W ]');
     set(gca,'fontsize',FS);

 subplot(3,1,3) % -----------------------------------------------------
     xP = 1e3.*t; yP = pL;
     plot(xP,yP,'m','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('p_L  [ W ]');
     set(gca,'fontsize',FS);

figure(4)  % ============================================================ 
   FS = 12;
   pos = [0.47 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
subplot(3,1,1) % -------------------------------------------------------
     xP = 1e3.*t; yP = uS;
     plot(xP,yP,'b','linewidth',2');
     grid on
     tm1 = '      R = ';
     tm2 = num2str(R/1e3, '%3.2f');
     tm3 = ' k\Omega   L = ';
     tm4 = num2str(L*1e3, '%3.2f');
     tm5 = ' mH   \tau = ';
     tm6 = num2str(tau*1e3, '%3.4f');
     tm7 = ' ms';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     title(tm,'fontweight','normal','fontsize',10);
     ylabel('u_S  [ J ]');
     set(gca,'fontsize',FS);
     set(gca,'yLim',[0 1e-6]);
     
subplot(3,1,2) % --------------------------------------------------
     xP = 1e3.*t; yP = uR;
     plot(xP,yP,'r','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('u_R  [ J ]');
     set(gca,'fontsize',FS);
  %   set(gca,'yLim',[0 1e-5]);
     
subplot(3,1,3)  % ------------------------------------------------------
     xP = 1e3.*t; yP = uL;
     plot(xP,yP,'m','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('u_L  [ J ]');
     set(gca,'fontsize',FS);
    % set(gca,'yLim',[0 1e-6]);


 figure(5)  % TEXT------------------------------------------------ 
   %FS = 14;
   pos = [0.27 0.05 0.35 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   x = 0; y = 100; dy = -15;
   tm1 = 'R  = ';
   tm2 = num2str(R/1e3,'%3.2f');
   tm3 = '  k\Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = 'L = ';
   tm2 = num2str(L*1e3,'%3.2f');
   tm3 = '  mH';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = '\tau = ';
   tm2 = num2str(tau*1e3,'%3.2e');
   tm3 = '  ms';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
 if flagV == 4
   y = y+dy;
   tm1 = 'f_S = ';
   tm2 = num2str(fS/1e3,'%3.2e');
   tm3 = '  kHz';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14); 
   
   y = y+dy;
   tm1 = 'f_{Peaks} = ';
   tm2 = num2str(f_Peaks/1e3,'%3.2e');
   tm3 = '  kHz';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14); 
   
   x = 40; y = 100-dy;  
   y = y+dy;
   tm1 = 'Z_R = ';
   tm2 = num2str(R,'%3.2f');
   tm3 = '  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
     
   y = y+dy;
   tm1 = 'Z_L = ( ';
   tm2 = num2str(real(ZL),'%3.2f');
   tm3 = ' + ';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = ' ';
   tm2 = num2str(imag(ZL),'%3.2f');
   tm3 = ' j )  \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x+10,y+8,tm);
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
   h_text = text(x+10,y+8,tm);
   set(h_text,'fontsize',14); 
   
   y = y+dy+10;
   tm1 = '|Z| = ';
   tm2 = num2str(abs(Z),'%3.2f');
   tm3 = ' \Omega';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = '\phi_R = ';
   tm2 = num2str(phi_ZR,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y+dy;
   tm1 = '\phi_L = ';
   tm2 = num2str(phi_ZL,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);

   y = y+dy;
   tm1 = '\phi_Z = ';
   tm2 = num2str(phi_Z,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   x = 78; y = 100;
   tm1 = 'At time t = ';
   tm2 = num2str(t(nP)*1e3,'%3.3f');
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
   tm1 = '\phi_R = ';
   tm2 = num2str(phi_vR,'%3.1f');
   tm3 = '^o';
   tm = [tm1 tm2 tm3];
   h_text = text(x,y,tm);
   set(h_text,'fontsize',14);
   
   y = y + dy;
   tm1 = '\phi_L = ';
   tm2 = num2str(phi_vL,'%3.1f');
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
 end
   axis off
   set(gca,'xLim',[0 100]);
   set(gca,'yLim',[0 100]);
   set(gca,'Position', [0.05 0.1 0.8 0.85]);
   
if flagV == 4
figure(6)   % phasors ------------------------------------------------
   xyLim = 12;
   xyTick = -xyLim:4:xyLim;
   
   FS = 14;
   pos = [0.50 0.15 0.28 0.43];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   %xP = [0 real(vR(nP))]; yP = [0 imag(vR(nP))];
   %plot(xP,yP,'r','linewidth',2');
  
   hold on
   arrow([0,0],[real(vR(nP)),imag(vR(nP)),0],'BaseAngle',35, 'EdgeColor','r',...
       'FaceColor','r',...
      'Length',350, 'linewidth',2);
   %x1 = [0 real(vL(nP))]; yP = [0 imag(vL(nP))];
  % plot(xP,yP,'m','linewidth',2');
   arrow([0,0],[real(vL(nP)),imag(vL(nP)),0],'BaseAngle',30, 'EdgeColor','m', ...
       'Length',50,'FaceColor','m','linewidth',2);
   
   arrow([0,0],[real(vS(nP)),imag(vS(nP)),0],'BaseAngle',30, 'EdgeColor','b',...
       'FaceColor','b',...
      'Length',30, 'linewidth',2);
   
   % xP = [0 real(vS(nP))]; yP = [0 imag(vS(nP))];
  % plot(xP,yP,'b','linewidth',2');
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
   box on
   legend(' v_R',' v_L',' v_S');
  
     tm1 = 'Phasors at time t  =   ';
     tm2 = num2str(1e3*t(nP), '%3.2f');
     tm3 = ' ms';
     tm = [tm1 tm2 tm3];
     title(tm,'fontweight','normal','fontsize',12);  
end

    
  