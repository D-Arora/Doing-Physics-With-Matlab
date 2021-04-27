% Laplace05C.m

% Ian Cooper
% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% SCRIPTS
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% SYSTEM: mass m  / spring k / dashpot-damping b.
% System driven through the dashpot by a sinusoidal displacement y of the
%      piston inside the dashpot.
% Laplace Transform (LT) used to solve the ODE describing the System
%   to give the displacement x(t)of a mass excited by an input signal f(t).

%   m*x'' + b*x' + k*x = f

% 1   Enter system parameters and 9initial conditions
% 2   Define the input signal (driving force) f(t).
% 3   LT of input signal --> F(s).
% 4   LT of ODE -->
%       m*(s^2*X - s*x(0) - x'(0)) + b*(s*X - x(0)) + k*X = F
%      (m*s^2 + b*s + k)*X - m*(s*x(0) - x'(0)) - b*x(0) - F = 0
%  5  Solve equation for X(s)
%  6  Inverse LT --> x(t)
%  7  Miscellaneous calculations: Bode plots and Nyquist plot

% The symbolic solution for the output x(t) is displayed in the
%   Command Window and a summary of results displayed in Figure(2).

% Numeric variables:
%    time tN / displacement xN / velocity vN / acceleration aN
%    input - driving force fN.

% The Script may need to be amended for different input functions

% 20210307     Matabl R3030b

clear
close all
clc

syms s t X

% SYSTEM PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% mass m / damping constant b / spring constant k
  m = 1;
  b = 1;
  k = 2;
  
% initial conditions: displacement x(0) and velocity v(0)  
  x0 = 1;
  v0 = 0.2;
  
% time scales: number of time steps, max simulation time [t/pi], / time 
  num = 501;
  tMax = 6;
  tN = linspace(0,tMax*pi, num);

  
% INPUT FUNCTION  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% driving frequency
  w = 1;
  
% amplitude
  A = 1;
  
% displacement of dashpot (numerical)
  yN = A.*cos(w*tN);
  
% Input signal (numeric)
  fN = -(b*A*w).*sin(w*tN);
     
% Input signal (symbolic)
  f = -(b*A*w)*sin(w*t);

% Text for input signal 
 txt_Input = 'f = -(b*A*w)*cos(w*t)'; 

% =================================================================== 
% LAPLACE TRANSFORMS AND SOLVING FOR OUTPUT SIGNAL X(s)  ---> x(t)
% MASS: displacement, velocity and acceleration
% L.T. of input signal  
  F = laplace(f,t,s);
  
% L.T. ODE equation 
  Z = (m*s^2 + b*s + k)*X - m*(s*x0 + v0) - b*x0 - F;

% Solve for X
  Sol_x = solve(Z, X)

% INVERSE L.T. to find x(t)
  sol_x = ilaplace(Sol_x,s,t)
  pretty(sol_x)
  
% Velocity calculation
  Sol_v = s*Sol_x - x0
  sol_v = ilaplace(Sol_v)  
  pretty(sol_v)
   
% Acceleration calculation
  Sol_a = s*s*Sol_x - s*x0 - v0
  sol_a = ilaplace(Sol_a)  
  pretty(sol_a)

% Extract numbers from symbolic expressions for output x v and a
  x = subs(sol_x,{t},{tN});
  xN = double(x);
  
  v = subs(sol_v,{t},{tN});
  vN = double(v);
  
  a = subs(sol_a,{t},{tN});
  aaN = double(a);
  
% Acceleration calculation from ODE  
  aN = -(b/m).*vN - (k/m).*xN + (fN./m) ;
  
  
% CALCULATIONS  ==========================================================
% Natural frequency == resonance frequency
  wN = sqrt(k/m);
% Natural period
  TN = 2*pi/wN;  
% Period of input signal
  T = 2*pi/w;
    
% Bode plots and Nyquist plot
% Input driving frequency
  wD = linspace(0.001,5,500);
% Transfer function
  G = b.*wD.*1i ./ (-m*wD.^2 + k + b.*wD.*1i);
% Ampitude of output signal
  Gmag = abs(G);
% Phase lag of output cf input [pi]    
  phi = angle(G)./pi;

% Transfer function at driving frequency   
  Gw = b.*w.*1i ./ (-m*w.^2 + k + b.*w.*1i);
  Gw_mag = abs(Gw);
  phi_w  = angle(Gw)/pi;

% Phase shift time of input to input
  dt = pi*phi_w/w;  
%   wR = sqrt(k/m - (b/(2*m))^2);
 
% Period of response - from peaks
%   [pks, locs] = findpeaks(xN);
%   if length(locs) > 2
%      P = ( tN(locs(end)) - tN(locs(2)) )/(length(locs)-2);
%   else
%      P = -1;
%   end
  

  
% GRAPHICS  ==============================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.2 0.35 0.45]);
  set(gcf,'color','w');
  FS = 14; LW = 2;
  box on
  
  xP = tN/pi;
  
subplot(2,2,1)  
  plot(xP,xN,'r','linewidth',LW)
  hold on
  plot(xP,yN,'b','linewidth',LW)
  set(gca,'xtick',0:2:tMax)
  xlabel('t / \pi');
  text(-2.7,0.25,'y','fontsize',FS,'Color',[0 0 1])
  text(-2.7,-0.25,'x','fontsize',FS,'Color',[1 0 0])
  ytickformat('%.1f')
  grid on
  txt = sprintf('m = %2.1f   b = %2.1f   k = %2.1f \n', m, b, k);
  title(txt,'fontweight','normal')
  set(gca,'fontsize',FS)
 
    
subplot(2,2,2) 
  yyaxis left
  plot(xP,vN,'r','linewidth',2)
  xlabel('t / \pi'); ylabel('v')
  grid on
  set(gca,'fontsize',FS)
  set(gca,'ycolor',[1 0 0])
  set(gca,'xtick',0:2:tMax)
  txt = sprintf('\\omega_N = %2.2f \n', wN);
  title(txt,'fontweight','normal')
  
  yyaxis right
  plot(xP,aaN,'m','linewidth',1.2)
%  ylim([-1 1])
  xlabel('t / \pi'); ylabel('a')
  set(gca,'ycolor',[1 0 1])
  
subplot(2,2,3)  
  plot(xN,vN,'r','linewidth',LW)
  hold on
  plot(x0,v0,'go','markerfacecolor','g','markersize',6)
  xlabel('x'); ylabel('v')
  grid on
  xlim([-1.1 1.1])
  ylim([-2.1, 2.1])
  set(gca,'xtick',-1:0.5:1)
  set(gca,'ytick',-2:1:2)
  xtickformat('%.1f')
  ytickformat('%.1f')
  set(gca,'fontsize',FS)
    
subplot(2,2,4)  
  plot(xP,fN,'b','linewidth',LW)
  xlabel('t / \pi'); ylabel('f')
  txt = sprintf('\\omega = %2.2f    T = %2.2f  \n', w, T);
  title(txt,'fontweight','normal')
  grid on
  set(gca,'fontsize',FS)  

  
figure(3)  % ========================================================
  set(gcf,'units','normalized');
  set(gcf,'position',[0.5 0.2 0.25 0.55]);
  set(gcf,'color','w');
  FS = 14; LW = 2;
  box on
  
subplot(3,1,1)
  plot(wD,Gmag,'r','linewidth',LW)
  hold on
  plot([w w], [0 Gw_mag],'k','linewidth',1)
  plot([0 w], [Gw_mag Gw_mag],'k','linewidth',1)
  plot(w,Gw_mag,'ro','markerfacecolor','r','markersize',8)
  
  txt = sprintf('\\omega_N = %2.2f   \\omega = %2.2f    |G| = %2.2f  \n',wN, w, Gw_mag);
  title(txt,'fontweight','normal')
  ylim([0 1.2])
  xlabel('\omega_D')
  ylabel('|G|')
  grid on
  set(gca,'fontsize',FS)
  
subplot(3,1,2)
  plot(wD,phi,'b','linewidth',LW)
  hold on
  plot([w w], [-1 phi_w],'k','linewidth',1)
  plot([0 w], [phi_w phi_w],'k','linewidth',1)
  plot(w,phi_w,'ro','markerfacecolor','r','markersize',8)
  
  ylim([-0.6 0.6])
  xlabel('\omega_D')
  ylabel('- \phi / \pi')
  txt = sprintf('\\phi = %2.2f \\pi  \\deltat = %2.2f   \n',-phi_w,dt);
  text(0.0,0.65,txt,'fontsize',14)
  
  if phi_w > 0
      txt = 'OUT leads IN';
      text(2.8,0.80,txt,'fontsize',12)
  end
  
  if phi_w < 0
      txt = 'OUT lags IN';
      text(2.8,0.80,txt,'fontsize',12)
  end
  
  grid on  
  set(gca,'fontsize',FS)
  
subplot(3,1,3)
  plot(real(G),imag(G),'k','linewidth',1)
  hold on
  plot([0 real(Gw)], [0 -imag(Gw)],'r','linewidth',1)
  plot([0 0.5], [0 0],'r','linewidth',1)
  plot(real(Gw),-imag(Gw),'ro','markerfacecolor','r','markersize',8)
  
  plot(real(G),-imag(G),'k','linewidth',1)
   
  xlabel('real(G)')
  ylabel('-imag(G)')
  grid on  
  axis square
  set(gca,'fontsize',FS) 
  
  
  
