% np001.m

% Physics of Neurones
% Simple model of membrane voltage: I_leak model

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/np001.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191105


tic

close all
clc
clear

global tau E


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% conductance [19e-3  S]
  G = 19e-3;
% Membrane capacitance [10e-6]
  C = 10e-6;
% Reverse potential / Nernst potential  [EL = -67e-3 V)
  E = -67e-3;
% Simulation time  [5]
  N = 10;
% Initial membrane voltage  [0 V]
  V0 = -100e-3;
  
  
% SETUP ===============================================================

% Time constant / max simulation time interval   [s]
  tau = C/G;
%  tMax = N * tau;
  tMax = 5e-3;
% Time interval for simulation
% tSpan = [0 tMax];
  tSpan = linspace(0,tMax,999);
  dt = tSpan(2)-tSpan(1);
% Initial conditions
  u0 = V0;
  
% Options
  options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% Solve ODE 
  [t,u] = ode45(@EqM, tSpan, u0, options);
  
% Membrane potential [V] 
  V =  u(:,1);  

% Analytical solution
  VA = E + (V0 - E).*exp(-t/tau);
  
% Currents  [A]
  IL = G.*(V - E);
  dVdt = gradient(V,dt);
  IC = C.*dVdt; 
  Inet = IL+IC;
  
  
% GRAPHICS ============================================================  
 
figure(1)
  pos = [0.05 0.56 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
  hold on
xP = t.*1e3; yP = V.*1e3;
  plot(xP,yP,'b','linewidth',LW)
  
  hold on
xP = t(1:100:end).*1e3; yP = VA(1:100:end).*1e3;
  plot(xP,yP,'ro','linewidth',LW)
  grid on

  ylabel('V_M  [ mV ]')
  xlabel('t  [ ms ]')
  tm = sprintf('G = %2.0f mS  C = %2.0f \\muF  E = %2.0f mV  \\tau = %2.2f ms' ...
      ,G*1e3,C*1e6,E*1e3, 1e3*tau);
  title(tm,'fontweight','normal')
  set(gca,'fontsize',12)  
  box on

figure(2)
  pos = [0.05 0.10 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  plot(t,IL.*1e3,'b','linewidth',LW)
  hold on
%   plot(t,IC.*1e3,'r','linewidth',LW)
%   plot(t,Inet,'k','linewidth',LW)

  grid on
%  legend('I_L','I_C','I_{net}')
  ylabel('I_L  [ mA ]')
  xlabel('t  [ ms ]')
  set(gca,'fontsize',12)
  
figure(3)
  pos = [0.35 0.10 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
    
  plot(V.*1e3,dVdt,'b','linewidth',LW)
  hold on
  Hplot = plot(V(end).*1e3,dVdt(end),'bo');
  set(Hplot,'markersize',8,'markerfacecolor','b');
  grid on
 
  ylabel('dV_M / dt  [ V/s ]')
  xlabel('V_M  [ mV ]')
  set(gca,'fontsize',12)
  box on
  
  
  figure(4)
  pos = [0.35 0.56 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
    
  plot(V.*1e3,IL.*1e3,'b','linewidth',LW)
  hold on
  Hplot = plot(V(end).*1e3,dVdt(end),'bo');
  set(Hplot,'markersize',8,'markerfacecolor','b');
  grid on
 
  ylabel('I_L  [ mA ]')
  xlabel('V_M  [ mV ]')
  set(gca,'fontsize',12)
  box on
  
  toc
  
  
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global tau E
  uDot = zeros(1,1);
  uDot(1) = -(1/tau) *  (u(1) - E);
end  
  



