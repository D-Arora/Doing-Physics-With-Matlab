% np001A.m

% Physics of Neurones
% Simple model of membrane voltage: I_leak / Na+ constant conductance model

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/np001.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191025

tic

close all
clc
clear

global  C GL GNa EL ENa


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% conductance [19e-3   74e-3  S]
  GL = 19e-3;  GNa = 74e-3;
% Membrane capacitance [10e-6]
  C = 10e-6;
% Reverse potential / Nernst potential  [EL = -67e-3 V   ENa = 60e-3)
  EL = -67e-3;   ENa = 60e-3;
% Simulation time  [13-3 s]
  tMax = 1e-3;
% Initial membrane voltage  [0 V]
  V0 = 100e-3;
  
  
% SETUP ===============================================================

% Time interval for simulation
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

% Currents  [A]
  IL  = GL.*(V - EL);
  INa = GNa.*(V - ENa);
  dVdt = gradient(V,dt);
  IC = C.*dVdt; 
  Inet = IL + INa + IC;
  
  
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
  grid on
  ylabel('V_M  [ mV ]')
  xlabel('t  [ ms ]')
  set(gca,'fontsize',12)   
  box on
  
figure(2)
  pos = [0.05 0.10 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  plot(t,IL.*1e3,'b','linewidth',LW)
  hold on
  plot(t,INa.*1e3,'r','linewidth',LW)
  plot(t,IC.*1e3,'m','linewidth',LW)
  plot(t,Inet,'k','linewidth',LW)

  grid on
  legend('I_L','I_{Na}','I_C','I_{net}','location','south', ...
      'orientation','horizontal')
  ylabel('I  [ mA ]')
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
  plot(V.*1e3,INa.*1e3,'r','linewidth',LW)
  
  Hplot = plot(V(end).*1e3,IL(end).*1e3,'bo');
  set(Hplot,'markersize',8,'markerfacecolor','b');
  Hplot = plot(V(end).*1e3,INa(end).*1e3,'ro');
  set(Hplot,'markersize',8,'markerfacecolor','r');
  
  grid on
 
  legend('I_L','I_{Na}','location','north','orientation','horizontal')
  ylabel('dV_M / dt  [ V/s ]')
  xlabel('V_M  [ mV ]')
  set(gca,'fontsize',12)
  box on 
  
  toc
  
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global  C GL GNa EL ENa
  
  uDot = zeros(1,1);
  uDot(1) = ( -GL*(u(1) - EL) - GNa*(u(1) - ENa) )/C;
end  
  



