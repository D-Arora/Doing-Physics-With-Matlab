% np002.m

% Physics of Neurones
% Simple model of membrane voltage% Leak + instataneous I-{Na,p} model


% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191021

tic


%close all
clc
clear

global C GL GNa EL ENa Vh k Iext


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Initial membrane voltage  [0 V]     unstable equilibrium  6.672903e-3
   V0 =  100e-3;  %6.672903e-3; 

% External current input  [ 0.60e-3 A]
   Iext = 0.9e-3;
  
  
% conductance [19e-3  74e-3  S]
  GL = 19e-3; GNa = 74e-3;
% Membrane capacitance [10e-6]
  C = 10e-6;
% Reverse potential / Nerest potential  [EL = -67e-3 V ENa = 60e-3]
  EL = -67e-3; ENa = 60e-3;
% Simulation time  [5e-3 s]
  tMax = 10e-3; 

  

  % V1/2 [V]   k [V]
    Vh = 19e-3; k = 9e-3;

  
  
% SETUP ===============================================================

% Time interval for simulation
% tSpan = [0 tMax];
  nT = 999;
  tSpan = linspace(0,tMax,nT);
  dt = tSpan(2)-tSpan(1);
% Initial conditions
  u0 = V0;
% Options
  options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Solve ODE 
  [t,u] = ode45(@EqM, tSpan, u0, options);
  
% Membrabe voltage [V] 
  V =  u(:,1);  

  m_inf = 1./( 1 + exp((Vh - V)/k) ) ;

% Currents  [A]
   IL   = GL*(V - EL);
   INa  = GNa.*m_inf.*(V - ENa);
   Inet = IL + INa;

   dvdt = gradient(V,dt);   


%%

% v = linspace(-100e-3,100e-3,199);
% mINF = 1./( 1 + exp((Vh - v)/k) ) ;
% figure(99)
% plot(v.*1e3,mINF,'b','linewidth',2)
% hold on
% plot([-100,Vh*1e3],[0.5 0.5],'b')
% plot([Vh*1e3 Vh*1e3],[0 0.5],'b')
% grid on 
% box on
% set(gca,'fontsize',12)
% xlabel('V_M   [ mV ]');
% ylabel('m_{inf}')
% tm = num2str(Vh*1e3,'%3.2f');
% text(Vh*1e3+5,0.05,tm,'fontsize',12)


%%

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
  box on  

  grid on

  ylabel('V  [ mV ]')
  xlabel('t  [ ms ]')
 % tm = sprintf('g = %2.0f mS  C = %2.0f \\muF  E = %2.0f mV  \\tau = %2.2f ms' ...
 %     ,g*1e3,C*1e6,E*1e3, 1e3*tau);
 % title(tm,'fontweight','normal')
  set(gca,'fontsize',12)   

figure(2)
  pos = [0.05 0.10 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  hold on
  
  plot(V.*1e3,Iext.*1e3.*ones(nT,1),'m','linewidth',1)
  Hplot = plot(V(end).*1e3,Iext(end).*1e3,'mo');
  set(Hplot,'markersize',8,'markerfacecolor','m');
  
  plot(V.*1e3,Inet.*1e3,'k','linewidth',1)
  Hplot = plot(V(end).*1e3,Inet(end).*1e3,'ko');
  set(Hplot,'markersize',8,'markerfacecolor','k');
  
  plot(V.*1e3,IL.*1e3,'b','linewidth',2)
  hold on
  Hplot = plot(V(end).*1e3,IL(end).*1e3,'bo');
  set(Hplot,'markersize',8,'markerfacecolor','b');
  
  plot(V.*1e3,INa.*1e3,'r','linewidth',2) ;
  Hplot = plot(V(end).*1e3,INa(end).*1e3,'ro');
  set(Hplot,'markersize',8,'markerfacecolor','r');
  
  tm = 'I_{L} (b)    I_{Na} (r)    I_{ext} (m)    I_{net} (k)';
  title(tm,'fontweight','normal')
  
  grid on
  box on
  
  ylabel('I  [ mA ]')
  xlabel('V  [ mV ]')
  set(gca,'fontsize',12)
%   
figure(3)
  pos = [0.35 0.10 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  hold on
  
  plot(V.*1e3,dvdt,'b','linewidth',LW)
  %plot(V(end).*1e3,IR(end).*1e3,'bo')
  grid on
 
  ylabel('dv/dt  [ V/s ]')
  xlabel('V  [ mV ]')
  set(gca,'fontsize',12)
  box on
  
figure(4)
  pos = [0.35 0.56 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  hold on
  
  plot(t.*1e3,Iext.*1e3.*ones(nT,1),'m','linewidth',1)
  plot(t.*1e3,Inet.*1e3,'k','linewidth',1)
  plot(t.*1e3,IL.*1e3,'b','linewidth',2)
  plot(t.*1e3,INa.*1e3,'r','linewidth',2) ;
  
  
  Hplot = plot(t(end).*1e3,Iext(end).*1e3,'mo');
  set(Hplot,'markersize',8,'markerfacecolor','m');
  
  Hplot = plot(t(end).*1e3,Inet(end).*1e3,'ko');
  set(Hplot,'markersize',8,'markerfacecolor','k');
  
  Hplot = plot(t(end).*1e3,IL(end).*1e3,'bo');
  set(Hplot,'markersize',8,'markerfacecolor','b');
  
  Hplot = plot(t(end).*1e3,INa(end).*1e3,'ro');
  set(Hplot,'markersize',8,'markerfacecolor','r');
    
  tm = 'I_{L} (b)    I_{Na} (r)    I_{ext} (m)    I_{net} (k)';
  title(tm,'fontweight','normal')
  
  grid on
  box on
  
  ylabel('I  [ mA ]')
  xlabel('t  [ ms ]')
  set(gca,'fontsize',12) 
  
  toc
  
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global C gL gNa EL ENa Vh k Iext
  uDot = zeros(1,1);
  m = 1/( 1 + exp((Vh - u(1))/k) );
  uDot(1) = ( Iext - gL*(u(1) - EL) - gNa*m*(u(1) - ENa) )/C;
end  
  



