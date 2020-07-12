% np002B.m

% Physics of Neurones
% [1D] Dynamic System    S.I. units used in Inputs and Calculations
%      Simple model of membrane potential and membrane currents
%      Leak current IL / fast sodium ion current INa model
%      For Leak current IL model, set GNa = 0
% State variable  membrane potential  x   V   VM
%      Comment/uncomment close all statement
%         for plot multiple plots in a Figure Window
% CELL #1 Input parameters
% CELL #2 Phase Portrait Plot  dV/dt vs V
% CELL #3 Time evoluation of the sate variable 
% CELL #4 Functions

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/np001.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191105    Matlab 2019b


tic

%close all
clc
clear

global C GL GNa EL ENa Vh k Iext 


% =====================================================================
% #1   INPUTS
% =====================================================================
% Initial membrane voltage  [0 V]     unstable equilibrium  6.672903e-3
   V0 =  -100e-3;  %6.672903e-3; 

% External current input  [ 0.60e-3 A]
   Iext = 1.00e-3;

% Simulation time  [2e-3 s]
  tMax = 20e-3; 
   
% Conductance [19e-3  74e-3  S]
%   Leak current model  GL = 0;
  GL = 19e-3; GNa = 74e-3;
%  GNa = 0;
  
% Membrane capacitance [10e-6]
  C = 10e-6;
  
% Reverse potential / Nerest potential  [EL = -67e-3 V ENa = 60e-3]
  EL = -67e-3; ENa = 60e-3;

% V1/2 [V]   k [V]
    Vh = 19e-3; k = 9e-3;

 
%% =====================================================================
%  #2   PHASE PORTRAIT PLOT xDot vs x 
%  =====================================================================
% State Variable x    [V]
    xMin = -60e-3;
    xMax = 60e-3;
    Nx = 9999;
    x = linspace(xMin,xMax, Nx);

% Define function
    mINF = 1./( 1 + exp((Vh - x)/k) ) ;
    xDot = ( Iext - GL.*(x-EL) - (GNa*mINF).*(x - ENa) )./C;
   
% Find eqilibrium points xDot = 0  
%   zCount     Number of equilibrium points  
%   zIndex     Indices for equilibrium points 
%   ep = -1    Stable equilibrium point
%   ep = + 1   Unstable equilibrium point

  zCount = 0;
  zIndex = -1; zLen = -1;
 
  for c = 1:Nx-1
    zSign = xDot(c)*xDot(c+1);
      if zSign <= 0
         zCount = zCount + 1;
         zIndex(zCount) = c;
      end
  end
  
  if zIndex > 0
  zLen = length(zIndex);
  if zLen > 0
  % Eigenvalues
    disp('Eigenvalues')
    disp(x(zIndex).*1e3);
    
    ep = zeros(zLen,1);
      for c =  1 : zLen
        ep(c) = sign( xDot(zIndex(c)+1) - xDot(zIndex(c)-1) );
      end
  end
  end
figure(1)  % Phase Portrait Plot
  pos = [0.05 0.55 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
 
  hold on
  xP = x.*1e3; yP = xDot;
  plot(xP,yP,'b','linewidth',LW)
  
for c = 1 : zLen
  if ep(c) == -1; col = 'r'; end
  if ep(c) ==  1; col = 'g'; end
  Hplot = plot(xP(zIndex(c)),yP(zIndex(c)),'o');
  set(Hplot,'markersize',8,'markerfacecolor',col,'markeredgecolor',col);
end

  box on
  grid on
  xlabel('V_m  [ mV ]  ')
  ylabel('dV_M / dt   [ V/s ]')
  set(gca,'fontsize',14)
  
figure(2)
  pos = [0.05 0.05 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;

  plot(x.*1e3,mINF,'b','linewidth',2)
  hold on
  plot([-100,Vh*1e3],[0.5 0.5],'b')
  plot([Vh*1e3 Vh*1e3],[0 0.5],'b')
  grid on 
  box on
  set(gca,'fontsize',12)
  xlabel('V_M   [ mV ]');
  ylabel('m_{inf}')
  tm = num2str(Vh*1e3,'%3.2f');
  text(Vh*1e3+5,0.05,tm,'fontsize',12)  
  

%% =====================================================================
%  #3   TIME EVOLUTION OF STATE VARIABLE:  membrane potential
%  =====================================================================  

% Time interval for simulation
% tSpan = [0 tMax];
  nT = 199;
  tSpan = linspace(0,tMax,nT);
  dt = tSpan(2)-tSpan(1);
% Initial conditions
  u0 = V0;
% Options
%  options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Solve ODE 
  [t,u] = ode45(@EqM, tSpan, u0);
  
% Membrane voltage [V] 
  V =  u(:,1);  

  m_inf = 1./( 1 + exp((Vh - V)/k) ) ;

% Currents  [A]
   IL   = GL*(V - EL);
   INa  = GNa.*m_inf.*(V - ENa);
   Inet = IL + INa;


figure(3)   % t  vs  Vm
  pos = [0.35 0.56 0.25 0.35];
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

figure(4)   % Vm  vs  I
  pos = [0.35 0.05 0.25 0.35];
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
   
figure(5)   % t  vs  I
  pos = [0.65 0.56 0.25 0.35];
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

%% =====================================================================  
% #4    FUNCTIONS    
% =====================================================================  
function uDot = EqM(t,u)
  global C GL GNa EL ENa Vh k Iext 
  uDot = zeros(1,1);
  
  if t > 0.001
   Iext = 1e-3; 
  end
  
  m = 1/( 1 + exp((Vh - u(1))/k) );
  uDot(1) = ( Iext - GL*(u(1) - EL) - GNa*m*(u(1) - ENa) )/C;
  
  
end  
  



