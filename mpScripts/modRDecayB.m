% modRDecayB.m

% RADIOACTIVRE DECAY
% Using ode45 to solve the ODE for radioactive decay    

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/modRDecay.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ON LINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Ian Cooper  MatlabVisualPhysics@gmail.com

% 190925   Matlab 2019b

tic



close all
clc
clear all

global K

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Default values
%   tH1 = 20 N01 = 100 th2 = 30 N02 = 0 tMax = 10*tH

% Half-life  tH / Initial number of radioactive nuclei N0
   tH(1) = 20;
   N0(1) = 100;
   tH(2) = 300;
   N0(2) = 0;
   
 % Maximum simulation time
  tMax = 20*tH(1);
  
  
% SETUP ===============================================================

% Decay constant lambda  K
  K = log(2)./tH;
% Time interval for simulation
  tSpan = [0 tMax];
% Initial conditions
  u0 = N0;
% Options
  options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Solve ODE 
  [t,u] = ode45(@EqM, tSpan, u0, options);
  
% Undecayed nuclei as a function of time 
  N1 = u(:,1);  
  N2 = u(:,2);

  
% GRAPHICS ============================================================  

figure(1)
  pos = [0.05 0.06 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
 
  xP = t; yP = N1;
  plot(xP,yP,'b','linewidth',LW)
  hold on
  yP = N2;
  plot(xP,yP,'r','linewidth',LW)
  
  grid on
  set(gca,'ytick',0:12.5:100);
  ytickformat('%2.1f')
  
  ylabel('N')
  xlabel('t')
  
  tm1 = '  t_{1/2} = ';  tm2 = num2str(tH,'%2.2f   ');
  tm = [tm1 tm2];
  title(tm)
  
  set(gca,'fontsize',12)  

  
toc  

  
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global K
  uDot = zeros(2,1);
  uDot(1) = -K(1)*u(1);
  uDot(2) = -K(2)*u(2) + K(1)*u(1);
end  
  



