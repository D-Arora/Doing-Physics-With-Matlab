% modRDecayA.m

% RADIOACTIVRE DECAY
% Using ode45 to solve the ODE for radioactive decay    

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
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
clear

global K


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Default values
%   tH = 100 N0 = 100 tMax = 10*tH

% Half-life  tH / Initial number of radioactive nuclei N0
   tH = 25;
   N0 = 100;

 % Maximum simulation time
  tMax = 8*tH;
  
  
% SETUP ===============================================================

% Decay constant lambda  K
  K = log(2)/tH;
% Time interval for simulation
  tSpan = [0 tMax];
% Initial conditions
  u0 = N0;
% Options
  options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Solve ODE 
  [t,u] = ode45(@EqM, tSpan, u0, options);
  
% Undecayed nuclei as a function of time 
  N =  u(:,1);  

% Analytical solution
  NA = N0 .* exp(-K*t);
  

% GRAPHICS ============================================================  

figure(1)
  pos = [0.05 0.06 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
 
  xP = t; yP = N;
  plot(xP,yP,'b','linewidth',LW)
  
  hold on
  xP = t(1:10:end); yP = NA(1:10:end);
  plot(xP,yP,'ro','linewidth',LW)
  
  grid on
  set(gca,'xtick',0:25:200);
  set(gca,'ytick',0:12.5:100);
  ytickformat('%2.1f')
  
  ylabel('N')
  xlabel('t  [s]')
  
  tm1 = '  t_{1/2} = ';  tm2 = num2str(tH,'%2.2f');
  tm = [tm1 tm2];
  title(tm)
  
  set(gca,'fontsize',12)  

  
toc  

  
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global K
  uDot = zeros(1,1);
  uDot(1) = -K*u(1);
end  
  



