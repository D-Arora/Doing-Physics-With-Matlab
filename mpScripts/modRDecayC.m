% modRDecayC.m

% RADIOACTIVRE DECAY
% Using ode45 to solve the ODE for radioactive decay: nuclear decay chain   

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/modRDecay.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ON LINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Ian Cooper  matlabvisualphysics@gmail.com
% 191008   Matlab 2019b

tic



close all
clc
clear all

global K

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


% Half-life  tH / Initial number of radioactive nuclei N0
   tH(1) = 2.5e5;
   N0(1) = 10^23;
   tH(2) = 8.04e4;
   N0(2) = 0;
   tH(3) = 1.63e3;
   N0(3) = 0;
   
 % Maximum simulation time
  tMax = 2.5e6;
  
  
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
  N3 = u(:,3);
  

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
  yP = N3;
  plot(xP,yP,'m','linewidth',LW)
  
  grid on
% set(gca,'ytick',0:12.5:100);
  ytickformat('%2.1f')
  
  ylabel('N')
  xlabel('t   [ years ]')
  
  tm1 = '       t_{1/2} = ';  tm2 = num2str(tH,'%2.2e   ');
  tm = [tm1 tm2];
  title(tm)
  
  legend('U-234','Th-230','Ra-226','location','north','orientation','horizontal')
  set(gca,'fontsize',12)  


figure(2)
  pos = [0.35 0.06 0.25 0.55];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
 
subplot(3,1,1)
  xP = t; yP = N1;
  plot(xP,yP,'b','linewidth',LW)
  grid on
  ylabel('N')
  tm1 = 'U-234    t_{1/2} = ';  tm2 = num2str(tH(1),'%2.2e   ');
  tm = [tm1 tm2];
  title(tm)
  set(gca,'fontsize',12)  
  
subplot(3,1,2)
  xP = t; yP = N2;
  plot(xP,yP,'r','linewidth',LW)
  grid on
  ylabel('N')
  tm1 = 'Th-230    t_{1/2} = ';  tm2 = num2str(tH(2),'%2.2e   ');
  tm = [tm1 tm2];
  title(tm)
  set(gca,'fontsize',12)    
 
 subplot(3,1,3)
  xP = t; yP = N3;
  plot(xP,yP,'m','linewidth',LW)
  grid on
  ylabel('N')
  tm1 = 'Ra-226    t_{1/2} = ';  tm2 = num2str(tH(3),'%2.2e   ');
  tm = [tm1 tm2];
  title(tm)
  set(gca,'fontsize',12)  
  xlabel('t  [ years ]');
toc  
  

% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global K
  uDot = zeros(3,1);
  uDot(1) = -K(1)*u(1);
  uDot(2) = -K(2)*u(2) + K(1)*u(1);
  uDot(3) = -K(3)*u(3) + K(2)*u(2);
end  
  



