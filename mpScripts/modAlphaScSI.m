% modAlphaScSI.m

% ALPHA PARTICLE SCATTERING
% COULOMB'S LAW    
%  S.I units
%   Solving the equation of motion using the ode45 function
%  Variation of impact parameter y0

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/modAlphaSc.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ON LINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190919   Matlab 2019b

global G

% close all
clc
clear

tic

% CONSTANTS ===========================================================
% Mass alpha particle  [kg] 
  m = 4*1.7e-27;
% Fundamental charge   [c]
  e = 1.602e-19;
% Coulomb Constant [N.m^2.C^-2]
  k = 9.0e9;
% Alpha particle charge
  q = 2*e;
% Uranium nuclei charge  
  Q = 92*e;
% Distance scaling [m]
  S = 1e-15;
  
  
% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Default values
%   x0 = -5   vx0 = 2   y0 = 1.0   vy0 = 0
%   tMax = 8  nT = 181
%   Larger values of nT give more accurate results

% Matrix u (X and Y displacements and velocities)
%  u(1) = x   u(2) = vx   u(3) = y; u(4) = vy = 0
% Initial conditions: t = 0
% KE alpha particle
  KE = 40*1e6*e;
  
% Impact parameter >>>>>>>>>>>>>>>
   y0  = 4.8*S;
   
   x0  = -20*S;
   vx0 = sqrt(2*KE/m);
   vy0 = 0;
    
 % Max time interval for simulation / number of calculations
 % tMax = abs(2*x0/vx0);
  tMax = 1e-21;
  nT = 15531; 

% Grid dimensions
  Grid = 10;
  xGrid = [-10 6];
  yGrid = [0 16];
  xtickGrid = -10:2:6;
  ytickGrid = 0:2:16;

% Plot figure 2   (1) yes /  (2)  no
  flag2 = 1;
  if flag2 == 1; close all; end
  
% Setup  ==============================================================
  
% Coulomb equation constant k q Q / m
  G = k*q*Q/m;

  
% Initialize u matrix
   u0 = [x0 vx0 y0 vy0];
% Time interval: equal time increments for ode45  
  tSpan = linspace(0,tMax,nT);
 
  options = odeset('RelTol',1e-6,'AbsTol',1e-6);

   
% CALCULATIONS  =======================================================

  [t,u] = ode45(@EqM, tSpan, u0, options);
  
% X and Y components of displacement and velocity as functions of time t 
  x =  u(:,1);
  vx = u(:,2);
  y =  u(:,3);
  vy = u(:,4);

% Instantaneous values  -----------------------------------------------
% Radius of orbit
   R = sqrt(x.^2 + y.^2);
% Velocity
   v = sqrt(vx.^2 + vy.^2);
% Angular momentum  
  L = x.*vy - y.*vx;
% Energies: kinetic, potential, total  
  K = 0.5*m.*v.^2;
  U = m*G./R;
  E = K + U;
  
% Exit slope of trajectory
  theta = atan2d((y(end)-y(end-150)),(x(end)-x(end-150)));


% GRAPHICS ============================================================

 figure(1)
  pos = [0.05 0.06 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
 
  col = [165 42 42]./255;
  Hplot = plot(0,0,'o');
  set(Hplot,'markersize',8,'markerfacecolor',col, ...
      'markeredgecolor',col)
  hold on
  
 % plot(u(:,1),u(:,3),'bo')
  plot(u(:,1)./S,u(:,3)./S,'b','linewidth',LW)
  Hplot = plot(u(1,1)./S,u(1,3)./S,'o');
  set(Hplot,'markerfacecolor','g','markeredgecolor', 'g', ...
      'markersize',5)
   Hplot = plot(u(end,1)./S,u(end,3)./S,'o');
   set(Hplot,'markerfacecolor','r','markeredgecolor', 'r', ...
       'markersize',5)
 

  xlabel('x  [fm]')
  ylabel('y  [fm]')
  set(gca,'fontsize',12)
 
%  xlim(xGrid)
%  ylim(yGrid)
%  set(gca,'xtick',xtickGrid);
%  set(gca,'ytick',ytickGrid);
 
 %axis equal
 axis equal
 
 tm1 = '\theta = '; tm2 = num2str(theta,'%2.1f'); tm3 = '  deg  ';
 tm4 = '   R_{min} =  '; tm5 = num2str(min(R)*1e15,'%2.1f');
 tm6 = '  fm';
 tm = [tm1 tm2 tm3 tm4 tm5 tm6];
 title(tm) 
      
 
 grid on

 
if flag2 == 1 
    
figure(2)  %  ---------------------------------------------------------
  pos = [0.35 0.06 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
       
subplot(4,1,1);
  plot(t,x./S,'b','linewidth',LW)
  hold on
  plot(t,y./S,'r','linewidth',LW)
  plot(t,R./S,'m','linewidth',LW)
  grid on
  ylabel('x, y, R  [ fm ]')
  legend('x','y','R','orientation','horizontal','location','southeast','box','off')
  set(gca,'fontsize',12)
  
  
subplot(4,1,2);
  plot(t,vx,'b')
  hold on
  plot(t,vy,'r','linewidth',LW)
  plot(t,v,'m','linewidth',LW)
  ylabel('v  [ m/s ]')
  legend('v_x','v_y','v','orientation','horizontal','location','southeast','box','off')
  grid on
  set(gca,'fontsize',12)
  
subplot(4,1,3);
  plot(t,L,'m','linewidth',LW)
 % ylim([-2*abs(max(L)) 2*abs(max(L))])
  ylim([-1 1])
  grid on
  ylabel('L  [ N.s ]')
  set(gca,'fontsize',12)
   
subplot(4,1,4);
  plot(t,K./(1e6*e),'r','linewidth',LW)
  hold on
  plot(t,U./(1e6*e),'b','linewidth',LW)
  plot(t,E./(1e6*e),'m','linewidth',LW)
  grid on
  xlabel('time t [ s ]')
  ylabel('K, U, E  [ MeV ]')
  legend('K','U','E','orientation','horizontal','location','east','box','off');
  set(gca,'fontsize',12)

end
 toc

 
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global G
  uDot = zeros(4,1);
  R3 = (u(1)^2 + u(3)^2)^1.5;
  G = 6.2499;
  uDot(1) = u(2);
  uDot(2) = G*u(1)/R3;
  uDot(3) = u(4);
  uDot(4) = G*u(3)/R3;
end

