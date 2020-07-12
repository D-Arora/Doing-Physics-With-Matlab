% modAlphaSc.m

% ALPHA PARTICLE SCATTERING
% COULOMB'S LAW    
%  Arbituary units kqQ = 1
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


% close all
clc
clear

tic


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Default values
%   x0 = -5   vx0 = 2   y0 = 1.0   vy0 = 0
%   tMax = 8  nT = 181
%   Larger values of nT give more accurate results

% Matrix u (X and Y displacements and velocities)
%  u(1) = x   u(2) = vx   u(3) = y; u(4) = vy = 0
% Initial conditions: t = 0
% Impact parameter >>>
   y0  = 0;
   
   x0  = -5;
   vx0 = 2;
   vy0 = 0;
    
 % Max time interval for simulation / number of calculations
  tMax = 6;
  nT = 531; 

% Grid dimensions
  Grid = 5;
  xGrid = [-Grid Grid];
  yGrid = xGrid;
  tickGrid = -Grid:Grid;

% Plot figure 2   (1) yes /  (2)  no
  flag2 = 1;
  if flag2 == 1; close all; end
  
% Setup  ==============================================================
  
% Mass of object  
  m = 1;

% Initialize u matrix
   u0 = [x0 vx0 y0 vy0];
% Time interval: equal time increments for ode45  
  tSpan = linspace(0,tMax,nT);
 
  options = odeset('RelTol',1e-6,'AbsTol',1e-4);

   
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
  U = m./R;
  E = K + U;
  

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
  plot(u(:,1),u(:,3),'b','linewidth',LW)
  Hplot = plot(u(1,1),u(1,3),'o');
  set(Hplot,'markerfacecolor','g','markeredgecolor', 'g', ...
      'markersize',5)
   Hplot = plot(u(end,1),u(end,3),'o');
   set(Hplot,'markerfacecolor','r','markeredgecolor', 'r', ...
       'markersize',5)
 

  xlabel('x')
  ylabel('y')
  set(gca,'fontsize',12)
 
 xlim(xGrid)
 ylim(yGrid)
 set(gca,'xtick',tickGrid);
 set(gca,'ytick',tickGrid);
 axis square
 

 grid on

 
if flag2 == 1 
    
figure(2)  %  ---------------------------------------------------------
  pos = [0.35 0.06 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
       
subplot(4,1,1);
  plot(t,x,'b','linewidth',LW)
  hold on
  plot(t,y,'r','linewidth',LW)
  plot(t,R,'m','linewidth',LW)
  grid on
  ylabel('displacement')
  legend('x','y','R','orientation','horizontal','location','southeast','box','off')
  set(gca,'fontsize',12)
  
  
subplot(4,1,2);
  plot(t,vx,'b')
  hold on
  plot(t,vy,'r','linewidth',LW)
  plot(t,v,'m','linewidth',LW)
  ylabel('velocity')
  legend('v_x','v_y','v','orientation','horizontal','location','southeast','box','off')
  grid on
  set(gca,'fontsize',12)
  
subplot(4,1,3);
  plot(t,L,'m','linewidth',LW)
 % ylim([-2*abs(max(L)) 2*abs(max(L))])
  ylim([-1 1])
  grid on
  ylabel('ang momentum')
  set(gca,'fontsize',12)
   
subplot(4,1,4);
  plot(t,K,'r','linewidth',LW)
  hold on
  plot(t,U,'b','linewidth',LW)
  plot(t,E,'m','linewidth',LW)
  grid on
  xlabel('time t [a.u.]')
  ylabel('energy')
  legend('K','U','E','orientation','horizontal','location','east','box','off');
  set(gca,'fontsize',12)

end
 toc

 
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  uDot = zeros(4,1);
  R3 = (u(1)^2 + u(3)^2)^1.5;
  
  uDot(1) = u(2);
  uDot(2) = u(1)/R3;
  uDot(3) = u(4);
  uDot(4) = u(3)/R3;
end

