% mecSatellite.m

% Motion of an object of mass m in the gravitational field of object mass M
%     m << M
% Universal Law of Gravitation:
%   solving the equation of motion using the ode45 function

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/mecSatelliteA.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ON LINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod5new.htm

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190919   Matlab 2019b


close all
clc
clear

tic

global GM


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Default values
%   x0 = 1   vx0 = 0   y0 = 0   vy0 = 1
%   tMax = 6.3  nT = 81
%   Larger values of nT give more accurate results

% Matrix u (X and Y displacements and velocities)
%  u(1) = x   u(2) = vx   u(3) = y; u(4) = vy (vy > 0)
% Initial conditions: t = 0
  x0  = 1;
  vx0 = 0;
  y0  = 0;
  vy0 = 1;
    
 % Max time interval for simulation / number of calculations
  tMax = 6.3;
  nT = 81; 


% Setup  ==============================================================
  
% Earth's radius
  R_E = 6.38e6;  
% Earth's mass  
  M_E = 5.98e24; 
% Universal gravitation constant 
  G = 6.67e-11;   
% Mass of object  
  m = 1;
% Use arbitrary units: Set G*M = 1  
%  GM = G*M_E;
  GM = 1;
  
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
  U = - GM*m./R;
  E = K + U;
  
% Orbital Period (linear interpolation)
  flagT = 0;
  Tindex = find(y<0,1);
  if Tindex > 0
    flagT = 1;
    mT = (y(Tindex) - y(Tindex-1))/(t(Tindex) - t(Tindex-1));
    bT = y(Tindex) - mT*t(Tindex);
    T = -2*bT/mT;
    dT = (t(Tindex)- t(Tindex-1))/2;
  end
  
% Circular orbit ------------------------------------------------------ 
  x0 = u0(1); y0 = u0(3);
  R0 = sqrt(x0^2 + y0^2);
% Orbital velocity 
   vCorb = sqrt(GM/R0);
% Period  
   TCorb = 2*pi*R0/vCorb;

% Escape velocity 
   vESC = sqrt(2)*vCorb;
   

% GRAPHICS ============================================================

 figure(1)
  pos = [0.05 0.06 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
 
  col = [165 42 42]./255;
  Hplot = plot(0,0,'o');
  set(Hplot,'markersize',18,'markerfacecolor',col, ...
      'markeredgecolor',col)
  hold on
  
  plot(u(:,1),u(:,3),'bo')
 % plot(u(:,1),u(:,3),'b')
  Hplot = plot(u(1,1),u(1,3),'o');
  set(Hplot,'markerfacecolor','g','markeredgecolor', 'g', ...
      'markersize',10)
   Hplot = plot(u(end,1),u(end,3),'o');
   set(Hplot,'markerfacecolor','r','markeredgecolor', 'r', ...
       'markersize',10)
   
  tm1 ='x     ';
  tm2 = '  v_{x0} = ';  tm3 = num2str(vx0,'%2.2f');
  tm4 = '  v_{y0} = ';  tm5 = num2str(vy0,'%2.2f');
  tm = [tm1 tm2 tm3 tm4 tm5];
  xlabel(tm)
  ylabel('y')
  if flagT == 1
      tm1 = 'T = '; tm2 = num2str(T,'%2.2f'); tm3 = '  \pm  ';
      tm4 = num2str(dT,'%2.2f');
      tm = [tm1 tm2 tm3 tm4];
      title(tm) 
  end
 set(gca,'fontsize',12)
 axis equal
 grid on
 
figure(2)  %  ---------------------------------------------------------
  pos = [0.35 0.06 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
       
subplot(4,1,1);
  plot(t,x,'b','linewidth',LW)
  hold on
  plot(t,y,'r','linewidth',LW)
  plot(t,R,'m','linewidth',LW)
  grid on
  ylabel('displacement')
  legend('x','y','R','orientation','horizontal','location','north','box','off')
  set(gca,'fontsize',12)
  
  tm1 = 'R_0 = ';         tm2 = num2str(R0,'%2.2f');
  tm3 = '  v_{Corb} = ';  tm4 = num2str(vCorb,'%2.2f');
  tm5 = '  T_{Corb} = ';  tm6 = num2str(TCorb,'%2.2f');
  tm7 = '  v_{ESC} = ';   tm8 = num2str(vESC,'%2.2f');
  tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8];
  title(tm) 
  
subplot(4,1,2);
  plot(t,vx,'b')
  hold on
  plot(t,vy,'r','linewidth',LW)
  plot(t,v,'m','linewidth',LW)
  ylabel('velocity')
  legend('v_x','v_y','v','orientation','horizontal','location','north','box','off')
  grid on
  set(gca,'fontsize',12)
  
subplot(4,1,3);
  plot(t,L,'m','linewidth',LW)
  ylim([0 2])
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
  legend('K','U','E','orientation','horizontal','location','north','box','off');
  set(gca,'fontsize',12)
  
 toc

 
% FUNCTIONS  ==========================================================  
  
function uDot = EqM(t,u)
  global GM
  uDot = zeros(4,1);
  R3 = (u(1)^2 + u(3)^2)^1.5;
  
  uDot(1) = u(2);
  uDot(2) = -GM*u(1)/R3;
  uDot(3) = u(4);
  uDot(4) = -GM*u(3)/R3;
end

