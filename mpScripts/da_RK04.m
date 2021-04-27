% da_RK04.m


% RUNGE_KUTTA METHOD
% Solving Initial Value Ordinary Differential Equations 
% Sets of First-Order Equations
% Higher order ODE

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/RK.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
%  200522     Matlab 2019b

clear 
close all
clc


global m b k

% INPUT SECTION  ======================================================
% >>>  Number of grid points 
  N = 201;
  
% >>>  t variable
  tMin = 0;
  tMax = 4.8;
  
% >>>  Initial value for x and y variable
  x = zeros(1,N); y = zeros(1,N); 
  x(1) = 0;
  y(1) = 0.1;
  
% Damped harmonic oscillator input parameters
% >>> Mass  [1.0]
  m = 1.0;
% >>> Period  [1.0]
  T = 1.0;
% >>> Damping constant  [0 to 20]
  b = 25;
% Spring constant
  k = 4*pi^2/(m*T^2);
% Critical damping constant
  b_critical = 2*m*sqrt(k/m);
  

% SETUP  ==============================================================  
% t interval and step size
  t = linspace(tMin,tMax,N);
  h = t(2) - t(1);
% xDot & yDot variables as a function of t
%      at x(nX) & y(nY) (calls function XDOT & YDOT)   
  nX = 1; nY = 1;
  xDot = XDOT(t,x(nX), y(nY));
  yDot = YDOT(t,x(nX),y(nY));
  
  
% Compute y and yDdot values as a function of time =====================
for n = 1 : N-1
   kx1 = h*XDOT(t(n), x(n), y(n));
   ky1 = h*YDOT(t(n), x(n), y(n));
   
   kx2 = h*XDOT(t(n) + h/2, x(n) + kx1/2, y(n) + ky1/2);
   ky2 = h*YDOT(t(n) + h/2, x(n) + kx1/2, y(n) + ky1/2);
   
   kx3 = h*XDOT(t(n) + h/2, x(n) + kx2/2, y(n) + ky2/2);
   ky3 = h*YDOT(t(n) + h/2, x(n) + kx2/2, y(n) + ky2/2);
   
   kx4 = h*XDOT(t(n) + h,   x(n) + kx3,   y(n) + ky3);
   ky4 = h*YDOT(t(n) + h,   x(n) + kx3,   y(n) + ky3);
   
   x(n+1) = x(n) + (kx1+2*kx2+2*kx3+kx4)/6;
   y(n+1) = y(n) + (ky1+2*ky2+2*ky3+ky4)/6;
end




% xDOT & yDot variables (calls functions XDOT YDOT) 
  xDot = XDOT(t,x(1),y(1));
  yDot = YDOT(t,x(1),y(1));
 

% [2D] plots  =========================================================   
   Y = linspace(1,3,N);
  [tt, yy] = meshgrid(t,Y);
  YYDot = YDOT(tt,yy);
  
  
% GRAPHICS  ===========================================================
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.10 0.10 0.25 0.6]);
   set(gcf,'color','w');

subplot(3,1,1)
   xP = t; yP = x;
   plot(xP,yP,'b','linewidth',2)
   grid on
   xlabel('t')
   ylabel('v = dy/dt')
   txt = sprintf('m = %2.1f  k = %2.1f   b = %2.1f  b_{critical} = %2.1f  ', ...
       m,k,b, b_critical);
   title(txt,'fontweight','normal')
   set(gca,'fontsize',12)
   
subplot(3,1,2)   
   yP = y;
   plot(xP,yP,'r','linewidth',2);
   xlabel('t')
   ylabel('y')
   grid on
   set(gca,'fontsize',12)
   
subplot(3,1,3)
   xP = x; yP = y;
   plot(xP,yP,'b','linewidth',2)
   hold on
   Hplot = plot(xP(1),yP(1),'go');
   set(Hplot,'markersize',8,'markerfacecolor','g')
   Hplot = plot(xP(N),yP(N),'ro');
   set(Hplot,'markersize',8,'markerfacecolor','r')
   grid on
   xlabel('v = dy/dy')
   ylabel('y')
  % txt = sprintf('y = %3.1f  ',y(nY));
  % title(txt,'fontweight','normal')
   set(gca,'fontsize',12)   
   

   
% FUNCTION  ===========================================================

function f = XDOT(t,x,y)
  global m b k
  f = -(k/m).*y - (b/m).*x;
end


function g = YDOT(t,x,y)
    g = x;
end





