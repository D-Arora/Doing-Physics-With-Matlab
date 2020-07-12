% math_ode_03.m

% Solving second order differential equations with the ode45 solver:
%    van der Pol Oscillator.

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190228 / Matlab version R2018b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/math_ODE_A.htm


close all
clear
clc

tic

% INPUTS  =============================================================
% default values and units  [ ]

% coefficient mu   [2]
   mu = 2;
% Driving force amplitude   [0]
   A = 0;
% Driving force period [10] 
   Tin = 10;
% Initial conditions [y = 0 m   v = 0 m/s]
   u0 = [0; 0.5];
% Time span  [0 20 s]
   tSpan = [0 50];

% For equal time increments and specify the number of time points use:  
%    tSpan = linspace(0, 20, 1001);

% Relative tolerance for ODE solver  [1e-6]
  RelTol = 1e-6;
  
  
% CALCULATIONS  ======================================================= 

% Constants for differential equation
   K(1) = mu;
   K(2) = A;
   K(3) = Tin;
   
  options = odeset('RelTol',RelTol); 
  [t, SOL] = ode45(@(t,u) FNode(t,u,K), tSpan, u0, options);

% Displacement  [m]
   y = SOL(:,1);
% Velocity  [m/s]
   v = SOL(:,2);
% Acceleration  [m/s^2]
  nmax = length(t);
  a = zeros(nmax,1);
  a(1) = (v(2) - v(1))/(t(2)- t(1));
  for n = 2 : nmax
     a(n) = (v(n)- v(n-1))/((t(n) - t(n-1)));
  end

  
% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.29 0.69];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
subplot(3,1,1)   
xP = t; yP = y;   
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('y  [ m ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',12)
   
subplot(3,1,2)   
   xP = t; yP = v;   
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('v  [ m/s ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',12)
      
subplot(3,1,3)   
   xP = t; yP = a;   
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('a  [ m/s^2 ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',12)   
   
   
figure(2)
   pos = [0.35 0.05 0.29 0.60];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
subplot(2,1,1)   
xP = y; yP = v;   
   plot(xP,yP,'b','linewidth',2)
   hold on
   Hplot = plot(xP(1),yP(1),'go');
   set(Hplot,'markersize',10,'markerfacecolor','g')
   grid on
   xlabel('y  [ m ]')
   ylabel('v  [ m/s ]')
   set(gca,'fontsize',12)
   
subplot(2,1,2)   
   xlim([0 100])
   ylim([20 130])
   
   t1 = 'initial displacement y(0) =  ';
   t2 = num2str(u0(1),'%2.4f  m \n');
   tm = [t1 t2];
   text(10,90,tm,'fontsize',14)
   
   t1 = 'initial velocity  v(0)  =  ';
   t2 = num2str(u0(2),'%2.4f  m.s^{-1} \n');
   tm = [t1 t2];
   text(10,75,tm,'fontsize',14)
   
   t1 = 'model coefficent  \mu  =  ';
   t2 = num2str(mu,'%2.4f   \n');
   tm = [t1 t2];
   text(10,60,tm,'fontsize',14)
   
   t1 = 'driving force amplitude  A  =  ';
   t2 = num2str(A,'%2.4f  m.s^{-2} \n');
   tm = [t1 t2];
   text(10,40,tm,'fontsize',14)
   
   t1 = 'driving force period  T_{in}  =  ';
   t2 = num2str(Tin,'%2.4f  s \n');
   tm = [t1 t2];
   text(10,20,tm,'fontsize',14)
      
   tm = 'd^2y/dt^2  =  - y + \mu (1- y^2) dy/dt + A cos(2\pi t / T_{in})';
   text(10,130,tm,'fontsize',14)
   
   t1 = 'Number of time points  =  ';
   t2 = num2str(length(t),'%4.0f  \n');
   tm = [t1 t2];
   text(10,115,tm,'fontsize',14)
   axis off
   
toc
   
% FUNCTIONS  ==========================================================

function du = FNode(t,u,K)
  y = u(1);
  yDot = u(2);
  du = zeros(2,1);
 % ODE: First derivative dy/dt
   du(1) = yDot;
 % ODE: Second derivative  d2y/dt2
    du(2) = -y + K(1)*(1-y^2)*yDot + K(2)*cos(2*pi*t/K(3)) ;
end