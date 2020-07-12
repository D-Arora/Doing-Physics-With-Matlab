% math_ode_04.m

% Solving second order differential equations with the ode45 solver:
%    Damped, driven harmonic oscillator.

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190228 / Matlab version R2018b

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/math_ODE_A.htm

close all
clear
clc

tic

% INPUTS  =============================================================
% default values and units  [ ]

% mass  [ m = 10 kg]
   m = 0.506606; 
% spring constant [k = 20 N/m]  
   k = 20;
% Damping constant  [b = 2  kg/s]
   b = 1;
% Amplitude of driving force  [A = 10 N]
   A = 2;
% frequency of driving force  [fD = 1  hz]
   fD = 1;
% Initial conditions [y = 0 m   v = 0 m/s]
   u0 = [1; 0];
% Time span  [0 20 s]
  tSpan = [0 10];

% For equal time increments and specify the number of time points use:  
%    tSpan = linspace(0, 20, 1001);

% Relative tolerance for ODE solver  [1e-6]
  RelTol = 1e-6;
  
  
% CALCULATIONS  ======================================================= 

% Resonance frequency and driving frequency  [Hz and rad/s]
  w0 = sqrt(k/m);
  f0 = w0 / (2*pi);
  wD = 2*pi*fD;
%   
% Constants for differential equation
%   yDotDot = -(k/m)y - (b/m)yDot + (A/m)sin(wDt);
%   yDotDot =   K(1)y +  K(2)yDot +  K(3)sin(K(4)t);
   K(1) = -k/m;
   K(2) = -b/m;
   K(3) = A/m;
   K(4) = wD;
 
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
   
   t1 = 'm  =  ';
   t2 = num2str(k,'%2.4f  kg \n');
   tm = [t1 t2];
   text(10,90,tm,'fontsize',14)
   
   t1 = 'k  =  ';
   t2 = num2str(m,'%2.4f  N.m^{-1} \n');
   tm = [t1 t2];
   text(50,90,tm,'fontsize',14)
   
   t1 = 'f_0  =  ';
   t2 = num2str(f0,'%2.4f  Hz \n');
   tm = [t1 t2];
   text(10,75,tm,'fontsize',14)
   
   
   t1 = 'b  =  ';
   t2 = num2str(b,'%2.4f  kg.s^{-1} \n');
   tm = [t1 t2];
   text(10,60,tm,'fontsize',14)
   
   t1 = 'A  =  ';
   t2 = num2str(A,'%2.4f  N \n');
   tm = [t1 t2];
   text(10,45,tm,'fontsize',14)
   
   t1 = 'f_D  =  ';
   t2 = num2str(fD,'%2.4f  Hz \n');
   tm = [t1 t2];
   text(50,45,tm,'fontsize',14)
   
   tm = 'm d^2y/dt^2  =  - \omega^2 y - b dy/dt + A sin(\omega_D t)';
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
    du(2) = K(1)*y + K(2)*yDot + K(3)*sin(K(4)*t);
end