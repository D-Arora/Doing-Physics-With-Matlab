% cnsEpiletor.m

%  Epiletor Model
%    Solved using ode45 
%    Solution Variables: y(:,1)    x1
%                        y(:,2)    y1
%                        y(:,3)    z
%                        y(:,4)    x2
%                        y(:,5)    y2
%                        y(:,6)    u

% Most of the parameters are specified in the INPUT SECTION
% Simulation time: tSpan
% Initial conditions:  y0 = [x1(0) y1(0) z(0) x2(0) y2(0) u(0)]
% Model constants  k
                  
% For different models, it may be necessary to make changes to the Script


% Ian Cooper
% email: matlabvisualphysics@gmail.com
% DOING PHYSICS WITH MATLAB 
%     https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%     https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cnsHindmarshA.m.pdf
% Date: 211216 / Matlab version: R2021b
  
close all
clear
clc

tic


% INPUTS  ==============================================================

% Initial conditions 
  x10 = 0; y10 = -5; z0 = 3; x20 = 0; y20 = 0; u0 = 0.1*x10;

  y0 = [x10; y10; z0; x20; y20; u0];

% Time domain
  tMax = 200;  N = 5001;
  tSpan = linspace(0,tMax,N);


% CALCULATION SECTION ===================================================

% Solve differential equations using ode45
%  opts = odeset('RelTol',1e-10);
  [t,y] = ode45(@(t,y) FNode(t,y), tSpan,y0);

  x1 = y(:,1); y1 = y(:,2); z = y(:,3);
  x2 = y(:,4); y2 = y(:,5); u = y(:,6);


% GRAPHICS  ===========================================================

 figure(1)   % Phase space plot: v vs w   ---------------------------
   pos = [0.05 0.05 0.5 0.50];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 14;  % fontsize
   hold on
   box on

  subplot(3,2,1)
    xP = t; yP = x1;
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('x_1')
    xlabel('t')
    set(gca,'fontsize',FS)

 subplot(3,2,2)
    xP = t; yP = y1;
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('y_1')
    xlabel('t')
    set(gca,'fontsize',FS)

 subplot(3,2,3)
    xP = t; yP = x2;
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('x_2')
    xlabel('t')
    set(gca,'fontsize',FS)

 subplot(3,2,4)
    xP = t; yP = y2;
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('y_2')
    xlabel('t')
    set(gca,'fontsize',FS)  

 subplot(3,2,5)
    xP = t; yP = z;
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('z')
    xlabel('t')
    set(gca,'fontsize',FS) 

 subplot(3,2,6)
    xP = t; yP = u;
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('u')
    xlabel('t')
    set(gca,'fontsize',FS)   

% --------------------------------------------------------------------
 figure(2)   % Phase space plot: v vs w  
   pos = [0.60 0.05 0.25 0.250];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 14;  % fontsize
   hold on
   box on

  xP = x1; yP = y1;
  plot(xP,yP,'b','linewidth',2)
  grid on
  ylabel('y_1')
  xlabel('x_1')
  set(gca,'fontsize',FS)   

% --------------------------------------------------------------------
 figure(3)   % Phase space plot: v vs w  
   pos = [0.60 0.4 0.25 0.350];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 14;  % fontsize
   hold on
   box on

  zP = z; yP = x1; xP = y1;
  plot3(xP,yP,zP,'b','linewidth',2)
  grid on
  ax = gca;
  ax.BoxStyle = 'full';
  zlabel('z')
  ylabel('x_1')
  xlabel('y_1')
  set(gca,'fontsize',FS)     
  view(-28,45)

  % --------------------------------------------------------------------
 figure(4)   % Phase space plot: v vs w  
   pos = [0.05 0.55 0.25 0.250];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 14;  % fontsize
   hold on
   box on

   xP = t; yP = x1+x2;
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('x_1 + x_2')
   xlabel('t')
   set(gca,'fontsize',FS)   

% ====================================================================    
  toc
  

% FUNCTIONS ===========================================================

function dydt = FNode(t,y)

  x1 = y(1); y1 = y(2); z = y(3);
  x2 = y(4); y2 = y(5); u = y(6);
  
  % x1_dot  
  Iext1 = 3.1 + 0.025*randn(1,1);   % white noise added
  dydt(1) = y1 - fn1(x1,x2,z) - z + Iext1;

  % y1_dot
  c = 1; d = 5;
  dydt(2) = c - d*x1^2 - y1 ;
  
  % z_dot
  r = 3.5e-4; s = 4; x0 = -1.6;
  if z >= 0
    dydt(3) = r * ( s*(x1 - x0) - z ) ;
  else
   dydt(3) = r * ( s*(x1 - x0) - z - 0.1*z^7 );
  end   

  % x2_dot   
  Iext2 = 0.45 + 0.25*randn(1,1);   % white noise added
  tau2 = 10;
  dydt(4) = -y2 + x2 - x2^3 + Iext2 + 2*u - 0.3*(z - 3.5);

  % y2_dot
  dydt(5) =(1/tau2) * ( -y2 +fn2(x2) );

  % u_dot
  gamma = 0.01;
  dydt(6) = -gamma*(u - 0.1*x1);

  dydt = dydt';

end  

function s12 = fn1(x1,x2,z)
  a = 1; b = 3;  m = 0;    % m = 0.5
 
  if x1 < 0
    s12 = a*x1^3 - b*x1^2;
  else
    s12 = -(m - x2 + 0.6*(z - 4)^2) * x1;
  end

end

function s2 = fn2(x2)
  a2 = 6;
  s2 = 0;
  if x2 >= -0.25; s2 = a2*(x2 + 0.25); end

end




