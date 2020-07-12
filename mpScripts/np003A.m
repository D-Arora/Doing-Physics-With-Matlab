% np003A.m

% Physics of Neurones
% [1D] Dynamical Models: EXAMPLES
%   CELL 1: Phase Portrait Plot
%            Define function in CELL 1
%   CELL 2: Time Evolution of State Variable: Solve equation of motion
%            Define function in CELL 4 at end of Script
%   CELL 3: Bifurcation Diagram 
%   CELL 4: Define function for oder45

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/np001.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191102

close all
clc
clear

global b
tic

% =====================================================================
% #1   PHASE PORTRAIT PLOT xDot vs x 
% =====================================================================

% Bifurcation Prameter
  b = 1;
% x grid 

  xMin = -2;
  xMax = 2;
  Nx = 9999;
  x = linspace(xMin,xMax, Nx);
  dx = x(2)-x(1);
   
% Define function
   xDot = b.*x - x.^3; 
%  xDot = -1 + x.^2;
%  xDot = x - x.^3;


% Find equilibrium points xDot = 0  
%   zCount     Number of equilibrium points  
%   zIndex     Indices for equilibrium points 
%   ep = -1    Stable equilibrium point
%   ep = + 1   Unstable equilibrium point

  zCount = 0;
  
  for c = 1:Nx-1
    zSign = xDot(c)*xDot(c+1);
      if zSign <= 0
         zCount = zCount + 1;
         zIndex(zCount) = c;
      end
  end
  
% Equilibrium Points 
  zLen = length(zIndex);
  if zLen > 0
    disp('Equilibrium Points')
    disp(x(zIndex));
% Sign of tangent at equilibrium point
    ep = zeros(zLen,1);
      for c =  1 : zLen
        ep(c) = sign( xDot(zIndex(c)+1) - xDot(zIndex(c)-1) );
        lambda(c) = ( xDot(zIndex(c)+1) - xDot(zIndex(c)-1) )/(2*dx);
      end
  end
    disp('  ')
    disp('Eigenvalues Lambda')
    disp(lambda);
  
figure(1)  % Phase Portrait Plot
  pos = [0.05 0.05 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
 
  hold on
  xP = x; yP = xDot;
  plot(xP,yP,'b','linewidth',LW)
  
for c = 1 : zLen
  if ep(c) == -1; col = 'r'; end
  if ep(c) ==  1; col = 'g'; end
  Hplot = plot(xP(zIndex(c)),yP(zIndex(c)),'o');
  set(Hplot,'markersize',8,'markerfacecolor',col,'markeredgecolor',col);
end

  box on
  grid on
  xlabel('x  ')
  ylabel('dx / dt ')
  set(gca,'fontsize',14)
  
  tm = sprintf('bifurcation parameter   b = %2.1f ',b);
  title(tm,'fontweight','normal')

 
%% =====================================================================
%  #2   TIME EVOLUTION OF STATE VARIABLE
%  =====================================================================

% Initial conditions
   u0 = -2;
   
% Time domain
   tMax = 5;
   tSpan = linspace(0,tMax,9999);
   dt = tSpan(2)-tSpan(1);
  
% Solve ODE: State variable u
   options = odeset('RelTol',1e-6,'AbsTol',1e-6);
   [t,u] = ode45(@EqM, tSpan, u0, options);  


figure(2)  
  pos = [0.05 0.56 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
  
  hold on
  xP = t; yP = u;
  plot(xP,yP,'b','linewidth',LW)
  
  ylim([xMin xMax])
  box on
  grid on
  xlabel('t  ')
  ylabel('x ')
  set(gca,'fontsize',14)

  tm = sprintf('bifurcation parameter   b = %2.1f ',b);
  title(tm,'fontweight','normal')
  
%% =====================================================================
%  #3   BIFURCATION DIAGRAM
%  =====================================================================
 
% Bifurcation parameter  a
N = 599;
aMin = -2.5; aMax = 2.5;
a1 = linspace(aMin, 0, N);
x1 = zeros(N,1);

a2 = linspace( 0,aMax, N);
x2 =  zeros(N,1);
x3 =  sqrt(a2);
x4 = -sqrt(a2);

figure(3)  
  pos = [0.45 0.56 0.25 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  LW = 2;
  
  hold on
  xP = a1; yP = x1;
  plot(xP,yP,'r','linewidth',LW)
  xP = a2; yP = x2;
  plot(xP,yP,'g','linewidth',LW)
  xP = a2; yP = x3;
  plot(xP,yP,'r','linewidth',LW)
  xP = a2; yP = x4;
  plot(xP,yP,'r','linewidth',LW)
  
  box on
  grid on
  xlabel('bifurcation parameter  b  ')
  ylabel('equilibium  x ')
  set(gca,'fontsize',14)


%%

toc

 
%% =====================================================================
%  #4   FUNCTIONS
%  =====================================================================

function uDot = EqM(t,u)
  global b
  uDot = zeros(1,1);
 
   uDot(1) = b*u(1) - u(1)^3;
%  uDot(1) = -1 + u(1)^2;
%  uDot(1) = u(1) - u(1)^3;
   
end  
