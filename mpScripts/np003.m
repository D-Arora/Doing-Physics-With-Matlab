% np003.m

% Physics of Neurones
% [1D] Dynamical Models: EXAMPLES
%   CELL 1: Phase Portrait Plot
%            Define function in CELL 1
%   CELL 2: Time Evolution of State Variable: Solve equation of motion
%            Define function in CELL 3 at end of Script

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

tic

% =====================================================================
% #1   PHASE PORTRAIT PLOT xDot vs x 
% =====================================================================
% x grid 
  xMin = -2;
  xMax = 2;
  Nx = 9999;
  x = linspace(xMin,xMax, Nx);

% Define function
%  xDot = -1 + x.^2;
   xDot = x - x.^3;
   
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

  zLen = length(zIndex);
  if zLen > 0
  % Equilibrium Points
    disp('Equilibrium Points')
    disp(x(zIndex));
    
    ep = zeros(zLen,1);
      for c =  1 : zLen
        ep(c) = sign( xDot(zIndex(c)+1) - xDot(zIndex(c)-1) );
      end
  end

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
  ylabel('x_{dot} ')
  set(gca,'fontsize',14)


 
%% =====================================================================
%  #2   TIME EVOLUTION OF STATE VARIABLE
%  =====================================================================

% Initial conditions
   u0 = 2;
   
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

  
%%
toc

 
%% =====================================================================
%  #3   FUNCTIONS
%  =====================================================================

function uDot = EqM(t,u)
  
  uDot = zeros(1,1);
  
%  uDot(1) = -1 + u(1)^2;
   uDot(1) = u(1) - u(1)^3;
end  
