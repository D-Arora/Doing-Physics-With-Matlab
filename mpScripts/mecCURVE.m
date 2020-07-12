% mecCURVE.m

% Simulation of a car going around a curved road.

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% VISUAL PHYSICS ONLINE
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod5new/mod52B.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170517

clear 
close all
clc


% =====================================================================





% SETUP: initial values ===============================================
% Time / time increment / time steps 
  tMax = 40; 
  nT = 120;
  t = linspace(0,tMax,nT)';
  dt = t(2) - t(1);
% Radius of curve
  R = 50.*ones(nT,1);
% (x,y) coordinates of car
  x = R(1).*ones(nT,1);
  y = zeros(nT,1);
% tangential velocity of car
  vT = zeros(nT,1);
% mass of car
  m = 1000;
% coefficent of friction
  mu = 0.6;
% acceleration due to gravity
  g = 9.81;
% tangentiual acceleration
  aT = 0.5; 
% angular acceleration (alpha)  [rad/s^2]
  A = zeros(nT,1);
% angular speed (omega)  [rad/s]
  W = zeros(nT,1);
% Angular displacement (theta)  [rad]
  T = zeros(nT,1);

%Initial condictions  ================================================
   
alpha(1) = aT/R(1); 
R(1) = 50;
x(1) = R(1); y(1) = 0;
vT(1) = 0;
A(1) = aT/R(1);
W(1) = 0;
T(1) = 0;

% Path of car calculation  ============================================

for ct = 2 : nT-1
  A(ct) = aT/R(ct);
  W(ct) = W(ct-1) + A(ct-1)*dt;
  T(ct) = T(ct-1) + W(ct-1)*dt + 0.5*A(ct-1)*dt^2; 
  
  vT(ct) = R(ct)*W(ct);
  x(ct) =  R(ct)*cos(T(ct));
  y(ct) =  R(ct)*sin(T(ct));
  
  R(ct+1) = vT(ct)^2/(mu*g);
    if R(ct+1) < R(1)
       R(ct+1) = R(1);   
    end   
   
end


% GRAPHICS ============================================================

   figure(1)
   pos = [0.1 0.2 0.25 0.35];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   box on
   
  % xP = t; yP = y;
   xP = x(1:nT-1); yP = y(1:nT-1);
   plot(xP,yP)
   axis([-80 80 -80 80])
   axis equal
