% sp_mod601.m
% Motion of a particle in uniform gravitational or electric field
% 18 aug 2017
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/


clear all
close all
clc


% INPUTS =================================================================
m = 10;      % mass of particle
u = 10;      % initial velocity (magnitude)
A = 75;      % intial velocity (angle) 
tMin = 0;    % time
tMax = 5;
Nt = 2000;
g = 5;      % acc due to gravity


flagF = 1;   % 1: gravitational field / 2: electric field

q = 1;       % particle's charge
V = 100;     % potential difference between charged plates;
d = 1e-3;    % distance bewtween plates


% Calculations ===========================================================
   t = linspace(tMin,tMax,Nt);   % time
% Inital velocityies
   ux = u*cosd(A);
   uy = u*sind(A);
% X & Y accelerations
  ax = 0;  
  if flagF == 1; ay = -g; end
  if flagF == 2
    E = V / d;   % electric field strength
    F = q*E;     % electric force on charged particle
    ay = F / m;   % accleration of particle
  end
    ayG = ay .* ones(Nt,1);
    axG = ax .* ones(Nt,1);
      
% Displacements
    sx = ux .* t;
    sy = uy .* t - (0.5*g) .* t.^2;
 
% velocites
     vx = ux .* ones(Nt,1);
     vy = uy + ay.*t;
% Energies
     KE = (0.5*m) .* (vx.^2 + vy'.^2);
     PE = -(m.*ay) .* sy';
     E = KE + PE;
     
% GRAPHICS ===============================================================
 
 figure(1)   % Trajectroy  ----------------------------------------------
   fs = 14;
   set(gcf,'Units','Normalized');
   set(gcf,'Position',[0.02 0.04 0.18 0.33]);
   xP = sx; yP = sy;
   plot(xP,yP,'b','linewidth',2);
   xlabel('s_x  [ m ]');
   ylabel('s_y   [ m ]');
   grid on
   axis equal
   axis([0 2*6 -2 2*6])
   set(gca,'xTick',0:2:2*6);
   set(gca,'yTick',0:2:2*6);
   set(gca,'fontSize',fs);
   
figure(2)   % Time a v s  ----------------------------------------------
   fs = 14;
   set(gcf,'Units','Normalized');
   set(gcf,'Position',[0.35 0.04 0.3 0.75]);
   subplot(3,1,1)
     xP = t; yP = sx;
     plot(xP,yP,'r','linewidth',2);
     hold on
     xP = t; yP = sy;
     plot(xP,yP,'b','linewidth',2);
     xlabel('time t  [ s ]');
     ylabel('s   [ m ]');
     legend('X','Y');
     grid on
     set(gca,'fontSize',fs);   
     
   subplot(3,1,2)
     xP = t; yP = vx;
     plot(xP,yP,'r','linewidth',2);
     hold on
     xP = t; yP = vy;
     plot(xP,yP,'b','linewidth',2);
     xlabel('time t  [ s ]');
     ylabel('v   [ m/s ]');
     legend('X','Y');
     grid on
     set(gca,'fontSize',fs);   
     
   subplot(3,1,3)
     xP = t; yP = axG;
     plot(xP,yP,'r','linewidth',2);
     hold on
     xP = t; yP = ayG;
     plot(xP,yP,'b','linewidth',2);
     xlabel('time t  [ s ]');
     ylabel('a   [ m/s^2 ]');
     legend('X','Y');
     grid on
     set(gca,'yLim',[-10.5, 0.5]);
     set(gca,'fontSize',fs);   
 
 figure(3)   % Energies  ----------------------------------------------
   fs = 14;
   set(gcf,'Units','Normalized');
   set(gcf,'Position',[0.6 0.04 0.3 0.25]);
   xP = t; yP = KE;
   plot(xP,yP,'b','linewidth',2);
   hold on
   xP = t; yP = PE;
   plot(xP,yP,'r','linewidth',2);
   xP = t; yP = E;
   plot(xP,yP,'m','linewidth',2);
   xlabel('t  [ s ]');
   ylabel('energy   [ J ]');
   grid on
   legend('KE','PE','E');
   
   set(gca,'fontSize',fs);
        