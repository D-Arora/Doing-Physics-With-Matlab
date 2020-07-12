% chaos07AB.m

% Runge-Kutta Solutions to Equations of Motion
%  Modelling the motion of a damped pendulum with an applied driving force
% BIFURCATION DIAGRAM
% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180813 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/mec_chaosA.pdf

clear 
close all
clc

global c

tic

% CONSTANTS and DEFAULT VALUES ==========================================
   % angular displaecment (theta)   x   [rad]
   % angular velocity     (omega)   v   [rad/s]
   
 
% Differential equation constants
% Equation of motion coefficients  
%     a = -( c1*sin(x) + c2*v + c3*cos(c4*t) )
%     v = dx/dt
%     c(1) --> nautal frequency for free samll amplitude oscillations 
%     c(2) --> strength of damping
%     c(3) --> strength of external driving force (amplitude)
%     c(4) -->  angular frequenct of driving force w_ext 
%     c(2), c93) and c(4) values entered in INPUT section
      c = zeros(4,1);
   
% INPUTS ===========================================================

% Time domain 
   nT = 5000;
% Max time interval
   tMax = 50;
            tMin = 0;
            t = linspace(tMin,tMax,nT);
            h = t(2) - t(1);
% Time interval for phase space plot start tS and finsih tF             
% Iindices nS and nF for plotting figure 3  phase space plot  
   tS = 0;
   tF = t(end);
           nS = find(t >= tS,1); nF = find(t >= tF,1); nR = nS:nF;       
% nautral frequency
   c(1) = 9*pi^2;         
% Damping constant   
   c(2) =  3*pi/2;
% Angular frequecny of driving force
   c(4) = 2*pi;
% Initial position of pendulum  xA [rad]  vA [rad/s] 
             x = zeros(nT,1);
             v = zeros(nT,1); 
   x(1) = 0;
   v(1) = 0.000;

% Strength of driving force c(3) = Fs*c1
% nS       number of calculations for varying strength of driving force
% Fs      strength of driving force factor
% theta   angular displacent at the end of the simulation
%         -pi <= theta <= +pi
   nS = 2000;
   Fs = linspace(0.6,1.8,nS);
   theta = zeros(nS,1);

% CALCULATIONS ========================================================

for cs = 1 : nS
    c(3) = Fs(cs)*c(1);
% Runga-Kutta Solution of differential equation
   for cc = 1 : nT-1
       [k1, k2, k3, k4] = coeff(t(cc),h,x(cc),v(cc));
  
       x(cc+1) = x(cc) + h*(v(cc) + (k1 + k2 +k3)/6);  
      
       v(cc+1) = v(cc) + (k1 + 2*k2 + 2*k3 + k4)/6;
   end
   
   theta(cs) = x(end);
   while theta(cs) > +pi; theta(cs) = theta(cs) - pi; end
   while theta(cs) < -pi; theta(cs) = theta(cs) + pi; end
end 
  
    
    
% GRAPHICS ==============================================================
     
   fs = 14;
figure(1)   
   pos = [0.32 0.05 0.29 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
            
   xP = Fs; yP = theta./pi;
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',2,'markerFaceColor','b','markerEdgeColor','b')
  
   ylim([-1.1,1.1])
   xlabel('strength (F_{driving})')
   ylabel('\theta  [ rad ]')
   grid on
   set(gca,'fontsize',fs)

toc


% FUNCTIONS   =========================================================

% Runga-Kutta coefficients
function [k1, k2, k3, k4] = coeff(t,h,x,v)
  k1 = h*fn(t,x,v);
  k2 = h*fn(t+h/2, x+h*v/2, v+k1/2);
  k3 = h*fn(t+h/2, x+h*v/2+h*k1/4 ,v+k2/2);
  k4 = h*fn(t+h,   x+h*v+h*k2/2,   v+k3);  
end

% Equation of motion
function  y = fn(t,x,v)
  global c
  y = -c(1)*sin(x) - c(2)*v + c(3)*cos(c(4)*t) ; 
end
