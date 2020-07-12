%chaos02.m

% Runge-Kutta Solutions to Equations of Motion
%  Modelling Duffing oscillators: free, viscous damping, and forced motions
% POINCASRE SECTIONS

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180827 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%     ../mphome.htm
%   Reference documentation and notes
%     http://www.physics.usyd.edu.au/teach_res/mp/doc/chaos02.htm

clear 
close all
clc

tic
   
% MAIN VARIABLES =======================================================
   % S.I. units unless stated otherwise
   % displacement   x   
   % velocity       v   
   % acceleration   a   
   % equation of motion coefficients  c (global)
   
global c
   
% INPUTS ============================================================

% Equation of motion coefficients  [0.1,1,-1,,0.38,1.4]
%     a = c(1)*v + c(2)*x + c(3)*x^3 * c(4)*cos(c(5)*t)
%     c(1) damping coefficient
%     c(2) x coefficient / c(3) x^3 coefficent
%     c(4) amplitude / c(5) angular velocity of driving force 
   
   c(1) =  0.1;
   c(4) =  0.3425;
   c(5) =  1.4;
   
   c(2) =  1;
   c(3) = -1;
   
% mass of system  [1] 
   m = 1.0;
   
% Time domain
%    nT number of points before plotting
%    nP number of period iteriations / nS start plotting number
%    w_ext, T_ext drivingg force angular velocity & period
%    h time step
   nT = 501; 
   nP = 24000;
   nS = 1;
      w_ext = c(5);
      T_ext = 2*pi/w_ext; 
      tMin = 0;
      tMax = T_ext;
      t = linspace(tMin,tMax,nT);
      h = t(2) - t(1);
      
% Initialise arrays 
      x = zeros(nT,1);
      v = zeros(nT,1);
 
% Initial Conditions    
   x(1) = -1;
   v(1) = 1;
   

% OUPUT  ==============================================================   
        
figure(1)   % Poincare Section
   pos = [0.02 0.05 0.42 0.42];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xlim([-2 2])
   title('Poincare Section','fontweight','normal');
   xlabel('x  [ m ]');
   ylabel('v  [ m.s^{-1} ]');
   hold on
   box on
   set(gca,'fontsize',14)
   grid on
   
   
% CALCULATIONS ========================================================

for cP = 1:nP 
% Runga-Kutta Solution of differential equation
   for cc = 1 : nT-1
       [k1, k2, k3, k4] = coeff(t(cc),h,x(cc),v(cc));
  
       x(cc+1) = x(cc) + h*(v(cc) + (k1 + k2 +k3)/6);  
      
       v(cc+1) = v(cc) + (k1 + 2*k2 + 2*k3 + k4)/6;
   end
      xP = x(cc);  yP = v(cc);
      if cP > nS
        plot(xP,yP,'b.')
        xlim([-2 2])
        ylim([-1.5 1.5])
      end
      t = linspace(tMin,tMax,nT) + cP*T_ext;
      x(1) = x(cc);
      v(1) = v(cc);
end
    
toc


% FUNCTIONS   =========================================================

% Runga-Kutta coefficients
function [k1, k2, k3, k4] = coeff(t,h,x,v)
  k1 = h*fn(t,x,v);
  k2 = h*fn(t+h/2, x+h*v/2, v+k1/2);
  k3 = h*fn(t+h/2, x+h*v/2+h*k1/4 ,v+k2/2);
  k4 = h*fn(t+h,   x+h*v+h*k2/2,   v+k3);  
end

% Differential equation for acceleration (force / mass)

function  y = fn(t,x,v)
   global c
   y = -c(1)*v + c(2)*x + c(3)* x^3 + c(4)*cos(c(5)*t);

end
