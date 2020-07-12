% mec_chaos_01.m

% Runge-Kutta Solutions to Equations of Motion
%  Modelling the motion of a SIMPLE PENDULUM

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180813 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm


clear all
close all
clc

global k m b

% Setup for saving images (im)   =======================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_chaos_001.gif';
% Delay in seconds before displaying the next image  
   delay = 0.20;  
% Frame counter start
   nt = 1; 


% INPUTS [default values]  ==============================================

% mass of system   [1]
   m = 1;
% spring constant (global variable)  [pi^2] 
   k = pi^2;
% damping constant b
   b = 0;
% Time domain:  nT must be an ODD number [501 0 8]
   nT = 501;
   tMin = 0;
   tMax = 10;
   t = linspace(tMin,tMax,nT);
   h = t(2) - t(1);

% Initialise arrays and set initial conditions [10 0]   
   x = zeros(nT,1);
   v = zeros(nT,1);
   
   x(1) = -10;
   v(1) = 0;
      
% Set plot limits - will need to change for different models
%  Time plots  [-10 10]  [-40 40]
   xLimits = [-10 10];
   yLimits = [-40 40];
   
  
  
% CALCULATIONS ========================================================

% Runga-Kutta Solution of differential equation
   for c = 1 : nT-1
       [k1, k2, k3, k4] = coeff(t(c),h,x(c),v(c));
  
       x(c+1) = x(c) + h*(v(c) + (k1 + k2 +k3)/6);  
      
       v(c+1) = v(c) + (k1 + 2*k2 + 2*k3 + k4)/6;
   end
   
 % Acceleration  a
    a = gradient(v)/h;
 % Forces
    Fnet = m.*a;
    Fs   = -k.*x;
 % Kinetic Energy  K
    K = (0.5*m) .* v.^2;
 % Potential energy   U
   Us = 0.5.*k.*x.^2;
 % Total Energy
    E = K + Us; 
    
%     dU = zeros(nT,1);
%     dU(1) = 0.5*F(1)*(x(2)-x(1));
%     
%     for c = 2:nT-1
%       dU(c+1) = dU(c) + 0.5*(F(c)+F(c+1))*(x(c+1)-x(c));
%     end

    
% GRAPHICS ==============================================================
   
   fs = 14;
figure(1)   % t vs x v a energy------------------------------------------
   pos = [0.02 0.05 0.29 0.7];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
subplot(4,1,1)   
   xP = t; yP = x;
   plot(xP,yP,'b','lineWidth',2);
   title('displacement','fontweight','normal');
   xlabel('t')
   ylabel('x')
   grid on
   set(gca,'fontsize',fs)
subplot(4,1,2);
   xP = t; yP = v;
   plot(xP,yP,'b','lineWidth',2);
   title('velocity','fontweight','normal');
   xlabel('t');
   ylabel('v');
   grid on
   set(gca,'fontsize',fs)
subplot(4,1,3)
   xP = t; yP = a;
   plot(xP,yP,'b','lineWidth',2);
   title('acceleration','fontweight','normal');
   ylabel('a');
   xlabel('t');
   grid on   
   set(gca,'fontsize',fs)
 subplot(4,1,4)
   xP = t; yP = K;
   plot(xP,yP,'r','lineWidth',2);
   hold on
   xP = t; yP = Us;
   plot(xP,yP,'b','lineWidth',2);
   xP = t; yP = E;
   plot(xP,yP,'k','lineWidth',2);
   ylabel('energy');
   xlabel('t');
   grid on   
   set(gca,'fontsize',fs)
   legend('K','U','E','location','northoutside','orientation','horizontal'),

   
figure(2)   % x vs  v f energy------------------------------------------
   pos = [0.32 0.05 0.29 0.7];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
subplot(3,1,1)   
   xP = x; yP = v;
   plot(xP,yP,'b','lineWidth',2);
   hold on
   xP = x(1); yP = v(1);
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',8,'markerFaceColor','r','markerEdgeColor','r')
   title('phase plot','fontweight','normal');
   xlabel('x')
   ylabel('v')
   grid on
   set(gca,'fontsize',fs)
subplot(3,1,2);
   xP = x; yP = Fnet;
   plot(xP,yP,'r','lineWidth',2);
   hold on
   xP = x; yP = Fs;
   plot(xP,yP,'b','lineWidth',2);
   hold on
   legend('F_{net}','F_S','location','northoutside','orientation','horizontal'),
   xlabel('x');
   ylabel('F');
   grid on
   set(gca,'fontsize',fs)
subplot(3,1,3)
   xP = x; yP = K;
   plot(xP,yP,'r','lineWidth',2);
   hold on
   xP = x; yP = Us;
   plot(xP,yP,'b','lineWidth',2)
   xP = x; yP = E;
   plot(xP,yP,'k','lineWidth',2)
   ylabel('energy');
   xlabel('x');
   grid on   
   set(gca,'fontsize',fs)
   legend('K','U','E','location','northoutside','orientation','horizontal'),
 
   
    
figure(4)   % t vs x   animation --------------------------------------
   pos = [0.42 0.45 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   for c = 1:5:nT
       
      subplot('Position',[0.15 0.36 0.8 0.6])
      xP = x; yP = v;
      plot(xP,yP,'b','lineWidth',2);
      xP = x(c); yP = v(c);
      hold on
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','r','markerEdgeColor','r')
      xlabel('x')
      ylabel('v')
      grid on
      set(gca,'fontsize',fs)
      xlim(xLimits)
      ylim(yLimits)
      hold off
      
      subplot('Position',[0.15 0.1 0.8 0.1])
      xP = x(c); yP = 1;
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','r','markerEdgeColor','r')
      xlabel('t')
      ylabel('x')
      grid on
      set(gca,'fontsize',fs)
      axis off
      xlim([-10, 10])
      pause(0.2)
      
      if f_gif > 0 
         frame = getframe(4);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
         if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
         nt = nt+1;
       end
      
   end
   

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
  global k m b
% Simple Harmonic Motion 
%    y = (- k * x) / m;       
    
% Damped Harmonic Motion    
   y = (- k * x - b * v)/m;
end
