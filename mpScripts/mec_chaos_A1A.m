% mec_chaos_01.m

% Runge-Kutta Solutions to Equations of Motion
%  Modelling the motion of a SIMPLE PENDULUM

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180813 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm


clear all
close all
clc


% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_chaos_A001.gif';
% Delay in seconds before displaying the next image  
   delay = 0.20;  
% Frame counter start
   nt = 1; 

% VARIABLES   ===========================================================
   % angular displaecment (theta)   xA   [rad]
   % angular velocity     (omega)   vA   [rad/s]
   % anular acceleration  (alpha)   aA   [rad/s^2]
   

   
% INPUTS [default values]  ==============================================
% acceleration due to gravity g   [m.s^-2]
   g = 9.8;
% length of pendulum L   [m]
   L = 9.8;
% mass of pendulum bob   m   [1]
   m = 1.0;
   
% Time domain:  nT must be an ODD number [501 0 8]
   nT = 501;
   tMin = 0;
   tMax = 50;
   t = linspace(tMin,tMax,nT);
   h = t(2) - t(1);

% Initialise arrays and set initial conditions [-pi/6  -1.93185]   
   xA = zeros(nT,1);
   vA = zeros(nT,1);
   
   xA(1) = -pi/6;
   vA(1) = -1.93185;
        
  
% CALCULATIONS ========================================================

% Runga-Kutta Solution of differential equation
   for c = 1 : nT-1
       [k1, k2, k3, k4] = coeff(t(c),h,xA(c),vA(c));
  
       xA(c+1) = xA(c) + h*(vA(c) + (k1 + k2 +k3)/6);  
      
       vA(c+1) = vA(c) + (k1 + 2*k2 + 2*k3 + k4)/6;
   end
   
 % Acceleration  a
    aA = gradient(vA)/h;
 % Force (restoring force)  F [N]
    F   = -m*g*sin(xA);
 % Kinetic Energy  K   [J]
    K = (0.5*m*L^2) .* vA.^2;
 % Potential energy   U
    U = (m*g*L).*(1-cos(xA));
 % Total Energy
    E = K + U; 
 % XY coodinates of the pendulum
   X = L.*sin(xA);
   Y = -L.*cos(xA);
    
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
   xP = t;  yP = xA/pi;   %yP = mod(xA,2*pi)./pi;
   plot(xP,yP,'b','lineWidth',2);
   title('angular displacement','fontweight','normal');
   xlabel('t   [s ]')
   ylabel('\theta  [ \pi rad ]')
   grid on
   set(gca,'fontsize',fs)
subplot(4,1,2);
   xP = t; yP = vA;
   plot(xP,yP,'b','lineWidth',2);
   title('angular velocity','fontweight','normal');
   xlabel('t  [ s ]');
   ylabel('\omega  [ rad.s^{-1} ]');
   grid on
   set(gca,'fontsize',fs)
subplot(4,1,3)
   xP = t; yP = aA;
   plot(xP,yP,'b','lineWidth',2);
   title('angular acceleration','fontweight','normal');
   ylabel('\alpha  [ rad.s^{-2} ]');
   xlabel('t  [ s ]');
   grid on   
   set(gca,'fontsize',fs)
 subplot(4,1,4)
   xP = t; yP = K;
   plot(xP,yP,'r','lineWidth',2);
   hold on
   xP = t; yP = U;
   plot(xP,yP,'b','lineWidth',2);
   xP = t; yP = E;
   plot(xP,yP,'k','lineWidth',2);
   ylabel('energy  [ J ]');
   xlabel('t  [ s ]');
   grid on   
   set(gca,'fontsize',fs)
   legend('K','U','E','location','northoutside','orientation','horizontal'),

hold off

figure(2)   % x vs  v f energy------------------------------------------
   pos = [0.32 0.05 0.29 0.7];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
subplot(3,1,1)   
   xP = xA/pi; yP = vA;
   plot(xP,yP,'b','lineWidth',2);
   hold on
   xP = xA(1)/pi; yP = vA(1);
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',8,'markerFaceColor','r','markerEdgeColor','r')
   title('phase plot','fontweight','normal');
   xlabel('\theta  [ \pi rad ]')
   ylabel('\omega  [ rad.s^{-1} ]')
   grid on
   set(gca,'fontsize',fs)
subplot(3,1,2);
   xP = xA/pi; yP = F;
   plot(xP,yP,'b','lineWidth',2);
   hold on
   title('Restoring force','fontweight','normal');
   xlabel('\theta  [ \pi rad ]');
   ylabel('F  [ N ]');
   grid on
   set(gca,'fontsize',fs)
subplot(3,1,3)
   xP = xA/pi; yP = K;
   plot(xP,yP,'r','lineWidth',2);
   hold on
   xP = xA/pi; yP = U;
   plot(xP,yP,'b','lineWidth',2)
   xP = xA/pi; yP = E;
   plot(xP,yP,'k','lineWidth',2)
   ylabel('energy  [ J ]');
   xlabel('\theta  [ \pi rad ]');
   grid on   
   set(gca,'fontsize',fs)
   legend('K','U','E','location','northoutside','orientation','horizontal'),

hold off   


figure(4)   % t vs x   animation --------------------------------------
   pos = [0.62 0.05 0.30 0.54];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   for c = 1:5:nT
      
      subplot('position',[0.1 0.1 0.8 0.45]) 
      xP = xA./pi; yP = vA;
      plot(xP,yP,'b','lineWidth',2);
      xP = xA(c)/pi; yP = vA(c);
      hold on
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','m','markerEdgeColor','m')
      xP = xA(1)/pi; yP = vA(1);
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','r','markerEdgeColor','r')
      
      xlabel('\theta  [ \pi rad ]')
      ylabel('\omega  [ rad.s^{-1}')
      grid on
      set(gca,'fontsize',12)
      xlim([min(xA)/pi max(xA)/pi])
   %   set(gca,'xtick',-0.5:0.1:0.5)
      hold off
   
   subplot('position',[0.2 0.62 0.6 0.35])    
   xP = 0; yP = 0;
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',8,'markerFaceColor','k','markerEdgeColor','k')
   xP = [0 X(c)]; yP = [0 Y(c)];
   plot(xP,yP,'k','lineWidth',1);
   hold on
   xP = X(c); yP = Y(c);
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',12,'markerFaceColor','m','markerEdgeColor','m')
   xlim([-10 10])
   ylim([-10 10])
   grid on
   set(gca,'fontsize',fs)
%    xlabel('x  [ m ]')
%    ylabel('y  [ m ]')
   axis square   
   axis off
   hold off
   
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
hold off


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
  
   y = -sin(x);        

end
