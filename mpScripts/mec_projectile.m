clear
%close all
clc


% INPUTS ===============================================================
  g = 9.8;
  theta = 30;
  b = 0.5;
  x0 = 0; y0 = 0;
  v0 = 10;
  tSpan = [0 2];

% SETUP ===============================================================
  vx0 = v0*cosd(theta); vy0 = v0*sind(theta);
  s0 = [x0 y0 vx0 vy0];

  K(1) = g; K(2) = b*cosd(theta); K(3) = b*sind(theta);

  opts = odeset('event',@gg);

  [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0,opts); 
 

% Displacement  [m]
   x = sol(:,1); y = sol(:,2);
% Velocity  [m/s]
   vx = sol(:,3); vy = sol(:,4);
   v = sqrt(vx.^2 + vy.^2);
% accelertion
   ax = -K(2).*v.^2;
   ay = -K(1) - K(3).*v.^2;

   max(y)
     

% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.29 0.69];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
 subplot(3,1,1)   
    xP = x; yP = y;   
    plot(xP,yP,'b','linewidth',2)
    xlim([0 10])
   grid on
   hold on
   ylabel('y  [ km ]')
   xlabel('x  [ km]')
   set(gca,'fontsize',12)
   
subplot(3,1,2)   
   xP = t; yP = vx;   
   plot(xP,yP,'b','linewidth',2)
   hold on
   yP = vy;   
   plot(xP,yP,'r','linewidth',2)
   grid on
   ylabel('v  [ m/s ]')
   xlabel('t  [ h]')
   set(gca,'fontsize',12)
      
subplot(3,1,3)   
   xP = t; yP = ax;   
   plot(xP,yP,'b','linewidth',2)
   hold on
   yP = ay;
    plot(xP,yP,'r','linewidth',2)
%    yP = -1e-3.*ones(length(xP),1);
%     plot(xP,yP,'m','linewidth',2)
%     yP = gradient(u,t(2)-t(1));
%   plot(xP,yP,'k','linewidth',2)  
   grid on
   ylabel('a  [ m/s^2 ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',12)   
   

function sDot = EM(t,s,K)
 

   sDot(1) = s(3);
   sDot(2) = s(4);

 % v2 = sDot(1)^2 + sDot(2)^2;
  v2 = s(3)^2 + s(4)^2;

   sDot(3) = -K(2)*v2;
%   if sDot(3) < 0; sDot(3) = 0 ; end

   sDot(4) = -K(1) - K(3)*v2;

  sDot = sDot';


end


function [ggstop,isterminal,direction] = gg(t,s)
ggstop = s(2);
isterminal = 1;
direction = [];
end


