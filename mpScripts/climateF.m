clear
close all
clc


% global a b

a = 10; b = 1e-3;
tSpan = [0 5*3600];
s0 = [0 0 0 10];

K = [a b];

[t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 


% Displacement  [m]
   x = sol(:,1); y = sol(:,2);
% Velocity  [m/s]
   u = sol(:,3); v = sol(:,4);
% accelertion
   ax = b.*v;
   ay = a - b.*u;

 %  ay(ay < 0) = 0;

% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.29 0.69];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
 subplot(3,1,1)   
    xP = x./1e3; yP = y./1e3;   
    plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('y  [ km ]')
   xlabel('x  [ km]')
   set(gca,'fontsize',12)
   
subplot(3,1,2)   
   xP = t./3600; yP = u;   
   plot(xP,yP,'b','linewidth',2)
   hold on
   yP = v;   
   plot(xP,yP,'r','linewidth',2)
   grid on
   ylabel('v  [ m/s ]')
   xlabel('t  [ h]')
   set(gca,'fontsize',12)
      
subplot(3,1,3)   
   xP = t./3600; yP = ax;   
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
   

function sDot = EM(~,s,K)
%global a b
a = K(1); b = K(2);
sDot(1) = s(3);
sDot(2) = s(4);
sDot(3) = b*sDot(2);
%sDot(3) = -0.0050;
sDot(4) = a - b*sDot(1);
% if sDot(4) < 0; sDot(4) = 0; end
%sDot(4) = -0.05;
sDot = sDot';


end

