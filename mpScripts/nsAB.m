% ns20_001.m

clear
close all
clc

% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_A.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 
  

% CONSTANTS
  c(1) = 0.09; c(2) = 0.5; c(3) = 0.7;  c(4) = 1.0;
%  c(1) = 1.95;
  Iext = 1;
  
  dt = 0.01;
  dx = 0.5;
  
% Stability constant
  s = dt*c(4)/dx^2;
 % s = 0;

   % Time
  NT = 2000;
  t = 0:dt:dt*(NT-1);

% X position
  NX = 100;
  x = linspace(0,dx*NX,NX);
  
  u = zeros(NX,NT);
  v = zeros(NX,NT);
  
  u(:,1) = -2 + 2*rand(NX,1);
  v(:,1) = ( u(:,1) + c(3) ) / c(2);  
  

for q = 1:NT-1     
for p = 2:NX-1
    u(p,q+1) = u(p,q)   + dt*( u(p,q) - u(p,q)^3/3 - v(p,q) + Iext ) / c(1);
    u(p,q+1) = u(p,q+1) + s*(u(p-1,q) - 2*u(p,q) + u(p+1,q));
    v(p,q+1) = v(p,q)   + dt*(u(p,q+1) - c(2)*v(p,q) + c(3) );
end
    u(1,q+1) = u(1,q)   + dt*( u(1,q) - u(1,q)^3/3 - v(1,q) + Iext ) / c(1);
    u(1,q+1) = u(1,q)   + s*( u(NX-1,q) - 2*u(1,q) + u(3,q) )/4;
    v(1,q+1) = v(1,q)   + dt*(u(1,q+1) - c(2)*v(1,q) + c(3) );
    
    u(NX,q+1) = u(1,q+1);
    u(NX,q+1) = u(1,q+1);
    v(NX,q+1) = u(1,q+1);
end   


figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.25 0.25]);
  set(gcf,'color','w');

for q = 1:1:NT
   xP = x; yP = u(:,q)';
   plot(xP,yP,'b','linewidth',2)
   xlabel('x  [a.u.]')
   ylabel('v_{membrane}  [mV]')
   grid on
   set(gca,'fontsize',12)
   ylim([-2,2]);
   xlim([0 max(x)])
   pause(0.00000001)
   
   if rem(q,10) == 0
   if f_gif > 0 
         frame = getframe(1);
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
end


figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.32 0.05 0.25 0.25]);
   set(gcf,'color','w');
   xP = t; yP = u(52,:);
   plot(xP,yP,'b','linewidth',2)
   hold on
   xP = t; yP = u(10,:);
   plot(xP,yP,'r','linewidth',2)
    xP = t; yP = u(90,:);
   plot(xP,yP,'m','linewidth',2)
   xlabel('t  [ms]')
   ylabel('v_{membrane}  [mV]')
   grid on
   set(gca,'fontsize',12)
   

  % FUNCTIONS  ========================================================
%   function vM = V(N, v, u, dt, I, s)
%    vM = zeros(N,1);
%    x = 2:N-1;
%    vM(x) = v(x) + dt.*(0.04.*v(x).*v(x) + 5.*v(x) + 140 - u(x) + I) + s.*(v(x-1) - 2*v(x) + v(x+1));
%   end
% 
% function uM = U(N, v, u,dt, a, b)
%    uM = zeros(N,1);
%    x = 2:N-1;
%    uM(x) = u(x) + dt.*a.*(b*v(x) - u(x)) + s.*(u(x-1) - 2.*u(x) + u(x+1));
% end