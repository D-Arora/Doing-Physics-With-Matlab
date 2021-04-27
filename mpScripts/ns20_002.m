% ns20_001.m

clear
close all
clc

% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_A.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 
  

% CONSTANTS
% a = 0.02; b = 0.20; c = -65; d = 8; I = 10; D = 0.10;
 a = -0.02; b = -1; c = -60; d = 8; I = 80;  D = 0.10;
 s = 0.5;
  
% Time
  dt = 0.005;
  NT = 1000;
  t = 0:dt:dt*(NT-1);

  % X position
  dx = sqrt(dt*D/s);
  NX = 100;
  x = linspace(0,dx*NX,NX);
  
  v = c.*ones(NX,NT);
  v(:,1) = v(:,1) + rand(NX,1);
  u = b.*v;  
  
 % s = 0;
     
for p = 2: NX-1
for q = 1:NT-1
    v(p,q+1) = v(p,q) + dt*(0.04.*v(p,q)*v(p,q) + 5*v(p,q) + 140 - u(p,q) + I);
    v(p,q+1) = v(p,q+1) + s*(v(p-1,q) - 2*v(p,q) + v(p+1,q));
    
    u(p,q+1) = u(p,q) + dt*a*(b*v(p,q) - u(p,q));
    u(p,q)   = u(p,q+1) + s.*(u(p-1,q) - 2.*u(p,q) + u(p+1,q));
    
  if v(p,q+1) > 30
     v(p,q+1) = c;
     u(p,q+1) = u(p,q+1) + d;
  end
  
   v(NX,q+1) = v(NX-1,q+1);
   u(NX,q+1) = u(NX-1,q+1);
   v(1,q+1) = v(2,q+1);
   u(1,q+1) = u(2,q+1);
  
end
end

   
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.25 0.25]);
  set(gcf,'color','w');

for q = 1:1:NT
   xP = x; yP = v(:,q)';
   plot(xP,yP,'b','linewidth',2)
   xlabel('x  [a.u.]')
   ylabel('v_{membrane}  [mV]')
   grid on
   set(gca,'fontsize',12)
   ylim([-80,40]);
   xlim([0 max(x)])
   pause(0.00000001)
   
   if rem(q,5) == 0
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
   xP = t; yP = v(52,:);
   plot(xP,yP,'b','linewidth',2)
   hold on
   xP = t; yP = v(10,:);
   plot(xP,yP,'r','linewidth',2)
    xP = t; yP = v(90,:);
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