% nsAC.m

clear
close all
clc


global   uF vF dt Iext e g d


% INPUTS ==============================================================
  I0 = 1.6;  % External stimulus (not in paper)
  dx = 0.5;
  dt = 0.01;
  N = 101;
  NT = 3000;
  Du = 0.5; 
  Dv = 0;
  e = 0.1;  % epsilon
  g = 0.5;   % gamma
  d = 0.7;   % delta
  
% SETUP =============================================================== 
% Grid
  xG = linspace(0,dx*N,N);
  [x, y] = meshgrid(xG,xG);

% Time
  t = 0:dt:dt*(NT-1);
  
% Constants for 2nd derivates   constants < 1 for stability
  uF = (dt/(dx)^2)*Du;
  vF = (dt/(dx)^2)*Dv;

% Arrays
  u1 = zeros(N,N);
  u2 = u1; v1 = u1; v2 = v1;
  UT = zeros(3,NT);

  
% >>>> Change external stimjulus  
  Iext = zeros(N,N);
  Iext(38:48, 45:55) = -1 + 4.*rand(11,11); 
  Iext(52:62, 45:55) = -4 + 2.*rand(11,11);
  Iext(52:52, 45:55) =  6;
%  Iext(50:50, 45:55) =  -5.6;
% Iext(1:N,1:10) = 538;
% Iext(15:20,25:30) = 3;
% Iext(65:70,65:70) = 2.5;
% Iext = 1 + 0.8.*rand(N,N);
% >>>>>>>>>> 

% Potentials
% Current values / Next values / Initial values / boundary values
  u1 = -2 + 4.*rand(N,N);
  v1 = (u1 + d) ./ g;

% Boundary Conditions  
  u1(1,1) = u1(2,2);
  u1(1,N) = u1(2,N-1);
  u1(N,1) = u1(N-1,2);
  u1(N,N) = u1(N-1,N-1);
    
  v1(1,1) = v1(2,2);
  v1(1,N) = v1(2,N-1);
  v1(N,1) = v1(N-1,2);
  v1(N,N) = v1(N-1,N-1);
  

% Update u and v variables

figure(1)
   pos = [0.1 0.1 0.30 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
for k = 1:NT
  u2 = U(N, u1, v1, dt, uF, Iext, e);
  u1 = u2;
  
  v2 = V(N, u1, v1, dt, vF, g, d);  
  v1 = v2;
  
  for c1 = 2:N-1
      u1(1,c1) = u1(2,c1);
      u1(c1,1) = u1(c1,2);
      u1(N,c1) = u1(N-1,c1);
      u1(c1,N) = u1(c1,N-1);
      v1(1,c1) = v1(2,c1);
      v1(c1,1) = v1(c1,2);
      v1(N,c1) = v1(N-1,c1);
      v1(c1,N) = v1(c1,N-1);
  end
 
  % Single neurons
     UT(1,k) = u1(30,30);
     UT(2,k) = u1(50,50);
     UT(3,k) = u1(43,43);
  

  if rem(k,5) == 0
    pcolor(x,y,u1)
    shading interp
    axis square
    Hcolorbar = colorbar;
  %  set(Hcolorbar,'ylim',[-2.5 2.5])
    zlim([-2.5 2.5])
    txt = sprintf('t = %2.0f   max %2.1f   min  %2.1f', k, max(max(u1)), min(min(u1)));
    title(txt)
    set(gca,'fontsize',12)
    pause(0.01)
  end
end  


figure(2)
   pos = [0.4 0.1 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
   plot(t, UT(1,:),'linewidth',2)
   hold on
   plot(t, UT(2,:),'r')
   plot(t, UT(3,:),'m','linewidth',2)
   grid on
   xlabel('t');
   ylabel('v_{membrane}');
   set(gca,'fontsize',12)

  % FUNCTIONS  ========================================================
  function u2 = U(N, u1, v1, dt, uF, Iext, e)
   u2 = zeros(N,N);
   x = 2:N-1;
   y = 2:N-1;
   u2(x,y) = u1(x,y) + dt.*( (u1(x,y) - u1(x,y).^3./3 - v1(x,y))./e + Iext(x,y) );
   u2(x,y) = u2(x,y) + uF.*(u1(x+1,y) - 2*u1(x,y) + u1(x-1,y));
   u2(x,y) = u2(x,y) + uF.*(u1(x,y+1) - 2*u1(x,y) + u1(x,y-1));
  end

function v2 = V(N, u1, v1, dt, vF, g, d)
  v2 = zeros(N,N);
  x = 2:N-1;
  y = 2:N-1;
    
  v2(x,y) = v1(x,y) + dt.*( u1(x,y) - g.*v1(x,y) + d);
  v2(x,y) = v2(x,y) + vF.*(v1(x+1,y) - 2.*v1(x,y) + v1(x-1,y));
  v2(x,y) = v2(x,y) + vF.*(v1(x,y+1) - 2*v1(x,y) + v1(x,y-1));
end