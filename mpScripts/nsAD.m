% nsAC.m

clear
close all
clc


global   uF vF dt C k vr vt a b Iext 

tic

% INPUTS ==============================================================
 
  C = 100; vr = -60; vt = -40; k = 0.7; 
  a = 0.02; b = -1; c = -60; d = 8; 
  vPeak = 30;        

  dx = 0.5;
  dt = 0.01;
  N = 101;
  NT = 30;
  Dv = 10; 
  Du = 0.1;
  
  
% SETUP =============================================================== 
% Grid
  xG = linspace(0,dx*N,N);
  yG = xG;
  [x, y] = meshgrid(xG,xG);

% Time
  t = 0:dt:dt*(NT-1);
  
% Constants for 2nd derivates   constants < 1 for stability
  uF = (dt/(dx)^2)*Du;
  vF = (dt/(dx)^2)*Dv;

% Arrays
 
  
%  Iext = zeros(N,N);
   Iext = 84.9*ones(N,N);
%  Iext(38:68, 38:68) = 300 + 100.*rand(31,31); 
%  Iext(38:48, 45:55) = 80 + 40.*rand(11,11); 
%  Iext(52:62, 45:55) = 80 + 30.*rand(11,11);
%  Iext(52:52, 45:55) =  82;
%   
% Iext = 281.*ones(N,N) + 5.*rand(N,N);
  
% Potentials
% Current values / Next values / Initial values / boundary values
  v1 = zeros(N,N); 
  v2 = zeros(N,N); u2 = zeros(N,N);
  UT = zeros(3,NT);
  
  
  for nx = 1:N
      for ny = 1: N
        v1(nx,ny) = -63 + 30.*cos(2*pi*xG(nx)*yG(ny));
      end
  end
       u1 = b.*v1; 
  
% % Boundary Conditions  
%   v1(1,1) = v1(2,2);
%   v1(1,N) = v1(2,N-1);
%   v1(N,1) = v1(N-1,2);
%   v1(N,N) = v1(N-1,N-1);
%   
%   u1(1,1) = u1(2,2);
%   u1(1,N) = u1(2,N-1);
%   u1(N,1) = u1(N-1,2);
%   u1(N,N) = u1(N-1,N-1);
%     
%   v1(1,1) = -60; v(1,N) = -60;
%   v1(N,N) = -60; v(N,1) = -60;

% Update u and v variables

figure(1)
   pos = [0.1 0.1 0.30 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
      
for ct = 1:NT
  v2 = V(N, v1, u1, dt, vF, Iext, C, k, vr, vt);
  v1 = v2;
    
  u2 = U(N, v1, u1, dt, uF, a, b, vr);  
  u1 = u2;
  
  for cx = 1: N
      for cy = 1:N
        if v1(cx,cy) >= vPeak            % a spike is fired!
           v1(cx,cy) = c;                % membrane voltage reset
           u1(cx,cy) = u1(cx,cy) + d;    % recovery variable update
        end
       if v1(cx,cy) < -100; v1(cx,cy) = -60; end
      end
  end
  
  
  
 % for c1 = 1:N
      v1(:,1) = v1(:,2);
      u1(:,1) = u1(:,2);
      v1(:,N) = v1(:,N-1);
      u1(:,N) = u1(:,N-1);
      
      v1(1,:) = v1(2,:);
      u1(1,:) = u1(2,:);
      
%       u1(c1,1) = u1(c1,2);
%       u1(N,c1) = u1(N-1,c1);
%       u1(c1,N) = u1(c1,N-1);
%       v1(1,c1) = v1(2,c1);
%       v1(c1,1) = v1(c1,2);
%       v1(N,c1) = v1(N-1,c1);
%       v1(c1,N) = v1(c1,N-1);
%  end
 
  % Single neurons
     UT(1,ct) = v1(50,50);
     UT(2,ct) = v1(10,15);
     UT(3,ct) = v1(43,50);
  

  if rem(ct,1) == 0
    pcolor(x,y,v1)
  %   imagesc(xG,xG,u1)
    shading interp
    axis square
%     Hcolorbar = colorbar;
%     set(Hcolorbar,'ylim',[-60 40])
%    zlim([-60 40])
    txt = sprintf('t_{step} = %4.0f   max %2.1f   min  %2.1f', ct, max(max(v1)), min(min(v1)));
    title(txt)
    set(gca,'fontsize',12)
    set(gca,'YDir','normal')
    box on
    pause(0.1)
  end
end  

 
figure(2)
   pos = [0.5 0.1 0.25 0.25];
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
   
   
 toc  

% FUNCTIONS  ========================================================
  
  function v2 = V(N, v1, u1, dt, vF, Iext, C, k, vr, vt)
       
    v2 = zeros(N,N);
    x = 2:N-1;
    y = 2:N-1;
    v2(x,y) = v1(x,y) + (dt/2).*(k.*(v1(x,y) - vr).*(v1(x,y) - vt) - u1(x,y) + Iext(x,y))./C;
    v2(x,y) = v1(x,y) + (dt/2).*(k.*(v1(x,y) - vr).*(v1(x,y) - vt) - u1(x,y) + Iext(x,y))./C;
    v2(x,y) = v2(x,y) + vF.*(v1(x+1,y) - 2.*v1(x,y) + v1(x-1,y));
    v2(x,y) = v2(x,y) + vF.*(v1(x,y+1) - 2.*v1(x,y) + v1(x,y-1));
  end

 function u2 = U(N, v1, u1, dt, uF, a, b, vr)
   u2 = zeros(N,N);
   x = 2:N-1;
   y = 2:N-1;
    
   u2(x,y) = u1(x,y) + dt.*(a.*(b.*(v1(x,y) - vr) - u1(x,y)));
   u2(x,y) = u2(x,y) + uF.*(u1(x+1,y) - 2.*u1(x,y) + u1(x-1,y));
   u2(x,y) = u2(x,y) + uF.*(u1(x,y+1) - 2.*u1(x,y) + u1(x,y-1));
 end

