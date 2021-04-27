% nsAE.m

% Izhikevich [2D] model with diffusion

clear
close all
clc


global   uF vF dt a b Iext 

tic

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
  a = -0.02; b = -1; c = -45; d = 0; 
  vPeak = 30;        

  dx = 0.5;
  dt = 0.001;
  N = 101;
  NT = 300;
  Dv = 10; 
  Du = 0.1;
  
% Stimulus current
  Iext = 80.*ones(N,N);
  
  
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

% Current values (1) / Next values (2)
  v1 = zeros(N,N); v2 = zeros(N,N);
  u2 = zeros(N,N);

% Single neurons  
  UT = zeros(3,NT);
  
% Initial values  
  v1 = -63 + 30.*cos(2.*pi.*x.*y./2);
  u1 = b.*v1; 
  
  
% UPDATE VARIABLES: v, u ==============================================

figure(1)
   pos = [0.1 0.1 0.30 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
      
for ct = 1:NT
  v2 = V(N, v1, u1, dt, Iext, vF);
  v1 = v2;
    
  u2 = U(N, v1, u1, dt, a, b, uF);  
  u1 = u2;

  % Reset v and u
  for cx = 1: N
      for cy = 1:N
        if v1(cx,cy) >= vPeak            % a spike is fired!
           v1(cx,cy) = c;                % membrane voltage reset
           
  %        Reset u: used u2 and not u1 for reset
           u1(cx,cy) = u2(cx,cy) + d;    % recovery variable update
        end
     % if v1(cx,cy) < -100; v1(cx,cy) = -60; end
      end
  end
     
 
 % Boundary conditions
      v1(:,1) = v1(:,2);
      v1(:,N) = v1(:,N-1);
      v1(1,:) = v1(2,:);
      v1(N,:) = v1(N-1,:);
      v1(1,1) = v1(2,2);
      v1(N,N) = v1(N-1,N-1);
      
      u1(:,1) = u1(:,2);
      u1(:,N) = u1(:,N-1);
      u1(1,:) = u1(2,:);
      u1(N,:) = u1(N-1,:);
      u1(1,1) = u1(2,2);
      u1(N,N) = u1(N-1,N-1);
      
 
  % Single neurons
     UT(1,ct) = v1(50,50);
     UT(2,ct) = v1(10,15);
     UT(3,ct) = v1(43,50);
  

  if rem(ct,1) == 0
    pcolor(x,y,v1)
    shading interp
    axis square
    Hcolorbar = colorbar;
    set(Hcolorbar,'ylim',[-100 30])
    zlim([-100 30])
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

function v2 = V(N, v1, u1, dt, Iext, vF)
       
    v2 = zeros(N,N);
    x = 2:N-1;
    y = 2:N-1;
    v2(x,y) = v1(x,y) + dt*(0.04*v1(x,y)^2 + 5*v1(x,y) + 140 - u1(x,y) + Iext(x,y));
    v2(x,y) = v2(x,y) + vF.*(v1(x+1,y) - 2.*v1(x,y) + v1(x-1,y));
    v2(x,y) = v2(x,y) + vF.*(v1(x,y+1) - 2.*v1(x,y) + v1(x,y-1));
  end

 function u2 = U(N, v1, u1, dt, a, b, uF)
   u2 = zeros(N,N);
   x = 2:N-1;
   y = 2:N-1;
    
   u2(x,y) = u1(x,y) + dt.*(a.*(b.*(v1(x,y)) - u1(x,y)));
   u2(x,y) = u2(x,y) + uF.*(u1(x+1,y) - 2.*u1(x,y) + u1(x-1,y));
   u2(x,y) = u2(x,y) + uF.*(u1(x,y+1) - 2.*u1(x,y) + u1(x,y-1));
 end

