% ODE_004.m


close all
clear
clc

% ========================================================================
%  DIRECTION - SLOPE FIELDS     t y tt yy   S = dy/dt
%  ISOCLINES  dy/dt = yDot = constant    ts  ys
% ========================================================================

% Isoclines  slope dy/dt = m = constant
% Initial conditions for integral (solution) curve t(1)  y(1)
  m = 1.0;     
  t(1) = -1;       
  y(1) = 1.5803;    

% ========================================================================
% ODE FUNCTIONS 
  z = 1;
%   z = 1      S = y^3 - 3y - t 
%   z = 2      S = y^2 - 2t
%   z = 3      S = 1 + t - y
%   z = 4      S = 2y/t
%   z = 5      S = -y / ((t^2 + y^2)
%   z = 6      S = t - 2y
%   z = 7      S = y - y^2   AUTONOMOUS EQUATION
%   Z = 8      s = 3y - y^2   AUTONOMOUS EQUATION
% ========================================================================
% Direction Field - grid
  tMin = -4; tMax = 4;
  yMin = -4; yMax = 4; 
  N = 19;
  tG = linspace(-tMax,tMax,N);
  yG = linspace(-tMax,tMax,N);
  [tt,yy] = meshgrid(tG,yG);

  ns = 999; 
  ys = linspace(-tMax,tMax,ns); 
% ========================================================================

  switch z
      case 1
           S =  yy.^3 - 3.*yy - tt;
           ts = ys.^3 - 3.*ys - m;
           tn = ys.^3 - 3.*ys; yn = ys;
      case 2
           S = yy.^2 - 2.*tt;
           ts = (ys.^2 - m)./2;
           tn = ys.^2./2; yn = ys;
      case 3
           S = 1 + tt- yy;
           ts = m - 1 + ys; 
           tn = -1 + ys; yn = ys;
      case 4
          S = 2.*yy./tt;
          ts = 2.*ys./m;
          tn = 0.*ys; yn = ys;
      case 5
          S = -yy./(tt.^2 + yy.^2);
          % isoclines are circles
       %   m = m;
          theta = linspace(0,2*pi,299);
          R = 1/(abs(2*m));
          ts = R.*cos(theta);
          ys = R.*sin(theta) - 1/(2*m);
          tn = 0.*ys; yn = ys;
      case 6
          S = tt - 2.*yy;
          ts = m + 2.*ys;
          tn = 2.*ys; yn = ys;
      case 7
          S = yy.*(1-yy);
          ts = linspace(-tMax,tMax,ns); tn = ts;
         
          % solve quadratic
           s1 = (1 + sqrt(1-4*m))/2;
           s2 = (1 - sqrt(1-4*m))/2;
           ys = s1.*ones(ns,1); 
           ys2 = s2.*ones(ns,1); 

           yn = zeros(ns,1); 
           yn2 = ones(ns,1);
      case 8
           S = 3.*yy - yy.^2;
           ts = linspace(-tMax,tMax,ns); tn = ts;
           yn = 3.*ones(ns,1);
           yn2 = zeros(ns,1);

           % solve quadratic
           s1 = (-3 + sqrt(9 - 4*m))/2;
           s2 = (-3 - sqrt(9 - 4*m))/2;
           ys = s1.*ones(ns,1); 
           ys2 = s2.*ones(ns,1); 
  end


   L = sqrt(1+S.^2);


% ========================================================================
%  INTEGRAL OR SOLUTION CURVES
%      Solve ODE using Euler method
% ========================================================================
     dt = 1e-5;
     s = 0;
     c = 1;
     dy = 1e-4;
 while abs(s) < 4 && c < 1e6 % &&   dy > 1e-8
    switch z 
      case 1
        y(c+1) = y(c) + dt*( y(c).^3 - 3*y(c) - t(c) );
      case 2
        y(c+1) = y(c) + dt*( y(c).^2 - 2*t(c) );
      case 3
         y(c+1) = y(c) + dt*( 1 + t(c) - y(c) );
     case 4
         y(c+1) = y(c) + dt*( 2*y(c)/t(c) ); 
     case 5
         y(c+1) = y(c) + dt*( -y(c)/(t(c)^2 + y(c)^2) ); 
     case 6
         y(c+1) = y(c) + dt*( t(c) - 2*y(c) ); 
     case 7
         y(c+1) = y(c) + dt*( y(c)*(1-y(c)) ); 
        case 8
         y(c+1) = y(c) + dt*( 3*y(c) - y(c)^2 );    
    end
    
     dy = abs(y(c+1) - y(c));
     s = y(c+1);
     t(c+1) = t(c) + dt;
     c = c+1;
     
 end
     tLen = length(t);
 

% ========================================================================
% GRAPHICS 
% ========================================================================

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'Position', [0.1 0.1 0.3 0.35])
  set(gcf,'color','w');
  FS = 14;
  
% Slope (direction) Field
  q = quiver(tt,yy,1./L,S./L,0.6,'b');
  q. AutoScaleFactor = 1.2;
  hold on

% Isocline
  plot(ts,ys,'m','linewidth',2);
  plot(tn,yn,'r','linewidth',2);

  if z == 7 || z == 8
     plot(ts,ys2,'m','linewidth',2);
     plot(tn,yn2,'r','linewidth',2);  
  end

% Integral (solution curve)
  plot(t,y,'k','linewidth',2);
  Hp = plot(t(1),y(1),'go') ;
  set(Hp,'markerfacecolor','g','markersize',8)
  
  if max(abs(y)) < 4
    Hp = plot(t(end),y(end),'ro') ;
   set(Hp,'markerfacecolor','r','markersize',8)
  end

  box on
  grid on
  xlabel('t')
  ylabel('y')
  set(gca,'FontSize',FS)
  xlim([-4 4])
  ylim([-4 4])
  xticks(-4:1:4)
  yticks(-4:1:4)
  %axis tight


% ======================================================================
% PHASE PLOTS
%  ======================================================================
 
  y = linspace(-4,4,599);

  figure(2)
  set(gcf,'units','normalized');
  set(gcf,'Position', [0.5 0.1 0.3 0.35])
  set(gcf,'color','w');
  FS = 14;
  hold on
 switch z
     case 1
        t = -4:2:4;
       for c = 1 : length(t)
           yDot = y.^3 - 3.*y - t(c);
           xP = y; yP = yDot;
           plot(xP,yP,'LineWidth',1)
       end

     case 7
           yDot = y - y.^2;
           xP = y; yP = yDot;
           plot(xP,yP,'LineWidth',1)   
    case 8
           yDot = 3.*y - y.^2;
           xP = y; yP = yDot;
           plot(xP,yP,'LineWidth',1)
      

 
 end
   plot([-4 4],[0 0],'r','LineWidth',2)
   box on
   grid on
   xticks(-4:1:4)
   xlabel('y')
   ylabel('dy/dt')
   set(gca,'FontSize',FS)


