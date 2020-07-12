% chaos07A.m

% Runge-Kutta Solutions to Equations of Motion
%  Modelling the motion of a damped pendulum with an applied driving force

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180813 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/mec_chaosA.pdf

clear 
close all
clc

global c

% ANIMATION SETUP   ===================================================
% flagG = 0  animation not displayed / flagG = 1 animation displayed
   flagG = 1;
% f_gif = 0 animated gif NOT saved / f_gif = 1 file saved
   f_gif = 0;
% File name for animated gif   
   ag_name = 'ag_chaos_A001.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 

% CONSTANTS and DEFAULT VALUES ==========================================
   % angular displaecment (theta)   x   [rad]
   % angular velocity     (omega)   v   [rad/s]
   % anular acceleration  (alpha)   a   [rad/s^2]
 
% Differential equation constants
% Equation of motion coefficients  
%     a = -( c1*sin(x) + c2*v + c3*cos(c4*t) )
%     v = dx/dt
%     c(1) --> nautal frequency for free samll amplitude oscillations 
%     c(2) --> strength of damping
%     c(3) --> strength of external driving force (amplitude)
%     c(4) -->  angular frequenct of driving force w_ext 
%     c(2), c93) and c(4) values entered in INPUT section
      c = zeros(4,1);
% mass of pendulum bob   m  [kg]
   m = 1.0;
% acceleration due to gravity g   [m.s^-2]
   g = 9.8;
% nautral period and frequencies of oscillation: default value T_0 = 1 s
   T0 = 1; w0 = 2*pi/T0; f0 = 1/T0;
% length of pendulum L   [m]
   L = T0^2*g/(4*pi^2);
% defalut value for c(1)   
  c(1) = g/L;
   
% INPUTS  ============================================================

% Time domain 
   nT = 10000;
% Max time interval
   tMax = 25;
            tMin = 0;
            t = linspace(tMin,tMax,nT);
            h = t(2) - t(1);
% Time interval for phase space plot start tS and finsih tF             
% Iindices nS and nF for plotting figure 3  phase space plot  
   tS = 0;
   tF = t(end);
           nS = find(t >= tS,1); nF = find(t >= tF,1); nR = nS:nF;       
% Damping constant   
   c(2) =  0.2;
% Strength of driving force   
   c(3) =  15;
% Angular frequecny of driving force
   c(4) = (0.5)*w0;
% Initial position of pendulum  xA [rad]  vA [rad/s] 
             x = zeros(nT,1);
             v = zeros(nT,1); 
   x(1) = -pi/2;
   v(1) = 0.000;

%    k = 1.07875;
%    c(1) = 9*pi^2;
%    c(2) = 3*pi/2;
%    c(3) = k*c(1);
%    c(4) = 2*pi;
   
% CALCULATIONS ========================================================
% Driving frequency
    wExt = c(4); fExt = wExt/(2*pi); TExt = 1/fExt;
    
% time interval plotting range  
    XY = nS:nF;    

% Runga-Kutta Solution of differential equation
   for cc = 1 : nT-1
       [k1, k2, k3, k4] = coeff(t(cc),h,x(cc),v(cc));
  
       x(cc+1) = x(cc) + h*(v(cc) + (k1 + k2 +k3)/6);  
      
       v(cc+1) = v(cc) + (k1 + 2*k2 + 2*k3 + k4)/6;
   end
   
 % Acceleration  a
    a = gradient(v)/h;
 % Force (restoring force)  Fc [N]
    Fc   = -m*g*sin(x);
 % Kinetic Energy  K   [J]
    K = (0.5*m*L^2) .* v.^2;
 % Potential energy due to conservative force   Uc  [J]
    Uc = (m*g*L).*(1-cos(x));
 % Total Energy due to conservative forces [J]
    Ec = K + Uc; 
 % XY coodinates of the pendulum
    X = L.*sin(x);
    Y = -L.*cos(x);
 
 % find period of oscillation (from last two peaks)   
    [z1, z2] = findpeaks(x);
    z3 = length(z2);
    T = t(z2(z3)) - t(z2(z3-1))
    
    
    
% GRAPHICS ==============================================================
     
   fs = 12;
figure(1)   
   pos = [0.02 0.05 0.29 0.75];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = t(XY); 
   subplot(4,1,1)   
   yP = x(XY)/pi;   
   plot(xP,yP,'b','lineWidth',2);
   title('angular displacement','fontweight','normal');
   xlabel('t   [s ]')
   ylabel('\theta  [ \pi rad ]')
   grid on
   set(gca,'fontsize',fs)

   subplot(4,1,2);
   yP = v(XY);
   plot(xP,yP,'b','lineWidth',2);
   title('angular velocity','fontweight','normal');
   xlabel('t  [ s ]');
   ylabel('\omega  [ rad.s^{-1} ]');
   grid on
   set(gca,'fontsize',fs)

   subplot(4,1,3)
   yP = a(XY);
   plot(xP,yP,'b','lineWidth',2);
   title('angular acceleration','fontweight','normal');
   ylabel('\alpha  [ rad.s^{-2} ]');
   xlabel('t  [ s ]');
   grid on   
   set(gca,'fontsize',fs)
 
   subplot(4,1,4)
   yP = K(XY);
   plot(xP,yP,'r','lineWidth',2);
   hold on
   yP = Uc(XY);
   plot(xP,yP,'b','lineWidth',2);
   yP = Ec(XY);
   plot(xP,yP,'k','lineWidth',2);
   ylabel('energy  [ J ]');
   xlabel('t  [ s ]');
   grid on   
   set(gca,'fontsize',fs)
   legend('K','U_c','E_c','location','northoutside','orientation','horizontal'),

hold off

figure(2)   
   pos = [0.32 0.45 0.29 0.29];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
  
   xP = -pi:pi/200:+pi; yP = g*L*(1-cos(xP));
   yyaxis left
   plot(xP./pi,yP,'b','lineWidth',2);
   ylabel('V [ J.kg^{-1} ]')
   
   set(gca,'xcolor',[0 0 0])
   set(gca,'ycolor',[0 0 1])
   
   yyaxis right
   yP = -g*sin(xP);
   plot(xP./pi,yP,'r','lineWidth',2);
   set(gca,'ycolor',[1 0 0])
   ylabel('F_c / m [ N.kg^{-1} ]')
   hold on
   xP = [-1 1]; yP = [0 0];
   plot(xP,yP,'r','lineWidth',0.5);
   title('Conservative forces','fontweight','normal');
   xlabel('\theta [ \pi rad ]')
   grid on
   set(gca,'fontsize',fs)

hold off


figure(3)   
   pos = [0.32 0.05 0.29 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

            
   xP = x(XY)/pi; yP = v(XY);
   plot(xP,yP,'b','lineWidth',1);
   hold on
   xP = x(1)/pi; yP = v(1);
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',6,'markerFaceColor','r','markerEdgeColor','r')
   xP = x(end)/pi; yP = v(end);
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',6,'markerFaceColor','m','markerEdgeColor','m')
 % xlim([-1,5])
   tm1 = 'Phase Space:   t_S = ';
   tm2 = num2str(tS,'%3.1f \n');
   tm3 = ' s     t_F = ';
   tm4 = num2str(tF,'%3.1f \n');
   tm5 = ' s';
   tm = [tm1 tm2 tm3 tm4 tm5];
   title(tm,'fontweight','normal');
   xlabel('\theta  [ \pi rad ]')
   ylabel('\omega  [ rad.s^{-1} ]')
   grid on
   set(gca,'fontsize',fs)
hold off   

if flagG > 0
    
    
figure(4)   
   pos = [0.62 0.05 0.30 0.54];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   for cc = 1:20:nT
      
      subplot('position',[0.1 0.1 0.8 0.45]) 
      xP = x./pi; yP = v;
      plot(xP,yP,'b','lineWidth',2);
      xP = x(cc)/pi; yP = v(cc);
      hold on
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','m','markerEdgeColor','m')
      xP = x(1)/pi; yP = v(1);
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','r','markerEdgeColor','r')
      tm1 = 't =  ';
      tm3 = '  s';
      tm2 = num2str(t(cc),'%2.1f\n');
      tm = [tm1 tm2 tm3];
      title(tm,'fontweight','normal')
      xlabel('\theta  [ \pi rad ]')
      ylabel('\omega  [ rad.s^{-1}')
      grid on
      set(gca,'fontsize',12)
      xlim([min(x)/pi max(x)/pi])
   %   set(gca,'xtick',-0.5:0.1:0.5)
      hold off
   
   subplot('position',[0.2 0.62 0.6 0.35])    
   xP = 0; yP = 0;
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',8,'markerFaceColor','k','markerEdgeColor','k')
   xP = [0 X(cc)]; yP = [0 Y(cc)];
   plot(xP,yP,'k','lineWidth',2);
   hold on
   xP = X(cc); yP = Y(cc);
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',12,'markerFaceColor','m','markerEdgeColor','m')
   xlim([-1.1*L 1.1*L])
   ylim([-1.1*L 1.1*L])
   grid on
   set(gca,'fontsize',fs)
   axis square   
   axis off
   hold off
   
  pause(0.1)
      
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

end

% FUNCTIONS   =========================================================

% Runga-Kutta coefficients
function [k1, k2, k3, k4] = coeff(t,h,x,v)
  k1 = h*fn(t,x,v);
  k2 = h*fn(t+h/2, x+h*v/2, v+k1/2);
  k3 = h*fn(t+h/2, x+h*v/2+h*k1/4 ,v+k2/2);
  k4 = h*fn(t+h,   x+h*v+h*k2/2,   v+k3);  
end

% Equation of motion
function  y = fn(t,x,v)
  global c
  y = -c(1)*sin(x) - c(2)*v + c(3)*cos(c(4)*t) ; 
end
