%chaos01.m

% Runge-Kutta Solutions to Equations of Motion
%  Modelling Duffing oscillators: free, viscous damping, and forced motions

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180823 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

clear all
close all
clc


% ANIMATION SETUP ======================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
% To plot figure 4 flag4 = 1 otherwise flag4 = 0
  flag4 = 1;
  f_gif = 1;
   ag_name = 'ag_chaos_A001.gif';
% Delay in seconds before displaying the next image  
   delay = 0.10;  
% Frame counter start
   nt = 1; 




   
% MAIN VARIABLES =======================================================
   % S.I. units unless stated otherwise
   % displacement   x   
   % velocity       v   
   % acceleration   a   
   % equation of motion coefficients  c (global)
   
global c
   
% INPUTS [default values] ==============================================

% Equation of motion coefficients  [0,1,-1,,0,0]
%     a = c(1)*v + c(2)*x + c(3)*x^3 * c(4)*cos(c(5)*t)
%     c(1) damping coefficient
%     c(2) x coefficient / c(3) x^3 coefficent
%     c(4) amplitude / c(5) angular velocity of driving force 
   
   c(1) =  0.1;
   c(4) =  0.38;
   c(5) =  1.4;
   
   c(2) =  1;
   c(3) = -1;
   
% mass of system  [1] 
   m = 1.0;
   
% Time domain 
% Time interval for plots  start tS (index nS) / finish  tF (index nF)
% domain for plttting nR
   nT = 18001;
   tMin = 0;
   tMax = 3600;
   
   t = linspace(tMin,tMax,nT);
   h = t(2) - t(1);
   
   tS = 3300;
   tF = t(end);
   
   nS = find(t >= tS,1); nF = find(t >= tF,1); nR = nS:nF;
   
% Initialise arrays 
   x = zeros(nT,1);
   v = zeros(nT,1);
 
% Iitial Conditions   [x(1) = 1 v(1) = 1]   
   x(1) = 0;
   v(1) = 0;
        


   
% CALCULATIONS ========================================================

% Runga-Kutta Solution of differential equation
   for cc = 1 : nT-1
       [k1, k2, k3, k4] = coeff(t(cc),h,x(cc),v(cc));
  
       x(cc+1) = x(cc) + h*(v(cc) + (k1 + k2 +k3)/6);  
      
       v(cc+1) = v(cc) + (k1 + 2*k2 + 2*k3 + k4)/6;
   end
   
 % Acceleration  a   
    a = c(1).*v + c(2).*x + c(3).*x.^3 + c(4).*cos(c(5)*t)';
  % Net Force Fnet  
    Fnet = m*a;
 % Kinetic Energy  K
    K = 0.5*m*v.^2;
 % Potential energy for conservation forces Fc: x = 0 --> Uc = 0
   Fc =  m.*(x - x.^3);
   Uc = m*(x.^4./4 - x.^2./2);
 % Potential V   (potential energy per unit mass)     
   V = Uc./m;
 % Total Energy E 
    E = K + Uc; 
 %driving force
   if c(5) > 0
       w_ext = c(5);
       T_ext = 2*pi/w_ext;
       f_ext = 1/T_ext;
   end
    
% GRAPHICS ==============================================================
     
   fs = 14;
figure(1)   % time graphs ----------------------------------------------
   pos = [0.02 0.05 0.29 0.75];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
subplot(4,1,1)   
   xP = t(nR);  yP = x(nR);   
   plot(xP,yP,'b','lineWidth',2);
   title('displacement','fontweight','normal');
   xlabel('t   [s ]')
   ylabel('x  [ m ]')
   grid on
   set(gca,'fontsize',fs)
subplot(4,1,2);
   xP = t(nR); yP = v(nR);
   plot(xP,yP,'b','lineWidth',2);
   title('velocity','fontweight','normal');
   xlabel('t  [ s ]');
   ylabel('v  [ m.s^{-1} ]');
   grid on
   set(gca,'fontsize',fs)
subplot(4,1,3)
   xP = t(nR); yP = a(nR);
   plot(xP,yP,'b','lineWidth',2);
   title('acceleration','fontweight','normal');
   ylabel('a  (F_{net} / m)  [m.s^{-2}]');
   xlabel('t  [ s ]');
   grid on   
   set(gca,'fontsize',fs)
  subplot(4,1,4)
   xP = t(nR); yP = K(nR);
   plot(xP,yP,'r','lineWidth',2);
   hold on
   xP = t(nR); yP = Uc(nR);
   plot(xP,yP,'b','lineWidth',2);
   xP = t; yP = E;
   plot(xP(nR),yP(nR),'k','lineWidth',2);
   ylabel('energy  [ J ]');
   xlabel('t  [ s ]');
   grid on   
   set(gca,'fontsize',fs)
   legend('K','U_c','E','location','northoutside','orientation','horizontal'),

hold off

  figure(2)   % displacement graphs ---------------------------------------
   pos = [0.32 0.45 0.29 0.29];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
  
   xP = -2:0.10:2; yP = -(c(3).*xP.^4./4 + c(2).*xP.^2./2);
   yyaxis left
   plot(xP,yP,'b','lineWidth',2);
   ylabel('V [ J.kg^{-1} ]')
   
   set(gca,'xcolor',[0 0 0])
   set(gca,'ycolor',[0 0 1])
   
   yyaxis right
   yP = (c(3).*xP.^3 + c(2).*xP);
   plot(xP,yP,'r','lineWidth',2);
   set(gca,'ycolor',[1 0 0])
   ylabel('F_c / m [ N.kg^{-1} ]')
   hold on
   xP = [-2 2]; yP = [0 0];
   plot(xP,yP,'r','lineWidth',0.5);
   title('Conservative forces','fontweight','normal');
   xlabel('x [ m ]')
   grid on
   set(gca,'fontsize',fs)

hold off

figure(3)   % displacement graphs ---------------------------------------
   pos = [0.32 0.05 0.29 0.29];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
  
   xP = x(nR); yP = v(nR);
   plot(xP,yP,'b','lineWidth',1);
  % plot(xP,yP,'b.');
   hold on
   xP = x(1); yP = v(1);     % initial values
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',6,'markerFaceColor','r','markerEdgeColor','r')
   xP = x(end); yP = v(end); % final values
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',6,'markerFaceColor','m','markerEdgeColor','m')
   xlim([-2 2])
   tm1 = 'Phase Space:   t_S = ';
   tm2 = num2str(tS,'%3.0f \n');
   tm3 = ' s     t_F = ';
   tm4 = num2str(tF,'%3.0f \n');
   tm5 = ' s';
   tm = [tm1 tm2 tm3 tm4 tm5];
   title(tm,'fontweight','normal');
   xlabel('x  [ m ]');
   ylabel('v  [ m.s^{-1} ]');
   grid on
   set(gca,'fontsize',fs)

   
hold off   


if flag4 == 1
figure(4)   % t vs x   animation --------------------------------------
   pos = [0.62 0.05 0.30 0.54];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   for cc =  nS:3:nF  %nR %1:10:nT
      
      subplot('position',[0.15 0.6 0.7 0.35]) 
      xP = -2:0.10:2; yP = (xP.^4./4 - xP.^2./2);
      plot(xP,yP,'b','lineWidth',2);
      xP = x(cc); yP = V(cc);
      hold on
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','m','markerEdgeColor','m')
      xP = x(1); yP = V(1);
      Hplot = plot(xP,yP,'ro');
      set(Hplot,'markersize',12,'markerFaceColor','r','markerEdgeColor','r')
       xlabel('x  [ m ]')
      ylabel('V  [ J.kg^{-1} ]')
      grid on
      set(gca,'fontsize',12)
      xlim([-2 2])
      ylim([-0.5 2])
   %   set(gca,'xtick',-0.5:0.1:0.5)
      hold off
   
     subplot('position',[0.15 0.1 0.7 0.35])    
     xP = x; yP = v;
     plot(xP(nR),yP(nR),'b','linewidth',1);
     hold on
     xP = x(cc); yP = v(cc);
     Hplot = plot(xP,yP,'o');
     set(Hplot,'markersize',12,'markerFaceColor','m','markerEdgeColor','m')
     xP = x(1); yP = v(1);
     Hplot = plot(xP,yP,'o');
     set(Hplot,'markersize',12,'markerFaceColor','r','markerEdgeColor','r')
     xlim([-2 2])
     grid on
     set(gca,'fontsize',fs)
     xlabel('x  [ m ]')
     ylabel('v  [ m.s^{-1} ]')
     tm1 = 't =  ';
     tm3 = '  s';
     tm2 = num2str(t(cc),'%2.1f\n');
     tm = [tm1 tm2 tm3];
     title(tm,'fontweight','normal')
    hold off
%    
     pause(0.02)
      
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
   global c
   y = -c(1)*v + c(2)*x + c(3)* x^3 + c(4)*cos(c(5)*t);

end
