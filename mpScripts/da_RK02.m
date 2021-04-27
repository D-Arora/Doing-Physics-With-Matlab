% da_RK01.m


% RUNGE_KUTTA METHOD
% Solving Initial Value Ordinary Differential Equations 

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/RK.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
%  200522     Matlab 2019b

clear 
close all
clc


% >>>>>  Number of grid points 
  N = 200;
  
% >>>>>  t variable
  tMin = 0;
  tMax = 2;
  
% >>>>>  initial value for y variable
  y = zeros(1,N); 
  y(1) = 1;
  

% SETUP  ==============================================================  
% t interval and step size
  t = linspace(tMin,tMax,N);
  h = t(2) - t(1);
% yDot variable as a function of t at y(nY) (calls function YDOT)   
  nY = 1;
  yDot = YDOT(t,y(nY));
  
% Compute y and yDdot values as a function of time =====================
for n = 1 : N-1
   k1 = YDOT(t(n),y(n));
   k2 = YDOT(t(n) + 0.5*h, y(n) + k1*h/2);
   k3 = YDOT(t(n) + 0.5*h, y(n) + k2*h/2);
   k4 = YDOT(t(n) + 1.0*h, y(n) + k3*h);

   y(n+1) = y(n) + h*(k1+2*k2+2*k3+k4)/6;
end
  
% yDot variable (calls function YDOT)  
  yDot = YDOT(t,y(1));

% EXACT y VALUE =======================================================
  tE = linspace(tMin,tMax,20);
%  yE = 2.*exp(tE) - tE -1;  
   yE = -2 -2.*tE - tE.^2 + 3.*exp(tE);
   

% [2D] plots  =========================================================   
   Y = linspace(1,3,N);
  [tt, yy] = meshgrid(t,Y);
  YYDot = YDOT(tt,yy);
  
  
% GRAPHICS  ===========================================================
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.10 0.10 0.25 0.5]);
   set(gcf,'color','w');

subplot(2,1,1)
   plot(t,yDot,'b','linewidth',2)
   grid on
   xlabel('t')
   ylabel('dy / dt')
   txt = sprintf('y = %3.1f  ',y(nY));
   title(txt,'fontweight','normal')
   set(gca,'fontsize',12)
   
subplot(2,1,2)   
   xP = tE; yP = yE;
   plot(xP,yP,'or');
   hold on
   xP = t; yP = y;
   plot(xP,yP,'b','linewidth',2)
   grid on
   xlabel('t')
   ylabel('y')
   legend('exact','RK','location','north')
   set(gca,'fontsize',12)
   
figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.380 0.10 0.25 0.5]);
   set(gcf,'color','w');
   
subplot(2,1,1)
   pcolor(tt,yy,YYDot)
   shading interp
   colorbar
   title('dy / dt','fontweight','normal')
   xlabel('t')
   ylabel('y')
   set(gca,'fontsize',12)
 
 subplot(2,1,2)  
   surf(tt,yy,YYDot)
   shading interp
   colorbar
   title('dy / dt','fontweight','normal')
   xlabel('t')
   ylabel('y')
   set(gca,'fontsize',12)
   
   
% FUNCTION  ===========================================================

function f = YDOT(t,y)
  %  f = t + y;
     f = t.^2 + y;
end





% clear all 
% close all
% clc
% h = 0.01;  % set the step size
% x = [0:0.1:10];  % set the interval of x
% y = zeros(1,length(x));
% y(1) = 0;   % set the intial value for y
% n = length(x)-1;
% y_dot =@(x,y)(x^3/3); %insert function to be solved
% for i = 1:n
%     k1 = y_dot(x(i),y(i));
%     k2 = y_dot(x(i)+.5*h,y(i)+.5*k1*h);
%     k3 = y_dot(x(i)+.5*h,y(i)+.5*k2*h);
%     k4 = y_dot(x(i)+h,y(i)+k3*h);
%     y(i+1) = y(i)+((k1+2*k2+2*k3+k4)/6)*h;
%     
% end
% %[t,y_check] = ode45(y_dot,x,0);
% plot(x,y)
%title('Eulers Method')
%hold on
%plot(x,y_check,'r')
%title('ode45 Check')


% %%
% f = @(t,y) (sin(t));% (2 - exp(-4*t) - 2*y);
% h = 0.1;  % Define Step Size
% t_final = 20;
% t = 0:h:t_final;
% y = zeros(1,numel(t));
% y(1) = 1;  % y0
% % You know the value a t = 0, thats why you'll state with t = h i.e. i = 2
% for i = 2:numel(t)
%     k1 = h*f(t(i-1),y(i-1));
%     k2 = h*f(t(i-1)+h/2, y(i-1)+k1/2);
%     k3 = h*f(t(i-1)+h/2, y(i-1)+k2/2);
%     k4 = h*f(t(i-1)+h, y(i-1)+k3);
%     y(i) = y(i-1) + (k1+2*k2+2*k3+k4)/6;
%    % disp([t(i) y(i)]);
% end
% 
% 
% %%
% % It calculates ODE using Runge-Kutta 4th order method
% % Author Ido Schwartz
% % Originally available form: http://www.mathworks.com/matlabcentral/fileexchange/29851-runge-kutta-4th-order-ode/content/Runge_Kutta_4.m
% % Edited by Amin A. Mohammed, for 2 ODEs(April 2016)
% 
% clc                                               % Clears the screen
% clear
% clc
% 
% h=0.1;                                             % step size
% x = 0:h:1;                                         % Calculates upto y(1)
% y = zeros(1,length(x)); 
% z = zeros(1,length(x)); 
% y(1) = 3;                                          % initial condition
% z(1) = 1;                                          % initial condition
% % F_xy = @(t,r) 3.*exp(-t)-0.4*r;                  % change the function as you desire
% F_xyz = @(x,y,z) z;                                  % change the function as you desire
% G_xyz = @(x,y,z) 6*y-z;
% 
% for i=1:(length(x)-1)                              % calculation loop
%     k_1 = F_xyz(x(i),y(i),z(i));
%     L_1 = G_xyz(x(i),y(i),z(i));
%     k_2 = F_xyz(x(i)+0.5*h,y(i)+0.5*h*k_1,z(i)+0.5*h*L_1);
%     L_2 = G_xyz(x(i)+0.5*h,y(i)+0.5*h*k_1,z(i)+0.5*h*L_1);
%     k_3 = F_xyz((x(i)+0.5*h),(y(i)+0.5*h*k_2),(z(i)+0.5*h*L_2));
%     L_3 = G_xyz((x(i)+0.5*h),(y(i)+0.5*h*k_2),(z(i)+0.5*h*L_2));
%     k_4 = F_xyz((x(i)+h),(y(i)+k_3*h),(z(i)+L_3*h)); % Corrected        
%     L_4 = G_xyz((x(i)+h),(y(i)+k_3*h),(z(i)+L_3*h));
% 
%     y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
%     z(i+1) = z(i) + (1/6)*(L_1+2*L_2+2*L_3+L_4)*h;  % main equation
% 
% end
% 
% plot(y,z)
