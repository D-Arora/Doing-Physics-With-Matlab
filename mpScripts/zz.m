

clear
close
clc

h = 0.2;
x(1) = 0;
u(1) = 2;
%ydot = 20;


for c = 1 : 2

   k1 = FF(x(c),u(c));
   k2 = FF(x(c)+0.5*h,u(c)+k1*h/2);
   k3 = FF(x(c)+0.5*h,u(c)+k2*h/2);
   k4 = FF(x(c)+1.0*h,u(c)+k3*h);
%     k1 = h*FF(x(c));
%     k2 = h*FF(x(c)+ 0.5*h*k1);
%     k3 = h*FF(x(c)+ 0.5*h*k2);
%     k4 = h*FF(x(c)+ h*k3);
    u(c+1) = u(c) + h*(k1+2*k2+2*k3+k4)/6;
    x(c+1) = x(c) + h;
end

x
u

figure(1)
plot(x,u);
hold on
plot(x,1+x.^3/3+0.2,'r')

function f = FF(x,u)

f = x^2 ;%-2*u+x+4;

 
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
