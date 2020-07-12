% osc_bungee.m

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 20140718


clear all
close all
clc
tic


% Input and Constants +++++++++++++++++++++++++++++++++++++++++++++++
xs = -9;                % initial position of jumper (m)  default = 0
ys = 0; 

L =   9;                % natural length of bungee cord (m)  default = 9
m =  80;                % mass of jumper (kg)   default = 80

D_1 = 4;                 % friction constant  default = 0
D_2 = 5;                 % drag constant    default = 0

g = 9.8;                % acceleration due to gravity (m/s^2)


% two-piece linear model for spring stiffness
k_1 = 160;               % default = 100
k_2 = 80;                % default = 100
e_1 = 5;                 % default = 5              

t_max = 25;              % max simulation time (s)   default = 25
nt = 100000;             % number of time steps
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ********************************************************************
%  Setup and initialize variables
% ********************************************************************

t = linspace(0,t_max,nt);           % time
dt = t(2)-t(1);                    % time step

% position
x = zeros(nt,1); y = zeros(nt,1);
r = zeros(nt,1); e = zeros(nt,1);

% velocity
vx = zeros(nt,1); vy = zeros(nt,1); v = zeros(nt,1);

% acceleration
ax = zeros(nt,1); ay = zeros(nt,1); a = zeros(nt,1);

% forces
F_Gx = zeros(nt,1); F_Gy = (-m*g) .* ones(nt,1);
Fx = zeros(nt,1); Fy = zeros(nt,1); F = zeros(nt,1);   % net force on jumper
F_E = zeros(nt,1);      % elastic restoring force
F_Ex = zeros(nt,1); F_Ey = zeros(nt,1);
F_D1x = zeros(nt,1); F_D1y = zeros(nt,1);
F_D2x = zeros(nt,1); F_D2y = zeros(nt,1);
FG = m*g;
% Energies
K = zeros(nt,1); U_E = zeros(nt,1); U_G = zeros(nt,1);
E = zeros(nt,1); 


% ********************************************************************
%  Calculations
% ********************************************************************

% time step 1 (t = 0)

if abs(xs) > L; xs = -L; end;     % bungee not extended at start

x(1) = xs; y(1) = ys;
r(1) = sqrt(x(1)^2 + y(1)^2);

Fx(1) =   0;
Fy(1) = -m*g;
F(1) = sqrt(Fx(1)^2 + Fy(1)^2);

ax(1) = Fx(1)/m; ay(1) = Fy(1)/m; a(1) = F(1)/m;

% time step 2 to nt
for c = 2 : nt;
    vx(c) = vx(c-1) + ax(c-1) * dt;
    vy(c) = vy(c-1) + ay(c-1) * dt;
    v(c)  = sqrt(vx(c)^2 + vy(c)^2);
 
    x(c) = x(c-1) + vx(c-1)*dt + 0.5 * ax(c-1) * dt^2;
    y(c) = y(c-1) + vy(c-1)*dt + 0.5 * ay(c-1) * dt^2;
    r(c) = sqrt(x(c)^2 + y(c)^2);

    e(c) = r(c) - L;
    
    if e(c) < 0, e(c) = 0; end;
    
   % % Gravitational force
   % F_Gx = 0; F_Gy = -m*g;       
    
    % Elastic force and potnetial energy
    if e(c) < 0                  
        F_Ex = 0; F_Ey = 0;
    else
        if e(c) <= e_1
           F_Ex(c) = -k_1*e(c)*x(c)/r(c);
           F_Ey(c) = -k_1*e(c)*y(c)/r(c);
           F_E(c) = sqrt(F_Ex(c).^2+F_Ey(c).^2);
           U_E(c) = 0.5*k_1*e(c)^2;
        else
           F_Ex(c) = -(k_2*e(c) + (k_1 - k_2)*e_1) * x(c) / r(c);
           F_Ey(c) = -(k_2*e(c) + (k_1 - k_2)*e_1) * y(c) / r(c);
           F_E(c) = sqrt(F_Ex(c).^2+F_Ey(c).^2);
           U_E(c) = 0.5*k_2*e(c)^2 + (k_1 - k_2)*e_1*e(c) - 0.5*(k_1 - k_2)* e_1^2;
        end
    end
    
    % Dissipative force
    F_D1x(c) = - D_1*vx(c);      F_D1y(c) = - D_1*vy(c);
    F_D2x(c) = - D_2*vx(c)*v(c); F_D2y(c) = - D_2*vy(c)*v(c);
    
    % Net force and acceleration
    Fx(c) =  F_Gx(c) + F_Ex(c) + F_D1x(c) + F_D2x(c) ;
    Fy(c) =  F_Gy(c) + F_Ey(c) + F_D1y(c) + F_D2y(c) ;
    F(c) = sqrt(Fx(c)^2 + Fy(c)^2);

    ax(c) = Fx(c)/m; ay(c) = Fy(c)/m; a(c) = F(c)/m;   
  
end
    K  = (0.5*m) .* v.^2;          % KE
    U_G = (m*g) .* y;              % PE: gravitation  
    E = K + U_G + U_E;             % total energy
    
 
%  Analytical & Numerical  calculations
v_A_ff = sqrt(2*g*L);           % free fall velocities
v_N_ff = v(min(find(e>0)));

% solution of quadratic equation to find e_max
a_q = 1;
b_q = (2/k_2)*((k_1 - k_2)*e_1 - m*g);
c_q = -((k_1-k_2)*e_1^2 + 2*m*g*L)/k_2;

[e_max1, e_max2] = eq_quadratic(a_q,b_q,c_q);

e_A_max = e_max1;               % max extension
h_A_drop = L+e_A_max;           % jump drop
h_N_drop = -min(y);


v_max = max(v);                 % magnitude max velocity  (m/s)
v_max1 = v_max*3.6;             % magnitude max velcoity  (km/h)
a_g = max(a)/g;                 % max acceleration a/g
F_max = max(F);                 % max net force on jumper
F_E_max = max(F_E);             % max elastic force - bungee cord

gforce_y = (ay+g) ./g;          % vertical g-force 


% ********************************************************************
%  Graphics
% ********************************************************************

figure(1) % 111111111111111111111111111111111111111111111111111111111
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.02 0.4 0.2 0.45]);
fs = 13;
% Plot: x, y trajectory
x_p = x;     % x data for plot
y_p = y;     % y data for plot
ns = nt /100;
title_x = 'x position  (m)';
title_y = 'y position  (m)';

lineWidth_p = 1.5;
set(gcf,'Name','Trajectory');
plot(x_p,y_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs);
ylabel(title_y,'fontsize',fs);

% Plot: end point for trajectory
color_p = [1 0 0] ;
size_p = 8;

hold on
h_plot_A2 = plot(x_p(end),y_p(end),'ro');
set(h_plot_A2,'MarkerEdgeColor',color_p,'MarkerFaceColor',color_p, ...
    'MarkerSize',size_p);

color_p = [0 0 1];
color_p1 = [1 1 1];
size_p = 4;

h_plot_A2 = plot(x_p(1:ns:end),y_p(1:ns:end),'ro');
set(h_plot_A2,'MarkerEdgeColor',color_p,'MarkerFaceColor',color_p1, ...
    'MarkerSize',size_p);
axis equal
set(gca,'Xlim',[-10 10]);

hold on
plot([0,x(end)],[0,y(end)],'k','linewidth',1);

set(gca,'fontsize',fs);

figure(2) % 222222222222222222222222222222222222222222222222222222222
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.1 0.1 0.5 0.6]);
set(gcf,'Name','Position / Velocity / Acceleration');
fs = 13;

subplot(3,2,1);
% Plot: x position and time
x_p = t;     % x data for plot
y_p = x;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'position  x  (m)';
lineWidth_p = 2;
color_p = [1 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(3,2,2);
% Plot: y position and time
x_p = t;     % x data for plot
y_p = y;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'position  y  (m)';
lineWidth_p = 2;
color_p = [1 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
hold on
lineWidth_p = 1;
color_p = [0 0 0];
y_p = -L .* ones(nt,1);
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p)
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on
set(gca,'Ylim',[1.1*min(y) 5]);

subplot(3,2,3);
% Plot: x velocity and time
x_p = t;     % x data for plot
y_p = vx;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'velocity  v_x  (m/s)';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(3,2,4);
% Plot: y velocity and time
x_p = t;      % x data for plot
y_p = vy;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'velocity  v_y  (m/s)';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(3,2,5);
% Plot: x acceleration and time
x_p = t;      % x data for plot
y_p = ax;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'acc  a_x  (m/s^2)';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(3,2,6);
% Plot: y acceleration and time
x_p = t;      % x data for plot
y_p = ay;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'acc  a_y  (m/s^2)';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);% 
grid on


figure(3)   % 3333333333333333333333333333333333333333333333333333333
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.3 0.05 0.5 0.8]);
set(gcf,'Name','Forces');
fs = 13;

subplot(5,2,1);
% Plot: x position and time
x_p = t;     % x data for plot
y_p = x;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'position  x  (m)';
lineWidth_p = 2;
color_p = [1 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on
% tm1 = 'Max force acting on jumper  F_{max} / F_G  =    ';
% tm2 = num2str(max(F)/FG,3);
% tm = [tm1 tm2];
% title(tm);

subplot(5,2,2);
% Plot: y position and time
x_p = t;     % x data for plot
y_p = y;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'position  y  (m)';
lineWidth_p = 2;
color_p = [1 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
hold on
lineWidth_p = 1;
color_p = [0 0 0];
y_p = -L .* ones(nt,1);
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p)
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on
set(gca,'Ylim',[1.1*min(y) 5]);

subplot(5,2,3);
% Plot: Gravitational force FGx
x_p = t;     % x data for plot
y_p = F_Gx;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'FG_x (N)';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(5,2,4);
% Plot: Gravitational force
x_p = t;        % x data for plot
y_p = F_Gy;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'FG_y  (N)';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
set(gca,'YLim',[-1.1*max(abs(F_Gy)), 1.1*max(abs(F_Gy))]);
grid on

subplot(5,2,5);
% Plot: Elastic restoring force x
x_p = t;     % x data for plot
y_p = F_Ex;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'FE_x ';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
grid on

subplot(5,2,6);
% Plot: Elastic restoring force y
x_p = t;      % x data for plot
y_p = F_Ey;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'FE_y';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(5,2,7);
% Plot: Dissipative  x  forces
x_p = t;         % x data for plot
y_p = F_D1x;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'FD_x';
lineWidth_p = 2;
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on

subplot(5,2,8);
% Plot: Dissipative  y  forces
x_p = t;         % x data for plot
y_p = F_D1y;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'FD_y';
lineWidth_p = 2;
color_p = [0 0 0];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);

grid on
hold on
y_p = F_D2y./FG;     % y data for plot
color_p = [1 0 0];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
set(gca,'fontsize',fs);
hh = legend('FD_1 \propto v','FD_2 \propto v^2','Location','North','Orientation','horizontal');
set(hh,'fontsize',10);

subplot(5,2,10);
 x_p = t;     % x data for plot
  y_p = gforce_y;     % y data for plot
  title_x = 'time  t  (s)';
  title_y = 'y: g-force';
  lineWidth_p = 2;
  color_p = [0 0 1];
  plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
  xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
  set(gca,'fontsize',fs);
  grid on
tm1 = 'g-force-{max}  =    ';
tm2 = num2str(max(gforce_y),3);
tm = [tm1 tm2];
title(tm,'fontweight','normal');


figure(4)   % 4444444444444444444444444444444444444444444444444444444
% Plot: energies K Ue Up E
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.6 0.05 0.3 0.7]);
set(gcf,'Name','Energies');
fs = 11;

subplot(2,1,1);
% Plot: x position and time
x_p = t;     % x data for plot
y_p = y;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'position  y  (m)';
lineWidth_p = 2;
color_p = [1 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
set(gca,'fontsize',fs);
grid on
hold on
lineWidth_p = 1;
color_p = [0 0 0];
y_p = -L .* ones(nt,1);
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p)
set(gca,'Ylim',[1.1*min(y) 5]);

subplot(2,1,2);
x_p = t;        % x data for plot
y_p = K;     % y data for plot
title_x = 'time  t  (s)';
title_y = 'energy   (J) ';
lineWidth_p = 1.5;
color_p = [1 0 0];
set(gcf,'Name','Energy');
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);

hold on

x_p = t;        % x data for plot
y_p = U_G;     % y data for plot
color_p = [0 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);

x_p = t;        % x data for plot
y_p = U_E;     % y data for plot
color_p = [1 0 1];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);

x_p = t;        % x data for plot
y_p = E;     % y data for plot
color_p = [0 0 0];
plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);

legend('E_K','E_G','E_E','E','Location','NorthOutside','Orientation','horizontal');
set(gca,'fontsize',fs);
grid on

% Plot: restoring force and extension
figure(5)   % 555555555555555555555555555555555555555555555555555555555555555555
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.02 0.05 0.2 0.25]);
set(gcf,'Name','Forces');
fs = 11;
  x_p = e;     % x data for plot
  y_p = F_E;     % y data for plot
  title_x = 'extension e  (m)';
  title_y = 'F_E  (N)';
  lineWidth_p = 2;
  color_p = [0 0 0];
  plot(x_p,y_p,'Color',color_p,'lineWidth',lineWidth_p);
  xlabel(title_x,'fontsize',fs); ylabel(title_y,'fontsize',fs);
  set(gca,'fontsize',fs);
  grid on

% ********************************************************************
%  Output results to screen
% ********************************************************************

disp('  ');
disp('    ');
disp('Input parameters');
fprintf('   Initial x coordinate of jumper, xs =  %0.2f  m  \n',xs); 
fprintf('   Initial y coordinate of jumper, ys =  %0.2f  m  \n',ys); 
fprintf('   Mass of jumper, m =  %0.2f  kg  \n',m); 
fprintf('   Weight of jumper, F_G = %0.0f  N    \n',FG); 
fprintf('   Natural length of bungee cord, L =  %0.2f  m  \n',L); 
fprintf('   Dissipative force constant, D_1 =  %0.2f  N.m^-1.s^-1  \n',D_1); 
fprintf('   Dissipative force constant, D_2 =  %0.2f  N.m^-2.s^-2  \n',D_2); 
fprintf('   Stiffness, k_1 =  %0.2f  N.m^1  \n',k_1); 
fprintf('   Stiffness, k_2 =  %0.2f  N.m^1  \n',k_2);
fprintf('   Stiffness: two piece linear transition, e_1 =  %0.2f  m  \n',e_1); 
disp('    ');
disp('Analytical (FD = 0) & numerical predictions    ');
fprintf('   A: max free fall velocity, v_ff = %0.2f  m/s  \n',v_A_ff); 
fprintf('   N: max free fall velocity, v_ff = %0.2f  m/s  \n',v_N_ff); 
disp('    ');
fprintf('   A: drop height, h_drop = %0.2f  m  \n',h_A_drop); 
fprintf('   N: drop height, h_drop = %0.2f  m  \n',h_N_drop); 
disp('    ');

disp('Output parameters    ');
fprintf('   max elongation factor  r_max / L = %0.2f    \n',max(r)/L); 
fprintf('   max velocity, v_max = %0.2f  m.s^-1  \n',v_max); 
fprintf('   max velocity, v_max = %0.2f  km.h^-1  \n',v_max1); 
fprintf('   max y-acceleration, ay_max / g = %0.2f m.s^-2   \n',max(ay));
fprintf('   max acceleration, a_max = %0.2f  m.s-2    \n',max(a));
fprintf('   max force on jumper, F_max = %0.0f  N    \n',F_max);
fprintf('   F_max / F_G  =  %0.2f      \n',F_max/FG);
fprintf('   max vertical g-force =  %0.2f      \n',max(gforce_y));
fprintf('   max elastic restoring force by bungee cord, F_E_max = %0.0f  N    \n',F_E_max);


