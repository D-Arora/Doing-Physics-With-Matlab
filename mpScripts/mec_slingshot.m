% mec_spring_energy.m
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% In the section PLANET: comment / uncomment to select
%   stationary or moving planet

clear all
close all
clc


% Input and Constants +++++++++++++++++++++++++++++++++++++++++++++++
K = 80;                 % constant of proportionality
step = 250;
t_max = 8;              % max simulation time
nt = 5000;              % number of time steps


% ********************************************************************
%  Setup and initialize variables
% ********************************************************************
v_0x = 3;        % initial velcoity: x cpt.
v_0y = 0;        % initial velcoity: y cpt.
x0 = -10;                    % initial displacement: x cpt.
y0 = 3;                      % initial displacement: y cpt.

t = linspace(0,t_max,nt);                % time
dt = t(2) - t(1);                        % time step 

% PLANET --------------------------------------------------------------
% Select planet stationary or moving by commenting / uncommenting code
    vP =  0;    xP0 = 2;   % Planet Stationary vP = 0
    vP = -2;    xP0 = 10;  % Moving planet
% --------------------------------------------------------------------        
xP = xP0 + vP * t;
yP = zeros(1,nt);


% ********************************************************************
%  Calculations
% ********************************************************************

x(1) = x0; y(1) = y0;
xd(1) = x(1) - xP(1);

Rd(1) = sqrt(xd(1)^2 + y(1)^2);
vx(1) = v_0x; vy(1) = v_0y;

ax(1) = -(K / abs(Rd(1))^3) * xd(1);
ay(1) = -(K / abs(Rd(1))^3) * y(1);

x(2) = x(1) + vx(1) * dt + 0.5 * ax(1) * dt^2;
y(2) = y(1) + vy(1) * dt + 0.5 * ay(1) * dt^2;

xd(2) = x(2) - xP(2);
Rd(2) = sqrt(xd(2)^2 + y(2)^2);

vx(2) = vx(1) + ax(1) * dt;
vy(2) = vy(1) + ay(1) * dt;

% time steps 3, 4, 5, ....
for c = 3 : nt
    x(c) = -K*dt^2*xd(c-1)/Rd(c-1)^3 + 2*x(c-1) - x(c-2);
    y(c) = -K*dt^2*y(c-1)/Rd(c-1)^3 + 2*y(c-1) - y(c-2);
    xd(c) = x(c) - xP(c);
    Rd(c) = sqrt(xd(c)^2 + y(c)^2);
end

% velocity ----------------------------------------------------------
R = x + j .* y;
temp_p = zeros(1,nt+2); temp_m = zeros(1,nt+2); temp = zeros(1,nt+2);
temp_p(1:end-2) = R; temp_m(3:end)= R; temp = temp_p - temp_m;

v(1) = abs((R(2)-R(1))/dt); v(nt) = abs((R(end-1)-R(end-2))/dt);
v(2:nt-1)= abs(temp(3:end-2)./(2*dt));

% energy
ET = -100*ones(1,nt);
EK = v.^2;
EP = ET- EK;

% ********************************************************************
%  Graphics
% ********************************************************************

figure(1) % 111111111111111111111111111111111111111111111111111111111

% Plot: x, y trajectory
x_p = real(R);     % x data for plot
y_p = imag(R);     % y data for plot
ns = nt /50;
title_x = 'x position  (a.u.)';
title_y = 'y position  (a.u.)';

lineWidth_p = 1.5;
%set(gcf,'Name','Orbit');
%axis_dim = [-4 4 -4 4];
%color_face = [1 0.6 0.3];
%h_rectangle = rectangle('Position', [-1,-1,2,2],'Curvature',[1 1]);
%set(h_rectangle,'FaceColor',color_face,'EdgeColor',color_face);
%hold on

plot(x_p,y_p,'lineWidth',lineWidth_p);
hold on
hp = plot(x_p(1:step:end),y_p(1:step:end),'o','MarkerSize',4);
set(hp,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 0]);
xlabel(title_x);
ylabel(title_y);
axis([-10 10 -10 10])
axis equal
box on
%axis(axis_dim);
hold on
%color_p = [0 0 1];
%color_p1 = [1 1 1];
%size_p = 4;
%axis off
%h_plot_A2 = plot(x_p(1:ns:end),y_p(1:ns:end),'ro');
%set(h_plot_A2,'MarkerEdgeColor',color_p,'MarkerFaceColor',color_p1, ...
%    'MarkerSize',size_p);
hp = plot(xP(1:step:end),yP(1:step:end),'o','MarkerSize',4);
set(hp,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [1 0 0]);

set(gca,'Xtick',[-10 -5 0 5 10]);
set(gca,'Ytick',[-10 -5 0 5 10]);
set(gca,'XLim',[-10 10]);
set(gca,'YLim',[-10 10]);

figure(2); % 22222222222222222222222222222222222222222222222222222222222
subplot(2,1,1);
plot(t,v,'linewidth',3);
grid on
xlabel('time  t  (a.u.)');
ylabel('speed of spacecraft (a.u.)');

subplot(2,1,2);
plot(t,abs(R),'linewidth',3);
grid on
xlabel('time  t  (a.u.)');
ylabel('separation distance  (a.u.)');

figure(4); % 44444444444444444444444444444444444444444444444444444444444
plot(t,EK,'r','linewidth',3);
hold on
plot(t,EP,'b','linewidth',3);
plot(t,ET,'k','linewidth',3);
grid on
xlabel('time  t  (a.u.)');
ylabel('energy (a.u.)');
legend('E_K','E_P','E');

figure(3) % 3333333333333333333333333333333333333333333333333333333333333
cc = 1;
hp = plot(x_p(100),y_p(100),'o','MarkerSize',4);
set(hp,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 0]);
hold on
hp = plot(xP(100),yP(100),'o','MarkerSize',6);
set(hp,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [1 0 0]);
axis([-10 10 -10 5])
axis equal
box on
axis([-10 10 -10 5])
axis off

M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,256,'dither');  %RGB to indexed images
im(1,1,1,1) = 0;

for c = 100 : 120: nt
hp = plot(x_p(c),y_p(c),'o','MarkerSize',4);
set(hp,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 0]);
hold on
hp = plot(xP(c),yP(c),'o','MarkerSize',6);
set(hp,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [1 0 0]);
xlabel(title_x);
ylabel(title_y);
axis([-10 10 -10 5])
axis equal
box on
axis([-10 10 -10 5])
axis off
pause(0.1)
    M = getframe(gcf) ;
    cc = cc+1;
    im(:,:,1,cc) = rgb2ind(M.cdata,map,'dither');
end

% SAVE or NOT SAVE animated gif +++++++++++++++++++++++++++++++++++++++++

delay = 0.2;
ag_name = 'ag_slingshot.gif';   % file name for animated gif
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);
