% mec_fr_bv.m
% 20 apr 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
clear all
close all
clc

% Input Parameters -------------------------------------------------------
   m = 2;
   b = 10;
   v0 = -10;
   dt = 1e-4;
   tmax = 1.95;
   
% SETUP ------------------------------------------------------------------
   g = 9.8;
   k = (b/m);
   vT = m*g / b;
   t(1) = 0; t(2) = dt;
   v(1) = v0;
   a(1) = g -(b/m) * v(1);
   v(2) = v(1) + a(1) * dt;
   a(2) = g -(b/m) * v(2);
   v(2) = v(1) + 0.5*(a(1)+a(2)) * dt;
    
   flag = 1;
   c = 3;

% NUMERICAL APPROACH -----------------------------------------------------

while flag > 0,
    v(c) = v(c-2) + (2*dt*g) * (1-v(c-1)/vT);
    t(c) = t(c-1) + dt;
    if t(c) > tmax
        flag = -10;
    end;
    c = c + 1;
end

x(1) = 0;
x(2) = 0.5*(v(1)+v(2)) * dt;

for cc = 3 : c-1
 x(cc) = x(cc-2) + 2*dt*v(cc-1);    
end

a = g - (b/m) .* v;

% ANALYTICAL ANALYSIS ---------------------------------------------------
vA = vT + (v0 - vT) .* exp(-k*t);
aA = g - k .* vA;
xA = vT .* t + (1/k) .* (v0-vT) .* (1-exp(-k*t));


figure(1)   % acceleration ------------------------------------------------
fs = 9;
set(gcf,'units','normalized','position',[0.1 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'acceleration a  [m.s^{-2}]';
xP = t; yP = a;
Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
   plot(xP,yP,'b','lineWidth',2)
hold on
yP = aA;
    plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
% set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
legend('N','A');

figure(2)   % velocity ----------------------------------------------------
fs = 9;
set(gcf,'units','normalized','position',[0.3 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'velocity  v  [m.s^{-1}]';
xP = t; yP = v;
Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
   plot(xP,yP,'b','lineWidth',2)
hold on
yP = vA;
    plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
% set(gca,'Xlim',Xrange);
% set(gca,'Ylim',Yrange);
legend('N','A');


figure(3)    % displacement ---------------------------------------------
fs = 9;
set(gcf,'units','normalized','position',[0.5 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'position  x  [m]';
xP = t; yP = x;
Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
   plot(xP,yP,'b','lineWidth',2)
hold on
yP = xA;
    plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
% set(gca,'Xlim',Xrange);
% set(gca,'Ylim',Yrange);
legend('N','A');

figure(4)    % displacement ---------------------------------------------
fs = 9;
set(gcf,'units','normalized','position',[0.7 0.2 0.18 0.22]);
tx = 'velocity  v  [m.s^{-1}]';
ty = 'position  x  [m]';
xP = v; yP = x;
Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
   plot(xP,yP,'b','lineWidth',2)
hold on
xP = vA; yP = xA;
    plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
% set(gca,'Xlim',Xrange);
% set(gca,'Ylim',Yrange);
legend('N','A');