% mec_mg_fr_v2.m

% 21 apr 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/

clear all
close all
clc

m = 1e-2;
b = 1e-4;                 % alpha
v0 = -12;
g = 9.8;

tMin = 0;
tMax = 10;
N = 500;

vT = sqrt(m*g/b);
t = linspace(tMin,tMax,N);
dt = t(2)-t(1);
x = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);
x(1) = 0;
v(1) = v0;
a(1) = g - (b/m) * v(1)^2;
v(2) = v(1) + a(1) * dt;
a(2) = g - (b/m) * v(2)^2;
v(2) = v(1) + 0.5 *(a(1)+a(2)) * dt;

x(2) = 0.5*(v(1)+v(2))*dt + 0.5*0.5*(a(1)+a(2))*dt^2;



for k = 1 : N-2
    v(k+2) = v(k) + 2*dt * (g - (b/m) * v(k+1)^2);
end

 for k = 1 : N-2
     x(k+2) = x(k) + 2*dt * v(k+1);
 end

a = g - (b/m) .* v.^2;


% ANALYTICAL ==============================================================

%K1 = exp( -sqrt(4*b*g) .*t );  
K1 = exp( -(2*g/vT) .*t );  
vA = (v0+vT) + (v0-vT) .* K1;
vA = vA ./ ((v0+vT) - (v0-vT) .* K1);
vA = vA .* vT;

xA = (vT^2/(2*g)) .* log((vT^2 - v0.^2)./(vT^2 - v.^2));
aA = g - (b/m) .* vA.^2;

% GRAPHICS ================================================================

figure(1)
col = 'b';
fs = 9;
set(gcf,'units','normalized','position',[0.2 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'velocity  v  [m.s^{-1}]';
xP = t; yP = v;

Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
Xrange = [0 1.0 * 2];
Yrange = [0 1.1 * 10];
   plot(xP,yP,col,'lineWidth',2)
hold on
yP = vA;
   plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
%set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);

figure(2)
fs = 9;
set(gcf,'units','normalized','position',[0.4 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'acceleration  a  [m.s^{-2}]';
xP = t; yP = a;

Xrange = [0 1.0 * max(xP)];
Yrange = [-1.1*max(yP) 1.1 * max(yP)];
Xrange = [0 1.0 * 2];
Yrange = [-1.1*10 1.1 * 10];
   plot(xP,yP,col','lineWidth',2)
hold on
yP = aA;
    plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
%set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);

figure(3)
fs = 9;
set(gcf,'units','normalized','position',[0.6 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'displacement  x  [m]';
xP = t; yP = x;

Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
Xrange = [0 1.0 * 2];
Yrange = [-20 1.1 * 40];
   plot(xP,yP,col','lineWidth',2)
hold on
yP = xA;
     plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
%set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
 
figure(4)
fs = 9;
set(gcf,'units','normalized','position',[0.8 0.2 0.18 0.22]);
tx = 'velocity  v   [m.s^{-1}]';
ty = 'displacement  x  [m]';
xP = v; yP = x;

Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
Xrange = [0 1.0 * 2];
Yrange = [0 1.1 * 10];
   plot(xP,yP,col','lineWidth',2)
hold on
 yP = xA;
     plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
%set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);