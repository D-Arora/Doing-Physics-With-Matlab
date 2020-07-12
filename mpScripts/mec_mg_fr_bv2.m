% mec_mg_fr_v2.m

% 21 apr 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

clear all
close all
clc

m = 1e-2;
b = 1e-4;                 % alpha
v0 = -12;
g = 9.8;

tMin = 0;
tMax = 3;
N = 500;

vT = sqrt(m*g/b);
t = linspace(tMin,tMax,N);
dt = t(2)-t(1);
x = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);
x(1) = 0;
v(1) = v0;
a(1) = g - (b/m) * v(1)^3/abs(v(1));
v(2) = v(1) + a(1) * dt;
a(2) = g - (b/m) * v(2)^3/abs(v(2));
v(2) = v(1) + 0.5 *(a(1)+a(2)) * dt;

x(2) = 0.5*(v(1)+v(2))*dt + 0.5*0.5*(a(1)+a(2))*dt^2;



for k = 1 : N-2
    v(k+2) = v(k) + 2*dt * (g - (b/m) * v(k+1)^3/abs(v(k+1)));
end

 for k = 1 : N-2
     x(k+2) = x(k) + 2*dt * v(k+1);
 end

a = g - (b/m) .* v.^3 ./ abs(v);

% ANALYTICAL ==============================================================

k = 1; flag = -1;
while flag < 0
  vA(k) = vT * tan(atan(v0/vT) + (g/vT)* t(k));  
  aA(k) = g + (b/m) * vA(k)^2;
  xA(k) = vT^2/(2*g) * log((vT^2+vA(k)^2)/(vT^2+v0^2)); 
  if vA(k) >= 0; flag = 0; end;
  k = k+1;
 
end
  ks = k; x0 = xA(k-1);
  
for k = ks : N
  v0 = 0;  
  K1 = exp( -(2*g/vT) * (t(k)- t(ks)));  
  vA(k) = (v0+vT) + (v0-vT) .* K1;
  vA(k) = vA(k) / ((v0+vT) - (v0-vT)* K1);
  vA(k) = vA(k) * vT;
  xA(k) = x0 + (vT^2/(2*g)) * log((vT^2 - v0^2)/(vT^2 - vA(k)^2));
  aA(k) = g - (b/m) * vA(k)^2;
end
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