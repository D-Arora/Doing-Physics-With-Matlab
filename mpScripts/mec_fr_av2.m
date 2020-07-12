% mec_fr_bv2.m

clear all
close all
clc


m = 2;
b = 5;
v0 = 10;
tMax = 2;

dt = 1e-7;
t(1) = 0; t(2) = dt;
v(1) = v0; v(2) = v0;

flag = 1;
k = -(2*dt*b/m);
c = 3;
while flag > 0,
    v(c) = k * v(c-1)^2 + v(c-2);
    t(c) = t(c-1) + dt;
    if t(c) > tMax
        flag = -10;
    end;
    c = c + 1;
end

x(1) = -v0 * 2* dt;
x(2) = -v0 * dt;

for cc = 3 : c-1
 x(cc) = x(cc-2) + 2*dt*v(cc-1);    
end
    

vA = v0 ./ (1 + (b*v0/m) .* t);
xA = (m/b) .* log(1+(b*v0/m) .* t);


figure(1)
fs = 9;
set(gcf,'units','normalized','position',[0.2 0.2 0.18 0.22]);
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
set(gca,'Xlim',Xrange);
set(gca,'Ylim',Yrange);


figure(2)
fs = 9;
set(gcf,'units','normalized','position',[0.4 0.2 0.18 0.22]);
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
set(gca,'Xlim',Xrange);
set(gca,'Ylim',Yrange);

figure(3)
fs = 9;
set(gcf,'units','normalized','position',[0.6 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'acceleration  a  [m.s_{-2}]';
xP = t; yP = b .* v.^2;
Xrange = [0 1.0 * max(xP)];
Yrange = [0 1.1 * max(yP)];
   plot(xP,yP,'b','lineWidth',2)
%hold on
%yP = xA;
 %   plot(xP,yP,'r')
grid on
xlabel(tx); ylabel(ty);
set(gca,'Xlim',Xrange);
set(gca,'Ylim',Yrange);

