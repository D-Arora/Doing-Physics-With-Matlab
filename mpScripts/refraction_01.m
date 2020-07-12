% refraction_01.m
clear all
close all
clc

wL1 = 8; wL2 = 4;
k1 = 2*pi/wL1; k2 = 2*pi/wL2;


x1 = 30; y1 = 0; x2 = 80; y2 = 100;

t0 = atand((y2-y1)/(x2-x1));
t1 = 90 - t0;
t2 = asind((wL2/wL1)*sind(t1));

m1 = tand(t0);
b1 = y1 - m1 * x1;

xp = 80; yp = 50;

p2 = t1-t2;
m2 = -tand(p2);
b2 = yp - m2 * xp;

xc = (b2-b1)/(m1 - m2);
yc = m1 * xc + b1;

dx = xp - xc;
d2 = dx / cosd(p2);
d1 = xc;


% normal
mn = -1/m1;
bn = yc - mn * xc;
yn1 = bn;
yn2 = mn * 100 + bn;


figure(1)
set(gca,'XLim',[0 100]);
set(gca,'YLim',[0 100]);
plot([x1 x2], [y1 y2]);
hold on
plot([0 xc], [yc yc],'r');
plot([xc xp], [yc yp],'r');
plot([xc 100], [yc yc],'k');
plot([0 100], [yn1 yn2],'k');

axis square

