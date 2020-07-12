% m2_hyperbola.m


close all
clear all
clc


x1 = linspace(4,10,200);
x2 = -x1;

y1 = 3 .* sqrt((x1/4).^2 - 1);

y2 = -y1;

figure(1)
plot(x1,y1);
hold on
plot(x2,y1);
plot(x1,y2);
plot(x2,y2);



