clear all
close all
clc

m = 60;;
g = 10;
k = 100;

e = linspace(0,20,500);

y1 = (m*g).*e;
y2 = -(0.5*k).*e.^2;
y = y1+y2;

emax = m*g/k

figure(1)
plot(e,y1)
hold on
plot(e,y2);
plot(e,y);
grid on
