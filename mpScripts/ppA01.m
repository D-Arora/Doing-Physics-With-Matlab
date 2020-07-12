% ppA01.m


clear all
close all
clc

hold on

N = 2000;
x = zeros(N,1);
y = zeros(N,1);
L = 20;
x(1) = -10;
y(1) = 9;

tMax = 5;
t = linspace(0,tMax,N);
dt = t(2)-t(1);

A = [1 4;2 -1];

x(2) = x(1) + dt*(A(1,1)*x(1) + A(1,2)*y(1)); 
y(2) = y(1) + dt*(A(2,1)*x(1) + A(2,2)* y(1));

flagS = 0; c = 3;
while flagS == 0
   x(c) = x(c-2) + 2*dt*(A(1,1)*x(c-1) + A(1,2)*y(c-1));
   y(c) = y(c-2) + 2*dt*(A(2,1)*x(c-1) + A(2,2)*y(c-1));
   if abs(x(c)) > L; flagS = 1; end
   if abs(x(c)) > L; flagS = 1; end
   if c > N-10; flagS = 1; end
   c = c+1;
end

XY = 1:c-1;

figure(1)
pos = [0.02 0.05 0.29 0.29];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
plot(t(XY),x(XY),'b')
hold on
plot(t(XY),y(XY),'r')

figure(2)
pos = [0.42 0.05 0.29 0.29];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
plot(x(XY),y(XY))
xlim([-L L])
ylim([-L L])
grid on



