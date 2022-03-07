% ODE_01.m


clear
close all
clc


tMin = 0;
tMax = 10;
tSpan = [tMin tMax];

y0 = [1,0];

m = 2;
b = 0.2;
k = 3;

K = [m b k];

[t,sol] = ode45(@(t,y) FNT(t,y,K),tSpan,[y0]);

y = sol(:,1);
%dTdt = sol(:,2);

figure(1)

xP = t; yP = y;
  plot(xP,yP)

grid on; box on
xlabel('t')
ylabel('y')
set(gca,'FontSize',12)


function dS = FNT(t,s,K)
  m = K(1); b = K(2); k = K(3);
  y = s(1);
  yDot = s(2);
  
  dS(1) = yDot;
  dS(2) = -(b/m)*yDot - (k/m)*y;

 dS = dS';

end

