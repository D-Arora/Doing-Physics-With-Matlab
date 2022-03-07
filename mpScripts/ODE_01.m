% ODE_01.m


clear
close all
clc


tMin = 0;
tMax = 1;
tSpan = [tMin tMax];

y0 = [1];

K = 0.0847;

[t,sol] = ode45(@(t,y) FNT(t,y,K),tSpan,[y0]);

T = sol(:,1);
%dTdt = sol(:,2);

figure(1)

xP = t; yP = T;
  plot(xP,yP)

grid on; box on
xlabel('t')
ylabel('y')
set(gca,'FontSize',12)


function yDot = FNT(t,y,K)

%yDot(1) = -K*y;
% yDot(1) = K*(20 - y)
%yDot(1) = t^4 - 2*y/t;
%yDot(1) = 375 - 3*y/(800+2*t);
 yDot(1) = 1+y*t;


yDot = yDot';
end

