% cemDiff01.m
% 19 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Numerical differentation
%    1st derivative dy/dx   forward   backward   central  difference methods
%    2nd derivative d2/dx2

close all
clear all
clc
tic

% SETUP ==================================================================
% number of data point sfor calculations
   N = 101;
% wavelength wL / wave number k / range x / function y
   wL = 20; 
   k = 2*pi/wL;

   x = linspace(0,100,N);
   dx = x(2)-x(1);
   y = sin(k*x);

% Gradient of function ===================================================
dydxA = k .* cos(k*x);         % Analytical (exact) 

dydxM = gradient(y,dx);        % Matlab gradient function

dydxF = zeros(1,N);            % Forward approximation
dydxF(N) = (y(N)-y(N-1))/(dx);
for n = 1: N-1
   dydxF(n) = (y(n+1)-y(n))/dx;
end

dydxB = zeros(1,N);            % Backward approximation
dydxB(1) = (y(2)-y(1))/(dx);
for n = 2: N
   dydxB(n) = (y(n)-y(n-1))/dx;
end

dydxC = zeros(1,N);            % Central difference approximation
dydxC(1) = (y(2)-y(1))/dx;
dydxC(N) = (y(N)-y(N-1))/dx;
for n = 2: N-1
   dydxC(n) = (y(n+1)-y(n-1))/(2*dx);
end

% Accuracy - ratio to exact value 
N
p = min(find(dydxA/k < 0.5));
pA = (dydxA(5)/k)/(dydxA(5)/k)
pM = (dydxM(5)/k)/(dydxA(5)/k)
pF = (dydxF(5)/k)/(dydxA(5)/k)
pB = (dydxB(5)/k)/(dydxA(5)/k)
pC = (dydxC(5)/k)/(dydxA(5)/k)

% SECOND DERIVATIVE ======================================================
d2ydx2A = -k^2 .* sin(k*x);
d2ydx2M = gradient(dydxM,dx);

d2ydx2C = zeros(1,N);            % Central difference approximation
d2ydx2C(1) = (dydxC(2)-dydxC(1))/dx;
d2ydx2C(N) = (dydxC(N)-dydxC(N-1))/dx;
for n = 2: N-1
   d2ydx2C(n) = (y(n+1)-2*y(n) + y(n-1))/(dx^2);
end

% del2 function: second derivative of y 

del2y = 4*del2(y,dx);

% GRAPHICS ===============================================================

figure(1)
% xP = x; yP = y;
% plot(xP, yP);
hold on

xP = x; yP = dydxA/k;
plot(xP, yP,'b','linewidth',2)

yP = dydxM/k;
plot(xP,yP,'r');
 
yP = dydxF/k;
plot(xP,yP,'k');

yP = dydxB/k;
plot(xP,yP,'c');

yP = dydxC/k;
plot(xP,yP,'m');

xlabel('x'); ylabel('(dy / dx)  / k');
set(gca,'yLim',[-1.1 1.1]);
legend('A', 'M', 'F', 'B', 'C');
box on


figure(2) % --------------------------------------------------------------

yP = d2ydx2A/k^2;
plot(xP,yP,'b','linewidth',2)

hold on
yP = d2ydx2M/k^2;
plot(xP,yP,'r') 

hold on
yP = d2ydx2C/k^2;
plot(xP,yP,'k') 


hold on
yP = del2y/k^2;
plot(xP,yP,'om')
xlabel('x'); ylabel('(d^2y / dx^2)  / k^2');

set(gca,'yLim',[-1.1 1.1]);
legend('A', 'Mgrad','C','Mdel2');
box on

toc