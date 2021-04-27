% Laplace05.m


% SYSTEM: mass m  / spring k / dashpot-damping b.
% Laplace Transform (LT) used to solve the ODE describing the System
%   to give the displacement x(t) of the mass excited by the input signal y(t).

%   m*x'' + b*x' + kx = y

% 1   Define the input signal (driving force) y(t).
% 2   LT of input signal --> Y(s).
% 3   Define forces acting on system and initial conditions x(0) and x'(0).
% 4   LT of ODE -->
%       m*(s^2*X - s*x(0) - x'(0)) + b*(s*X - x(0)) + k*X = Y
%      (m*s^2 + b*s + k)*X - m*(s*x(0) - x'(0)) - b*x(0) - Y = 0
%
%  5  Solve equation for X(s)
%  6  Inverse LT --> x(t)


clear
close all
clc

% Symbolic variabls
syms s t X 

% System parameters
  m = 1;
  b = 0;
  k = 1;
  x0 = 1;
  v0 = 0;

% Input signal (driving force)
%  y = cos(2*t);
   y = 4*exp(2*t); 
% LT of y(t) --> Y(s)
  Y = laplace(y,t,s)

% % LT of ODE -->  
% %  Z = (m*s^2 + b*s + k)*X - m*(s*x0 - v0) - b*x0 - Y
%    Z = (s^2 - 3*s +2)*X + 3*s - 14 - Y
%   Sol = solve(Z, X)
%   
%   
% 
% sol = ilaplace(Sol,s,t)


Z = (s^2 - 3*s + 2)*Y + 3*s -14 - Y 

%Y1 = s * Y + 3
%Y2 = s * Y1 - 5
%Sol = solve(Y2 - 3 * Y1 + 2 * Y - F, Y)
Sol = solve(Z, Y)

sol = ilaplace(Sol,s,t)
pretty(sol)

pretty(sol)
ezplot(sol,[0,10])
grid on
title('GRAPHIC DISPLAY OF SOLUTION OF ODE BY LAPLACE TRANSFORM')
xlabel('time'),ylabel('f(t)')
legend('laplace transform')


%%
syms s t Y
f = 2
Y = laplace(f,t,s)

Z = (s^2+1)*Y-1 - Y 

%Y1 = s * Y + 3
%Y2 = s * Y1 - 5
%Sol = solve(Y2 - 3 * Y1 + 2 * Y - F, Y)
Sol = solve(Z, Y)

sol = ilaplace(Sol,s,t)
pretty(sol)
ezplot(sol,[0,20])
grid on
title('GRAPHIC DISPLAY OF SOLUTION OF ODE BY LAPLACE TRANSFORM')
xlabel('time'),ylabel('f(t)')
legend('laplace transform')

