% Laplace05.m

clear
close all
%clc


syms s t Y
f = 4 * exp(2 * t)
F = laplace(f,t,s)

Z = (s^2 - 3*s + 2)*Y + 3*s -14 - F 

%Y1 = s * Y + 3
%Y2 = s * Y1 - 5
%Sol = solve(Y2 - 3 * Y1 + 2 * Y - F, Y)
Sol = solve(Z, Y)

sol = ilaplace(Sol,s,t)
pretty(sol)
ezplot(sol,[0,10])
grid on
title('GRAPHIC DISPLAY OF SOLUTION OF ODE BY LAPLACE TRANSFORM')
xlabel('time'),ylabel('f(t)')
legend('laplace transform')


%%
% syms s t Y
% f = 4*exp(2*t)
% F = laplace(f,t,s)
% 
% Z = (s^2+1)*Y-1 - F 
% 
% %Y1 = s * Y + 3
% %Y2 = s * Y1 - 5
% %Sol = solve(Y2 - 3 * Y1 + 2 * Y - F, Y)
% Sol = solve(Z, Y)
% 
% sol = ilaplace(Sol,s,t)
% pretty(sol)
% ezplot(sol,[0,10])
% grid on
% title('GRAPHIC DISPLAY OF SOLUTION OF ODE BY LAPLACE TRANSFORM')
% xlabel('time'),ylabel('f(t)')
% legend('laplace transform')

