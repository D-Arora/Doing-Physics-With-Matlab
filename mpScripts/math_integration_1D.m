% math_integration_1D.m

% One-dimensioanl integration Methods

% Users functiond
%     simpson_1d.m

% 24 aug 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm


clear all
close all
clc
format short

N = 99;                       % number of partitions

a = 0;                       % lower bound
b = pi/2;                    % uper bound

x = linspace(a,b,N);         % partion of the variable x

y1 = cos(x);                 % test integrand 1
y2 = exp(1i.*x);             % test integrand 2

% Closed rectangle rule -------------------------------------------------
tic
h = (b-a) / (N-1);           % dx x(2) - x(1) = h
disp('Closed rectangle rule')
Integral_1 = sum(y1) * h    % estimate of the integral
Integral_2 = sum(y2) * h    % estimate of the integral
disp(toc)
disp('  ');

% Open Midpoint Rule -------------------------------------------------
tic
disp('Open Midpoint rule')
clear x;
h = (b - a) / N;
c = 1 : (2*N-1);
x = a + (h/2) .* c;         % dx = x(2) - x(1) = h
Integral_1 = sum(y1) * h    % estimate of the integral
Integral_2 = sum(y2) * h    % estimate of the integral
disp(toc)
disp('  ');

% Trapezoidal Rule -------------------------------------------------------
tic
disp('Trapezoidal rule')
x = linspace(a,b,N);         % partion of the variable x
Integral_1 = trapz(x,y1)     % estimate of the integral
Integral_2 = trapz(x,y2)     % estimate of the integral
disp(toc)
disp('  ');

% Simpson's 1/3 rule -----------------------------------------------------
tic
disp('Simpson 1/3 rule')
Integral_1 = simpson1d(y1,a,b)     % estimate of the integral
Integral_2 = simpson1d(y2,a,b)     % estimate of the integral
disp(toc)
disp('  ');

% Graphics --------------------------------------------------------------
% figure(1)
% subplot(2,1,1)
% plot(x,y1)
% 
% subplot(2,1,2);
% plot(x,real(y2))
% hold on
% plot(x,imag(y2))

