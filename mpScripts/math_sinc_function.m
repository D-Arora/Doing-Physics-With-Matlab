% math_sinc_function.m

% calls  turningPoints.m
% calls  simpson1d.m

% UNNORMALIZED SINC FUNCTION
% 05 sep 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% ../mphome.htm


%% SINC FUNCTION  sin(x) / x **********************************************
% plot y = sinc(x)
clear all
close all
clc

xMin = -40;
xMax = 40;
N = 999;

x = linspace(xMin, xMax, N);

y = sin(x+eps) ./ (x + eps);

figure(1)
fs = 12;
plot(x,y,'linewidth',2)
grid on
xlabel('x','fontsize',fs);
ylabel('sinc(x)','fontsize',fs);
set(gca,'fontsize',fs);

%%   plot y = sin(x) against x / pi

clear all
close all
clc

xMin = -20;                  % range for x values
xMax = 20;
N = 999;                     % number of partitons

x = linspace(xMin, xMax, N);     % x values
 
y = sin(x+eps) ./(x+eps);         % y values   sinc(x)

yC = cos(x);                      % cosine function

dy_dx = gradient(y);              % gradient of sinc function

% turning points: zero crossings max and min
x_pi = x./pi;
xData = x_pi; 
yData = y .* y;
[indexMin indexMax] = turningPoints(xData, yData) 
disp('y = 0 at x/pi')
x_pi_zeros = x_pi(indexMin)

disp('  ');
disp('x_pi = max / min')
x_m = x_pi(indexMax)
disp('  ');
disp('y = max / min')
y_m = y(indexMax)
disp('  ')
disp('y*y = max / min')
y_m = y(indexMax).^2 
disp('  ')


figure(1)     % sinc
fs = 12;
xP = x./pi; 
yP = y;

plot(xP,yP,'linewidth',2)
hold on
yP = yC;
plot(xP,yP,'m');
grid on
xlabel('x / \pi','fontsize',fs);
ylabel('sinc(x)','fontsize',fs);
set(gca,'fontsize',fs);
set(gca,'Xtick',[-6 : 6]);
set(gca,'Xlim',[-6.5 6.5]);

figure(2)    %  sinc squared
fs = 12;

xP = x_pi;
yP = y.^2;
plot(xP,yP,'linewidth',2);
grid on
xlabel('x / \pi','fontsize',fs);
ylabel('(sinc(x))^2','fontsize',fs);
set(gca,'fontsize',fs);
set(gca,'Xtick',[-6 : 6]);
set(gca,'Xlim',[-6.5 6.5]);

figure(3)    % gradient
fs = 12;

xP = x_pi;
yP = dy_dx;
plot(xP,yP,'linewidth',2);
grid on
xlabel('x / \pi','fontsize',fs);
ylabel('d(sinc(x)) / dx','fontsize',fs);
set(gca,'fontsize',fs);
set(gca,'Xtick',[-6 : 6]);
set(gca,'Xlim',[-6.5 6.5]);

%%  integral of sinc and (sinc)^2

clear all
close all
clc

xMin = -200;                  % range for x values
xMax = 200;
N = 999;                     % number of partitons   must be ODD

x = linspace(xMin, xMax, N);     % x values
 
y = sin(x+eps) ./(x+eps);         % y values   sinc(x)


integral = simpson1d(y,xMin,xMax)/pi


