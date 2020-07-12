% dydx01.m
% 1 jun 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
clear all
close all
clc

%% y = constant
N = 500;
x1 = 0;
x2 = 10;
x = linspace(x1,x2,N);
y = 10.*ones(1,N);
dydx = zeros(1,N);

figure(1)
    col = 'b';
    fs = 9; lW = 2;
    set(gcf,'units','normalized','position',[0.1 0.3 0.22 0.42]);
    tx = 'x';
    ty = 'y';
    tM = 'function y = f(x)= 10';
    yRange = [-1 11];
    xP = x; yP = y;
    subplot(2,1,1)    
    plot(xP,yP,col,'lineWidth',lW);
    grid on
    box on
    xlabel(tx); ylabel(ty); title(tM);
    set(gca,'Ylim',yRange);
    
    col = 'r';
    tx = 'x';
    ty = 'dy/dx';
    tM = 'derivative dy/dx = 0';
    yRange = [-1 11];
    xP = x; yP = dydx;
    subplot(2,1,2)    
    plot(xP,yP,col,'lineWidth',lW);
    grid on
    box on
    xlabel(tx); ylabel(ty); title(tM);
    set(gca,'Ylim',yRange);
    
%% y = polynomial
clear all
close all
clc
N = 500;
x1 = -1;
x2 = 3;
x = linspace(x1,x2,N);
y = 2*x.^3 -5.*x.^2 - 8.*x + 12;
dydx = 6.*x.^2 - 10.*x -8;

figure(1)
    col = 'b';
    fs = 9; lW = 2;
    set(gcf,'units','normalized','position',[0.1 0.3 0.22 0.42]);
    tx = 'x';
    ty = 'y';
    tM = 'function y = 2x^3 - 5x^2 - 8x + 12';
    yRange = [-1 11];
    xP = x; yP = y;
    subplot(2,1,1)    
    plot(xP,yP,col,'lineWidth',lW);
    grid on
    box on
    xlabel(tx); ylabel(ty); title(tM);
    %set(gca,'Ylim',yRange);
    
    col = 'r';
    tx = 'x';
    ty = 'dy/dx';
    tM = 'derivative dy/dx = 6x^2 - 10x - 8';
    yRange = [-1 11];
    xP = x; yP = dydx;
    subplot(2,1,2)    
    plot(xP,yP,col,'lineWidth',lW);
    grid on
    box on
    xlabel(tx); ylabel(ty); title(tM);
    %set(gca,'Ylim',yRange);

%% y = f(x)
clear all
close all
clc
N = 500;
x1 = 1.2;
x2 = 10;
x = linspace(x1,x2,N);
%u = 3.*x.^6 + 4.*x.^(-0.5);
%v = 3.*x.^2 + 6.*x.^0.5 + 8;
%y = u .* v;
%dydx = 72.*x.^7 + 117.*x.^(11/2) + 144.*x.^5 + 18.*x.^0.5 - 0.*x.^(-1) - 16.*x.^(-3/2);

% u = 2.*x.^2 + 3.*x +5;
% y = u.^(7/2);
% dydx = (7/2).*u.^(5/2).*(4.*x + 3);

u = (2.*x + 1).^(0.5);
v = (2.*x - 1).^(-0.5);
y = u .* v;
dydx = (-(2.*x + 1).*(2.*x -1).^0.5 + (2.*x - 1).^(3/2)) ./ ((2.*x - 1).^2 .* (2.*x + 1).^0.5); 



figure(1)
    col = 'b';
    fs = 9; lW = 2;
    set(gcf,'units','normalized','position',[0.1 0.3 0.22 0.5]);
    tx = 'x';
    ty = 'y';
    tM = 'function y = f(x)';
    yRange = [-1 11];
    xP = x; yP = y;
    subplot(2,1,1)    
    plot(xP,yP,col,'lineWidth',lW);
    grid on
    box on
    xlabel(tx); ylabel(ty); title(tM);
    %set(gca,'Ylim',yRange);
    
    col = 'r';
    tx = 'x';
    ty = 'dy/dx';
    tM = 'derivative dy/dx ';
    yRange = [-1 11];
    xP = x; yP = dydx;
    subplot(2,1,2)    
    plot(xP,yP,col,'lineWidth',lW);
    
    grid on
    box on
    xlabel(tx); ylabel(ty); title(tM);
    %set(gca,'Ylim',yRange);
