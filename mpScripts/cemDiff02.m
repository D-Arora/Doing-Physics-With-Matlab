% cemDiff01.m
% 19 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% GRADIENT FUNCTION  of f(x,y)
%   meshgrid / [3D] plots: contourf and quiver

clear all
close all
clc

N = 101;
xMax = 2;
xMin = -xMax;
yMax= xMax;
yMin = -yMax;
x = linspace(xMin, xMax,N);
y = linspace(yMin, yMax,N);
dx = x(2)-x(1);   dy = y(2)-y(1);
[xx, yy] = meshgrid(x,y);

f = xx .* exp(-xx.^2 - yy.^2);
%f = 10 .* (2.*xx.*yy - 3.*xx.^2 -4.*yy.^2 - 18.*xx + 28.*yy +12); 

[delx, dely] = gradient(f,dx,dy);

figure(1);
set(gcf,'units','normalized','position',[0.2 0.2 0.3 0.5]);
surf(xx,yy,f);
shading interp;
xlabel('x  [km]'); ylabel('y  [km]'); zlabel('h  [m]');
view(10, 30);
colormap(summer)
set(gca,'fontsize',14)
axis square

figure(2) % -------------------------------------------------------------
set(gcf,'units','normalized','position',[0.6 0.2 0.35 0.4]);
contourf(xx, yy, f, 16);
shading interp
hold on
d = 1:10:N;
p1 = xx(d,d); p2 = yy(d,d); p3 = delx(d,d); p4 = dely(d,d);
h = quiver(p1, p2, p3, p4);
set(h,'color',[0 0 1],'linewidth', 2);
axis equal
set(gca,'xLim',[xMin, xMax]);
set(gca,'yLim',[yMin, yMax]);
colormap(summer)
xlabel('x  [km]'); ylabel('y  [km]');
set(gca,'fontsize',14)