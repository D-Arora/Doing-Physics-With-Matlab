clear all 
close all
clc



%%
num = 5;
kx = 3; ky = 3;
x = linspace(0,pi/2,num);
y = x;

[xx yy] = meshgrid(x,y);

f = sin(kx * xx) .* sin(ky * yy)

integral2 = simpson2d(f,x(1),x(end),y(1),y(end))


f1 = sin(kx * x);
integral1 = simpson1d(f1,x(1),x(end))

%%   EXAMPLE 1   f = x^2 *  y^3
clear all
close all
clc

num = 5;
xMin = 0;
xMax = 2;
yMin = 1;
yMax = 5;

x = linspace(xMin,xMax,num);
y = linspace(yMin,yMax,num);

[xx yy] = meshgrid(x,y);
f = xx.^2 .* yy.^3;

Ixy = simpson2d(f,x(1),x(end),y(1),y(end))

figure(1)
fs = 12;
tx = 'x';
ty = 'y';
set(gca,'fontsize',fs);
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.4 0.4]);
plot(xx,yy,'o')
axis([-1 3 0 6])
set(gcf,'Xlim',[-1 3]);
axis equal
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);

%% EXAMPLE 2 - catesian coordinates    
tic
clear all
close all
clc

num = 999;
xMin = -1;
xMax = 1;
yMin = -1;
yMax = 1;

a = 1;

x = linspace(xMin,xMax,num);
y = linspace(yMin,yMax,num);

[xx yy] = meshgrid(x,y);

f = real(sqrt(a^2 - xx.^2 - yy.^2));
f((xx.^2 + yy.^2) > a^2) = 0;

ax = xMin; bx = xMax;
ay = yMin; by = yMax;

format long

Ixy = simpson2d(f,ax,bx,ay,by)

figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.2]);
fs = 12;
tx = 'x';
ty = 'y';
set(gca,'fontsize',fs);
surf(xx,yy,f)
shading interp
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
toc

%% EXAMPLE 3    polar coordinates
clear all
close all
clc

num = 299;
xMin = 0;
xMax = 2;
yMin = 0;
yMax = 2*pi;

a = xMax;
h = 3;

x = linspace(xMin,xMax,num);
y = linspace(yMin,yMax,num);

[xx yy] = meshgrid(x,y);

xP = xx .* cos(yy);                    % cartesian coordinates for Q         
yP = xx .* sin(yy);

f = xx .* xx;

ax = xMin; bx = xMax;
ay = yMin; by = yMax;


format long

Ixy = simpson2d(f,ax,bx,ay,by)

figure(1)
fs = 12;
tx = 'xP';
ty = 'yP';
set(gca,'fontsize',fs);
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.4 0.45]);
plot(xP,yP,'o')
%axis([-1 3 0 6])
%set(gcf,'Xlim',[-1 3]);
axis equal
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);

figure(2)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.2]);
fs = 12;
tx = 'rho';
ty = 'phi';
set(gca,'fontsize',fs);
surf(xx,yy,f)
shading interp
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
title('f(rho,phi)','fontsize',fs);









