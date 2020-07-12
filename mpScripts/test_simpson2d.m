clear all 
close all
clc

num = 99;

%%
kx = 3; ky = 3;
x = linspace(0,pi/2,num);
y = x;

[xx yy] = meshgrid(x,y);

f = sin(kx * xx) .* sin(ky * yy)

integral2 = simpson2d(f,x(1),x(end),y(1),y(end))


f1 = sin(kx * x);
integral1 = simpson1d(f1,x(1),x(end))

%%
clear all
close all
clc

num = 9;
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
ty = 'y'
set(gca,'fontsize',fs);
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.4 0.4]);
plot(xx,yy,'o')
axis([-1 3 0 6])
set(gcf,'Xlim',[-1 3]);
axis equal
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);



