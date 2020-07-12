% math_vol_01.m
% 12 may 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% ../mphome.htm
clear all
close all
clc

% INPUTS =================================================================
% number of partitions
    N = 2000;
% xA   xB  lower and upper bounds of region to be rotated
   xA = -2;
   xB = 2;

% x1   x2  limits for plotting function y = f(x)
   x1 = -2.1;
   x2 = 2.1;

% x values for region   /   xf values for function 
   x  = linspace(xA,xB,N);
   xF = linspace(x1,x2,N);
  
% Region y   /   Function  yF
   k = 5;
   y  = 2+(x.^4 - k .* x.^2)/4;
   yF = 2+(xF.^4 - k .* xF.^2)/4;
   %y = k .* sin(2*pi*x/2);



% GRAPHICS ===============================================================

figure(1)    % Region   /    Function   ---------------------------------
    col = 'b';
    fs = 9;
    set(gcf,'units','normalized','position',[0.1 0.7 0.22 0.22]);
    tx = 'x';
    ty = 'y ';
    
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx); ylabel(ty);
    hold on
   
    xP = x; yP = y;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 0.8]);
   
   set(gca,'Ylim',[0 1.1 * max(y)]);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
    grid on;
    box on;

figure(2)   % [3D] plot -------------------------------------------------
set(gcf,'units','normalized','position',[0.35 0.7 0.22 0.22]);
[X,Y,Z] = cylinder(y,100);
z = xA + Z.*(xB-xA);
surf(z,Y,X);
%axis square
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
box on

figure(3)   % [3D] plot --------------------------------------------------
set(gcf,'units','normalized','position',[0.6 0.7 0.22 0.22]);
[X,Y,Z] = cylinder(y,100);
z = xA + Z.*(xB-xA);
surfl(z,Y,X);
axis square
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
colormap(copper)
hold on
plot3([xA xB],[0 0], [0 0]);
box on
view(-15,26);

figure(4)   % animation of the rotation ---------------------------------
set(gcf,'units','normalized','position',[0.1 0.4 0.22 0.22]);

% X-axis
plot3([xA xB],[0 0], [0 0],'k','linewidth',2);
set(gca,'Xlim',[xA xB]);
set(gca,'Ylim',[xA xB]);
set(gca,'Zlim',[xA xB]);
xlabel('x'); ylabel('y'); zlabel('z') 
hold on
grid on
box on
view(70,10);

% Function
xP = x; yP = y; zP = zeros(1,N);
plot3(xP,yP,zP,'b','linewidth',2);

t = 0 : pi/100 : 2*pi;
xP = x(end) .* ones(1,length(t)); yP = y(end) .* cos(t); zP = y(end).*sin(t);
h_plot3 = plot3(xP,yP,zP,'r','linewidth',1);
xP = x(1) .* ones(1,length(t)); yP = y(1) .* cos(t); zP = y(1).*sin(t);
h_plot3 = plot3(xP,yP,zP,'r','linewidth',1);

for t = 0:pi/6:2*pi;
xP = x; yP = y .* cos(t); zP = y.*sin(t);
h_plot3 = plot3(xP,yP,zP,'r','linewidth',1);
%plot3(xB,yP(end),zP(end),'ro');
pause(1)
set(gca,'Xlim',[xA xB]);
set(gca,'Ylim',[xA xB]);
set(gca,'Zlim',[xA xB]);
end





