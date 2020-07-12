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
    N = 2001;     % must be odd
% xA   xB  lower and upper bounds of region to be rotated
   xA = -3;
   xB = 3;

% x1   x2  limits for plotting function y = f(x)
   x1 = -3/2;
   x2 = 3/2;

% axis of rotation
   yR = 0;
% x values for region   /   xf values for function 
   x  = linspace(xA,xB,N);
   xF = linspace(x1,x2,N);
  
% Region y   /   Function  yF
  % y  = 2 .* sqrt(x);
   %yF = 2 .* sqrt(xF);
   %k = 5;
   %y  = 2+(x.^4 - k .* x.^2)/4;
   %yF = 2+(xF.^4 - k .* xF.^2)/4;
   %y = k .* sin(2*pi*x/2);
  y =  5/2 + sqrt(9/4 - xF.^2);
  yF = 5/2 - sqrt(9/4 - xF.^2);
  %y = ones(1,N) .* 2;
  %yF = y;
   
   yA = y(1); yB = y(end);
   
   
% Volume calculation by disk method
%fn = (y + abs(yR)).^2;
fn = (ones(1,N) .* abs(yR)).^2;
a = xA; b = xB;
vol_pie = simpson1d(fn,a,b);

disp('volume/pie');
disp(vol_pie);


% Setup for animated gif
ag_name = 'ag_volume.gif';   % file name for animated gif
delay = 1;                   % A scalar value 0 to 655 (inclusive), which 
                             % specifies the delay in seconds before
                             % displaying the next image.
   

% GRAPHICS ===============================================================

figure(1)    % Region   /    Function   ---------------------------------
    col = 'b';
    fs = 9;
    set(gcf,'units','normalized','position',[0.1 0.7 0.22 0.22]);
    tx = 'x';
    ty = 'y';
    
    hold on
    
    %xP = x; yP = y;
    %h_area = area(xP,yP);
    %set(h_area,'FaceColor',[0.8 0.8 1]);
    
    xP = xF; yP = y;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx); ylabel(ty);
    
    hold on
    
    col = 'r';
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx); ylabel(ty);
    title('Function y = f(x)');  
    
% axis of rotation
    xP = [xA xB]; yP = [yR yR];
    col = 'r';
    plot(xP,yP,col,'lineWidth',2);
    text(-1,-0.5,'axis of rotation','Color','r');
   
%set(gca,'Xtick',xValues);
   grid on;
   box on;
   axis equal
   set(gca,'Ylim',[-1 5]);
   set(gca,'XLim',[-3 3]);
   set(gca,'Xtick',[-3:3]);
   
figure(2)    % Region   /    Function   ---------------------------------
    col = 'b';
    fs = 9;
    set(gcf,'units','normalized','position',[0.35 0.7 0.22 0.22]);
    tx = 'x';
    ty = 'y';
        hold on
        
% region - reflection
    xP = xF; yP = y;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 1]);
    hold on
    
% function
    xP = xF; yP = y;
    plot(xP,yP,col,'lineWidth',2);    

    
% axis of rotation
    xP = [xA xB]; yP = [yR yR];
    col = 'r';
    plot(xP,yP,col,'lineWidth',2);
    %text(1,0.8,'axis of rotation','color','r');

    xlabel('x');
    ylabel('y'); 
    title('Rotated Region 1');
    
 grid on;
   box on;
   axis equal
   set(gca,'Ylim',[-1 5]);
   set(gca,'XLim',[-3 3]);
   set(gca,'Xtick',[-3:3]);

% ========================================================================

figure(3)    % Region   /    Function   ---------------------------------
    col = 'r';
    fs = 9;
    set(gcf,'units','normalized','position',[0.65 0.7 0.22 0.22]);
    tx = 'x';
    ty = 'y';
        hold on
        
% region - reflection
    xP = xF; yP = yF;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.6 0.6]);
    hold on
    
% function
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);    
   
% axis of rotation
    xP = [xA xB]; yP = [yR yR];
    col = 'r';
    plot(xP,yP,col,'lineWidth',2);
    %text(1,0.8,'axis of rotation','color','r');

    xlabel('x');
    ylabel('y'); 
    title('Rotated Region 2');
    
 grid on;
   box on;
   axis equal
   set(gca,'Ylim',[-1 5]);
   set(gca,'XLim',[-3 3]);
   set(gca,'Xtick',[-3:3]);   

% ========================================================================  
    
figure(6)   % [3D] plot -------------------------------------------------
set(gcf,'units','normalized','position',[0.35 0.4 0.22 0.22]);

[X,Y,Z] = cylinder(y-yR,100);
z = xA + Z.*(xB-xA);
surf(z,Y+yR,X);
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
box on
view(26,24);
axis off
title('outside surface');
set(gca,'Xlim',[-4 4]);
set(gca,'Ylim',[-4 4]);
set(gca,'Zlim',[-4 4]);


% =======================================================================

figure(7)   % [3D] plot -------------------------------------------------
set(gcf,'units','normalized','position',[0.65 0.4 0.22 0.22]);

[X,Y,Z] = cylinder(yF-yR,100);
z = xA + Z.*(xB-xA);
surf(z,Y+yR,X);
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
box on
view(26,24);
axis off
title('inside surface');
set(gca,'Xlim',[-4 4]);
set(gca,'Ylim',[-4 4]);
set(gca,'Zlim',[-4 4]);

% =======================================================================

figure(4)   % animation of the rotation ---------------------------------
set(gcf,'units','normalized','position',[0.1 0.4 0.22 0.22]);

% OUTTER SURFACE
% X Y Z axes
   col = 'r';  LW = 2;  Nt = 100; dx = 150;
   plot3([-3 3], [0 0], [0 0],col,'linewidth',LW);
   hold on
   col = 'k'; LW = 1;
   plot3([0 0],[-4 4],[0 0],col,'linewidth',LW);
   plot3([0 0], [0 0], [-4 4],col,'linewidth',LW);
   xlabel('x'); ylabel('y'); zlabel('z') 
   hold on
   

 % outter function y
     col = 'b'; LW = 1;
 for c = 1 : dx : N
   t = linspace(0, 2*pi, Nt);
   xP = xF(c) .* ones(1,length(t));
   yP = y(c)  .* cos(t);
   zP = y(c)  .* sin(t);
   plot3(xP,yP,zP,col,'linewidth',LW);

 end 

set(gca,'Xlim',[-4 4]);
set(gca,'Ylim',[-4 4]);
set(gca,'Zlim',[-4 4]);



% INNER SURFACE

 % inner function yF
     col = 'm'; LW = 1;
 for c = 1 : dx : N
   t = linspace(0, 2*pi, Nt);
   xP = xF(c) .* ones(1,length(t));
   yP = yF(c)  .* cos(t);
   zP = yF(c)  .* sin(t);
   plot3(xP,yP,zP,col,'linewidth',LW);

 end 

set(gca,'Xlim',[-4 4]);
set(gca,'Ylim',[-4 4]);
set(gca,'Zlim',[-4 4]);

grid on
box on
view(26, 24);

