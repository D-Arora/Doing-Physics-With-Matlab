% cemDiff03.m
% 18 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% DIVERGENCE OF A [2D] VECTOR FIELD
   

clear all
close all
clc

% INPUTS  ================================================================
% Number of grid points (integer)    N = 101 default value
   N = 101;
% Range for X and Y values  [minX maxX minY maxY]
   XY = [-50 50 -50 50];
   
% CALCULATIONS  ==========================================================   
% X and Y ranges
minX = XY(1); maxX = XY(2); minY  = XY(3); maxY = XY(4);
x = linspace(minX, maxX,N);
y = linspace(minY, maxY,N);
[xx, yy] = meshgrid(x,y);

% Define [2D] vector field V
%   sinusoidal function;
       wLx = 50; wLy = 40;
       Vxx = cos(2*pi*xx/wLx);
       Vyy = cos(2*pi*yy/wLy);
%   vector field
     % Vxx = xx.* yy;
     % Vyy = yy.^2;

% DIVERGENCE
     divV  = divergence(xx, yy, Vxx, Vyy);

% GRAPHICS ===============================================================

figure(1);  %  -----------------------------------------------------------
   set(gcf,'units','normalized','position',[0.2 0.2 0.3 0.5]);
   surf(xx,yy,divV);
   shading interp;
   colormap(summer)
   set(gca,'fontsize',14)
   axis square
   box on
   view(-26,52);
   title('divergence');
   xlabel('x'); ylabel('y'); zlabel('div');
   
figure(2) % -------------------------------------------------------------
   set(gcf,'units','normalized','position',[0.6 0.2 0.35 0.4]);
   pcolor(xx,yy,divV);
   shading interp
   colorbar
   hold on
   d = 1:8:N;
   p1 = xx(d,d); p2 = yy(d,d); p3 = Vxx(d,d); p4 = Vyy(d,d);
   h = quiver(p1, p2, p3, p4);
   set(h,'color',[0 0 1],'linewidth', 2);
   axis equal
   set(gca,'xLim',[minX, maxX]);
   set(gca,'yLim',[minY, maxY]);
   colormap(summer)
   xlabel('x'); ylabel('y');
   set(gca,'fontsize',14)
   title('vector field & divergence')
   box on
   
   