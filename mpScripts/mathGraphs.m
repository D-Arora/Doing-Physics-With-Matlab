% mathGraphs.m

% mscript for drawing a simple XY graphs
clc
clear all
close all

% INPUTS =================================================================

% Data
   xMin = 0;
   xMax = 0.6;
   N = 2000;
   x = linspace(xMin,xMax,N);
   limX = [0 1];
   tickX = 0:0.2:1;
% Function   
   y = 135 .* x;
   tickY = [0:20:140];
   limY = [0 150];
% X-axis and Y-axis labels  
   xL1 = 'draw  d  [ m ]';
   yL1 = 'force  F [ N ]';
   %5xL2 = 'cos^2(\theta)';
   %yL2 = yL1;
   
%  FIGURES  ==============================================================

figure(1)
   set(gcf,'units','normalized','position',[0.01 0.2 0.23 0.32]);
   fs = 14;
   
   xP = x;
   yP = y;
   
   plot(xP,yP,'b','lineWidth',1);
   set(gca,'xLim',limX);
   set(gca,'xTick',tickX);
   set(gca,'yLim',limY);
   set(gca,'yTick',tickY);
   grid on
   set(gca,'Fontsize',fs);
   xlabel(xL1);
   ylabel(yL1);
   
figure(2)
   set(gcf,'units','normalized','position',[0.01 0.2 0.23 0.32]);
   fs = 14;
   
   xP = linspace(xMin,xMax,N);
   yP = 10 .* cosd(xP).^2;
   xP = cosd(xP).^2;
   plot(xP,yP,'lineWidth',2);
   set(gca,'xTick',-1:0.25:1);
   grid on
   set(gca,'Fontsize',fs);
   xlabel(xL2);
   ylabel(yL2);   
   