%m5_expt533.m
clear all
close all
clc

% XY data
xP = [0 0.50 1.02 1.45 2.00 2.48 2.94 3.41 3.87 4.42 4.82 5.33 5.82 6.25 6.75 7.21 7.68 8.11 8.54];
yP = [0 0.80 1.49 2.16 2.72 3.21 3.59 4.02 4.26 4.36 4.58 4.61 4.61 4.52 4.37 4.15 3.81 3.47 3.14];
sf = 1/11.65;
xP = xP.*sf;
yP = yP.*sf;
tS = 0:length(xP)-1;

tf = 0.0251;
t = tS*tf;
u= 2.822;


% GRAPHICS ===============================================================

figure(1)   % -----------------------------------------------------------
   fs = 14;
   set(gcf,'Units','Normalized');
   set(gcf,'Position',[0.02 0.04 0.3 0.25]);
   hPlot = plot(xP,yP,'mo');
   set(hPlot,'markerfacecolor','m');
   set(gca,'fontSize',fs);
   xlabel('s_x   [ m ]');
   ylabel('s_y   [ m ]');
   grid on
   axis equal
   axis([0 0.8 0 0.5])

figure(2)   % ------------------------------------------------------------
  fs = 14;
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.02 0.4 0.3 0.25]);
  hPlot = plot(t,yP,'bo');
  set(hPlot,'markerfacecolor','b');
  set(gca,'fontSize',fs);
  xlabel('time t  [ s ]');
  ylabel('s_y   [ m ]');
  grid on
  axis([0 0.5 0 0.5])


  
figure(3)   % ------------------------------------------------------------
  fs = 14;
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.32 0.04 0.3 0.25]);
   X = tS(2:end);
   Y = yP(2:end)./X;
   hPlot = plot(X,Y,'bo');
   set(hPlot,'markerfacecolor','b');
   set(gca,'fontSize',fs);
   xlabel('time steps  t_S  [ steps ]');
   ylabel('s_y / t_S  [ m/steps ]');
   grid on
   
figure(4)   % ------------------------------------------------------------
  fs = 14;
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.32 0.4 0.3 0.25]);
   X = tS(2:end).*tf;
   Y = yP(2:end)./X;
   hPlot = plot(X,Y,'bo');
   set(hPlot,'markerfacecolor','b');
   set(gca,'fontSize',fs);
   xlabel('time t  [ s]');
   ylabel('s_y / t  [ m/s ]');
   grid on

figure(5)   % ------------------------------------------------------------
  fs = 14;
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.65 0.04 0.3 0.25]);
   X = tS(2:end).*tf;
   Y = 2.*yP(2:end)./X - u;
   
   hPlot = plot(X,Y,'bo');
   set(hPlot,'markerfacecolor','b');
   set(gca,'fontSize',fs);
   xlabel('time  t  [ s ]');
   ylabel('v_y  [ m/s ]');
   grid on

figure(6)   % ------------------------------------------------------------
  fs = 14;
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.65 0.4 0.3 0.25]);
   X = tS(2:end).*tf;
   Y = xP(2:end);
   hPlot = plot(X,Y,'ro');
   set(hPlot,'markerfacecolor','r');
   set(gca,'fontSize',fs);
   xlabel('time  t  [ s ]');
   ylabel('x  [ m ]');
   grid on
 