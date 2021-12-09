% climateC.m

clc
clear
close all



N = 999;
A = 50;
T = 100;
v = 0.25;

t = linspace(0,3*T,N);
w = 2*pi/T;

% xC = A.*cos(w*t - pi/2); 
% yC = A.*sin(w*t - pi/2);


xC = A.*sin(w*t); 
yC = -A.*cos(w*t);

xB = zeros(1,N); yB =  v.* t - A;

xBC = xB - xC; yBC = yB - yC;

vBCx =  -(A*w).*cos(w*t);
vBCy = v + (A*w).*sin(w*t);
vA = atan2(vBCy,vBCx);

aBCx = (A*w^2).*sin(w*t);
aBCy = (A*w^2).*cos(w*t);
aA = atan2(aBCy,aBCx);



% GRAPHICS *************************************************************

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.05 0.05 0.35 0.8500])
   set(gcf,'color','w');
   FS = 12;

subplot(3,2,1)
   plot(xC,yC,'b','linewidth',1)
   hold on
   plot(xB,yB,'k','linewidth',1)
   plot(xC(end),yC(end),'bo','markerfacecolor','b','MarkerSize',8)
   plot(xB(end),yB(end),'ko','markerfacecolor','k','MarkerSize',8)
   xlabel('x'); ylabel('y')
   xlim([-52 52]); ylim([-52 52])
  set(gca,'XTick',-50:25:50)
  set(gca,'yTick',-50:25:50)
  grid on
  axis square
  set(gca,'fontsize',FS)

subplot(3,2,2)
   plot(xBC,yBC,'b','linewidth',2)
   hold on
   plot(xBC(end),yBC(end),'bo','markerfacecolor','b','MarkerSize',8)
   xlim([-52 52]); ylim([-102 102])
   xticks(-50:25:50)
   yticks(-100:25:100)
   grid on
   axis square
   xlabel('x_{BC}'); ylabel('y_{BC}')
   set(gca,'fontsize',FS)

subplot(3,2,3)
  plot(t,xBC,'b','linewidth',2)
   hold on
   plot(t(end),xBC(end),'bo','markerfacecolor','b','MarkerSize',8)
   grid on
  ylabel('x_{BC}'); xlabel('t')
  set(gca,'fontsize',FS)

subplot(3,2,4)
   plot(t,yBC,'b','linewidth',2)
   hold on
   plot(t(end),yBC(end),'bo','markerfacecolor','b','MarkerSize',8)
   grid on
   ylabel('y_{BC}'); xlabel('t')
   set(gca,'fontsize',FS)

subplot(3,2,5)
 plot(t,(vA - aA)./pi,'b','linewidth',2')
 hold on
 grid on
 yticks(-1.5:0.5:1)
 ylabel('\theta'); xlabel('t')
 set(gca,'fontsize',FS)
%  plot(vA(end)/pi,aA(end)/pi,'ro','markerfacecolor','r','MarkerSize',8)
%  plot(vA(500)/pi,aA(500)/pi,'go','markerfacecolor','g','MarkerSize',8)

subplot(3,2,6)
   %vX = xC(300) + 50*cos(vA(300));
   %vY = yC(300) + 50*sin(vA(300));
   % plot(vX,vY,'bo','markerfacecolor','b','MarkerSize',8)
   plot(xC,yC,'b','linewidth',1)
   hold on
   R = 30;
   magV = R; angleV = vA(1); L = 10; W = 5; LW = 1; 

   for c = [250 375 500 625 750 875 999]
       plot(xC(c),yC(c),'bo','markerfacecolor','b','MarkerSize',8)
     % xP = [xC(c), xC(c) + R*cos(vA(c))];
     % yP = [yC(c), yC(c) + R*sin(vA(c))];
      zT = xC(c) + 1i*yC(c);
      angleV = vA(c); col = [1 0 0];
      DrawArrow(zT,magV,angleV,L,W,LW,col)
      angleV = aA(c); col = [1 0 1];
      DrawArrow(zT,magV,angleV,L,W,LW,col)
     % plot(xP,yP,'r','linewidth',2)
     % xP = [xC(c), xC(c) + R*cos(aA(c))];
     % yP = [yC(c), yC(c) + R*sin(aA(c))];
     % plot(xP,yP,'m','linewidth',2)
   end

   xlim([-75 75]); ylim([-75 75])
   xticks(-50:25:50)
   yticks(-50:25:50)
   grid on
   axis square
   set(gca,'fontsize',FS)
