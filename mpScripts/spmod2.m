% spMod2.m

%%  P20008
clear all
close all
clc

figure(1)   
   pos = [0.1 0.1 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = [0 3 5]; yP = [0 6 6];
   plot(xP,yP,'b','linewidth',2);
   set(gca,'xLim',[0 5]);
   set(gca,'yLim',[0 8]);
   set(gca,'xTick',0:5);
   set(gca,'yTick',0:8);
   ylabel('s  [ m ]');
   xlabel('t  [ s ]'); 
   set(gca,'fontsize',14)
   grid on
   
%%  P20007
clear all
close all
clc

figure(1)   
   pos = [0.1 0.1 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = [0 6 6 12]; yP = [10 6 4 0];
   plot(xP,yP,'b','linewidth',2);
   set(gca,'xLim',[0 13]);
   set(gca,'yLim',[0 10]);
   set(gca,'xTick',0:13);
   set(gca,'yTick',0:10);
   ylabel('p (3 kg object) [ kg.m.s^{-1} ]');
   xlabel('t  [ s ]');
   set(gca,'fontsize',14)
   grid on
   
   %%  P20011
   m = 1000;
   dt = 4.0;
   v1 = 30*(-cosd(45)-1i*sind(45));
   v2 = 10*(sind(30)-1i*cosd(30));
   p1 = m * v1;
   p2 = m * v2;
   dp = p2-p1;
   Favg = dp/dt;
   F = abs(Favg)
   theta = rad2deg(angle(Favg))