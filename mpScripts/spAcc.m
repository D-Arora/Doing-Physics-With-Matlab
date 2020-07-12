% spAcc.m


clear all
close all
clc

% =======================================================================
% Setup for saving images (im) 
% ======================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_freeFall.gif';
% Delay in seconds before displaying the next image  
   delay = 0;  
% Frame counter start
   nt = 1; 

% ========================================================================

g = -9.81;
u = 50;
 
tMax = -2*u/g;
N = 100;

t = linspace(0,tMax,N);
s = u.*t + (0.5*g).*t.^2;
v = u + g.*t;

a = g .* ones(N,1);


for c = 1: N
figure(1)
   fs = 14;
   set(gcf,'units','normalized','position',[0.1 0.1 0.5 0.8]);
   
   xP = t; yP = s;
   pos1 = [0.2 0.7 0.7 0.25];
   subplot('Position',pos1);
   plot(xP,yP,'b');
   set(gca,'fontsize',fs);
   set(gca,'xLim',[0 11]);
   set(gca,'yLim',[0 150]);
   grid on
   hold on
   
   yP = v;
   pos2 = [0.2 0.4 0.7 0.25];
   subplot('Position',pos2);
   plot(xP,yP,'b');
   set(gca,'fontsize',fs);
   grid on
   hold on
   
   yP = a;
   pos3 = [0.2 0.1 0.7 0.25];
   subplot('Position',pos3);
   plot(xP,yP,'b');
   set(gca,'fontsize',fs);
   set(gca,'yLim',[-10.5 0.5]);
   grid on
   xlabel('time  t  [ s ]');
   hold on
   
%    yP = s(1:N/10:N/2); xP = zeros(length(yP),1);
%    pos4 = [0.1 0.2 0.1 0.6];
%    subplot('Position',pos4);
%    plot(xP,yP,'o');
%    set(gca,'fontsize',fs);
%    %set(gca,'yLim',[-10.5 0.5]);
%    grid on
   
   % --------------------------------------------------------------------
   xP = t(1:c); yP = s(1:c);
   pos1 = [0.2 0.7 0.7 0.25];
   subplot('Position',pos1);
   plot(xP,yP,'r','lineWidth',3);
   set(gca,'fontsize',fs);
   set(gca,'xLim',[0 11]);
   set(gca,'yLim',[0 150]);
   set(gca,'xTick',0:1:11);
   grid on
   hold on
   grid on
   ylabel('s [ m ]');
   
   xP = t(1:c); yP = v(1:c);
   pos2 = [0.2 0.4 0.7 0.25];
   subplot('Position',pos2);
   plot(xP,yP,'r','lineWidth',3);
   set(gca,'fontsize',fs);
   grid on
   set(gca,'xLim',[0 11]);
   set(gca,'yLim',[-55 55]);
   set(gca,'yTick',-50:10:50);
   set(gca,'xTick',0:1:11);
   ylabel('v [ m.s^{-1} ]');
   
   xP = t(1:c); yP = a(1:c);
   pos3 = [0.2 0.1 0.7 0.25];
   subplot('Position',pos3);
   plot(xP,yP,'r','lineWidth',3);
   set(gca,'fontsize',fs);
   set(gca,'yLim',[-10.5 0.5]);
   set(gca,'xLim',[0 11]);
   set(gca,'xTick',0:1:11);
   grid on
   xlabel('t  [ s ]');
   ylabel('a [ m.s^{-2} ]');
   
   yP = s(c); xP = 0;
   pos4 = [0.02 0.2 0.1 0.6];
   subplot('Position',pos4);
   hPlot = plot(xP,yP,'o');
   set(gca,'yLim',[-5 150]);
   set(hPlot,'MarkerFacecolor','r','MarkerSize',25)
   hold on
   plot([-0.5 0.5],[0 0],'k');
   hold off
   axis off
   pause(0.001)
   
   if f_gif > 0 
       frame = getframe(1);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
     %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
    end
   
end