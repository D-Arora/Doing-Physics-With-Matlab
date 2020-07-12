% m31_simulations.m



%%  variation of v with FT
clear all
close all
clc

L = 0.4;
m = 3e-3;
mu = m/L;

FT = linspace(0,1000,500);

v = sqrt(FT ./mu);

figure(1)
   pos = [0.1 0.2 0.30 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   
   xP = FT;  yP = v;
   plot(xP, yP,'b','linewidth',2);
   
   xlabel('F_T  [ N ] ');
   ylabel('v  [ m.s^{-1} ]');
   box on 
   grid on
   set(gca, 'fontsize',14);
   
   hold on
   
   c = find(FT>799, 1 );
   xP = [FT(c), FT(c)];
   yP = [0, v(c)];
   plot(xP,yP,'r');
   xP = [0 FT(c)]; yP = [v(c), v(c)];
   plot(xP,yP,'r');
   
   
   t1 = '\mu = ';
   t2 = num2str(mu,'%2.2e \n');
   t3 = '  kg.m^{-1}';
   txy = [t1 t2 t3];
   hText = text(400,370,txy);
   set(hText,'fontsize',14);
   
   t1 = 'v = ';
   t2 = num2str(v(c),'%2.0f \n');
   t3 = '  m.s^{-1}';
   txy = [t1 t2 t3];
   hText = text(40,280,txy);
   set(hText,'fontsize',14);
   
   t1 = 'F_T = ';
   t2 = num2str(FT(c),'%2.0f \n');
   t3 = '  N';
   txy = [t1 t2 t3];
   hText = text(580,140,txy);
   set(hText,'fontsize',14);

   
%%  pa31_graphs  Q1  sine function plot ==================================
close all
clear all
clc

figure(1)
   pos = [0.1 0.2 0.30 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'Color','w');

   A = 7;
   wL = 60;
   k = 2*pi/wL;
   x = linspace(0,120,1000);
   y = A .* sin(k*x);
   
   plot(x,y,'b','linewidth',2);
   hold on
   y = A .* sin(k*x+pi);
  % plot(x,y,'r','linewidth',2);
   
   axis([0 120 -11 11]);
   set(gca,'xTick',0:20:120);
   box on
   grid on
   set(gca,'fontsize',14); 
   xlabel('x  [ mm ]');
   ylabel('s  [ mm ]');
  % legend('t','t+T/2');

  
 %%   pa31_graphs    Question 2
close all
clear all
clc

figure(1)
   pos = [0.1 0.2 0.30 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'Color','w');

   N = 1000;
   A = [20 30 10];
   wL = [20 60 30];
   k = 2*pi./wL;
   x = linspace(0,120,N);
   col(1) = 'b';
   col(2) = 'r';
   col(3) = 'g';
   yR = zeros(1,N);
   
   hold on
   for c = 1:3
      y = A(c) .* sin(k(c)*x);
      yR = yR + y;
      plot(x,y,'color',col(c),'linewidth',2);
   end
     plot(x,yR,'k','linewidth',3); 
   axis([0 60 -51 51]);
   set(gca,'xTick',0:10:120);
   set(gca,'yTick',-50:10:50);
   box on
   grid on
   set(gca,'fontsize',14);
   xlabel('x  [ mm ]');
   ylabel('s  [ mm ]');
  % legend('t','t+T/2');
  

  
%%  pa31_graphs.doxc   Question 3
clear all
close all
clc

x = 0:10:120;
y = zeros(1,length(x));

figure(1)
   pos = [0.1 0.2 0.30 0.1];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'Color','w');
   hold on
   
   hPlot = plot(x,y,'bo');
   set(hPlot,'markersize',8,'markerfacecolor','b');

   x = [0 17 27 30 33 43 60 77 87 90 93 103 120];
   y = ones(1,length(x));
   hPlot = plot(x,y,'ro');
   set(hPlot,'markersize',8,'markerfacecolor','r');
   set(hPlot,'markeredgecolor','r');
   axis off
   

%%
clear all
close all
clc


t = linspace(0,120,1000);
T = [50 25 12.5];
A = [10 5 2.5];

w = 2*pi ./ T;

y1 = A(1) .* sin(w(1)*t);
y2 = A(2) .* sin(w(2)*t);
y3 = A(3) .* sin(w(3)*t);

figure(1)
   pos = [0.1 0.2 0.30 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'Color','w');
   
   xP = t;  yP = y1;
   plot(xP, yP,'r','linewidth',1);
   hold on
   yP = y2;
   plot(xP, yP,'k','linewidth',1);
   yP = y3;
   plot(xP, yP,'m','linewidth',1);
   yP = y1+y2+y3;
   plot(xP, yP,'b','linewidth',2);
   
   
   xlabel('time  t  [s]');
   ylabel('wavefunction  s  [a.u.]');
   box on 
   grid on
   set(gca, 'fontsize',14);
   
   %% variation of v with mu
  clear all
close all
clc

L = 0.4;
FT = 800;

mu = linspace(5,200,500) .* 1e-3;

v = sqrt(FT ./mu);

figure(1)
   pos = [0.1 0.2 0.30 0.28];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   
   xP = mu;  yP = v;
   plot(xP, yP,'b','linewidth',2);
   
   xlabel('\mu  [ kg.m^{-1} ] ');
   ylabel('v  [ m.s^{-1} ]');
   box on 
   grid on
   set(gca, 'fontsize',14);
   
   hold on
   
%    c = find(FT>799, 1 );
%    xP = [FT(c), FT(c)];
%    yP = [0, v(c)];
%    plot(xP,yP,'r');
%    xP = [0 FT(c)]; yP = [v(c), v(c)];
%    plot(xP,yP,'r');
%    
%    
   t1 = 'F_T = ';
   t2 = num2str(FT,'%2.0f \n');
   t3 = '  N';
   txy = [t1 t2 t3];
   hText = text(0.100,250,txy);
   set(hText,'fontsize',14);
%    
%    t1 = 'v = ';
%    t2 = num2str(v(c),'%2.0f \n');
%    t3 = '  m.s^{-1}';
%    txy = [t1 t2 t3];
%    hText = text(40,280,txy);
%    set(hText,'fontsize',14);
%    
%    t1 = 'F_T = ';
%    t2 = num2str(FT(c),'%2.0f \n');
%    t3 = '  N';
%    txy = [t1 t2 t3];
%    hText = text(580,140,txy);
%    set(hText,'fontsize',14);


%%
% pa31_graphs.docx   transverse wave
clear all
close all
clc


figure(1)
   pos = [0.1 0.1 0.50 0.8];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'Color','w');
   N = 1000;
   dy = 2;
   x = linspace(0,100,N);
   wL = 20; k = 2*pi/wL;
   T = 80; w = 2*pi/T; t = 0;
         
   hold on
   for c = 0 : 8
      t = c * T/8; 
      y = c *dy + sin(k*x - w*t);
      plot(x,y,'k','linewidth',2);
      xP = x(51); yP = y(51);
      hPlot = plot(xP,yP,'bo');
      set(hPlot,'markersize',8,'markerfacecolor','b');
      xP = [0 x(end)]; yP = [c*dy c*dy];
      plot(xP,yP,'k');    
   end
   grid on
   axis off
   
%%
% pa31_graphs.docx   longitudinal wave
clear all
close all
clc


figure(1)
   pos = [0.1 0.1 0.50 0.8];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'Color','w');
   N = 90;
   dy = 2;
   x0 = linspace(0,90,N);
   wL = 30; k = 2*pi/wL;
   T = 80; w = 2*pi/T; 
   y = zeros(N,1); 
   A = 2.5;
   
   hold on
   for c = 0 : 8
      t = c * T/8; 
      y = y + dy;
      x = x0 + A .*sin(k*x0 - w*t);
      for cc = 1 : N
       xP = [x(cc) x(cc)]; yP = [y-dy/4 y+dy/4];
       plot(xP,yP,'k');
      end
      %set(hPlot,'markersize',4,'markerfacecolor','b');
      % set(hPlot,'markeredgecolor','b');
     
      xP = x(17); yP = y(17);
      hPlot = plot(xP,yP,'ro');
      set(hPlot,'markersize',8,'markerfacecolor','r');
      set(hPlot,'markeredgecolor','r');
   end
   grid on
   axis off
   
 