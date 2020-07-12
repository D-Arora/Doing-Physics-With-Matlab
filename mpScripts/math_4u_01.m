% math_4u_01.m

clear all
close all
clc

xMin = 0;
xMax= 2.5;
N = 500;

x = linspace(xMin,xMax,N);
y = 3.*x.^3 + 2.*x +1;

x1 = linspace(1,2,N);
y2 = 3.*x1.^3 + 2.*x1 +1;

figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.2 0.3]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x1;   yP = y2;
   h_area = area(xP,yP,'lineWidth',LW);
   set(h_area,'FaceColor',[0.5 0.5 0.8],'EdgeColor',[0.5 0.5 0.8]);
   hold on
   
   xP = x;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
   
   text(1.5,7,'A','fontsize',14);
   title('y = 3x^2 + 2x + 1');
   xlabel(tx); ylabel(ty);
   grid on; box on;
 
   
%%    math_trig.doc   sine function
clear all
close all
clc

tMin = -360;
tMax= 360;
N = 2500;

t = linspace(tMin,tMax,N);
y = sind(t);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
      
   xlabel(tx); ylabel(ty);
   grid on; box on;
   set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',tMin:45:tMax);
   
 %%    math_trig.doc   cosine function
clear all
close all
clc

tMin = -360;
tMax= 360;
N = 2500;

t = linspace(tMin,tMax,N);
y = cosd(t);
t1 = tMin:30:tMax;
y1 = cosd(t1);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
    xP = t1;   yP = y1;
   h_plot = plot(xP,yP,'o');
   set(h_plot,'MarkerFaceColor','b','MarkerEdgeColor','b');
   col = 'k';LW = 1.5;
   xP = [tMin tMax];   yP = [0 0];
   plot(xP,yP,col,'lineWidth',LW);
   
   col = 'k';LW = 0.5;
   xP = [0 0];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   
   col = 'k';LW = 0.5;
   xP = [90 90];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [180 180];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [270 270];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW);
   xP = [360 360];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   
   xP = [-360 360];   yP = [0.5 0.5];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [1/sqrt(2) 1/sqrt(2)];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [sqrt(3)/2 sqrt(3)/2];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-0.5 -0.5];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-1/sqrt(2) -1/sqrt(2)];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-sqrt(3)/2 -sqrt(3)/2];
   plot(xP,yP,col,'lineWidth',LW); 
   
   
   xlabel(tx); ylabel(ty);
   grid on; box on;
   
   text(15,1.2,'1st');
   text(15,1.1,'quadrant');
   text(15+90,1.2,'2nd');
   text(15+90,1.1,'quadrant');
   text(15+180,1.2,'3nd');
   text(15+180,1.1,'quadrant');
   text(15+270,1.2,'4th');
   text(15+270,1.1,'quadrant');
   
   set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',tMin:45:tMax);
   set(gca,'YLim',[-1.2 1.3]);  
   
%%    math_trig.doc   tan function
clear all
close all
clc

tMin = -360;
tMax= 360;
N = 2500;
Ylimit = 20;
t = linspace(tMin,tMax,N);
y = tand(t);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
    col = 'k';LW = 0.5;
   xP = [-360 360];   yP = [0 0];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [0 0];   yP = [-Ylimit Ylimit];
   plot(xP,yP,col,'lineWidth',LW);
   xP = [-360 360];   yP = [1 1];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-1 -1];
   plot(xP,yP,col,'lineWidth',LW); 
   
   %xlabel(tx);
   ylabel(ty);
   grid on; box on;
   set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',tMin:45:tMax);
   
   set(gca,'YLim',[-Ylimit Ylimit]);
   
   %%    math_trig.doc   sine and cosine functions
clear all
close all
clc

tMin = -360;
tMax= 360;
N = 2500;
Ylimit = 1.20;
t = linspace(tMin,tMax,N);
ys = sind(t);
yc = cosd(t);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.5 0.3]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = ys;
   plot(xP,yP,'lineWidth',LW);
   hold on
    col = 'r';LW = 2;
   xP = t;   yP = yc;
   plot(xP,yP,col,'lineWidth',LW); 
   
   col = 'k';LW = 0.5;
   xP = [-360 360];   yP = [0 0];
   plot(xP,yP,col,'lineWidth',LW); 
   
   col = 'k';LW = 0.5;
   xP = [0 0];   yP = [-Ylimit Ylimit];
   plot(xP,yP,col,'lineWidth',LW); 
   
   
   %xlabel(tx);
   ylabel(ty);
   grid on; box on;
   set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',tMin:45:tMax);
   
   set(gca,'YLim',[-Ylimit Ylimit]);
   
   %%    math_trig.doc   sine function
clear all
close all
clc

tMin = -360;
tMax= 360;
N = 2500;

t = linspace(tMin,tMax,N);
y = sind(t);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
      
   xlabel(tx); ylabel(ty);
   grid on; box on;
   set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',tMin:45:tMax);
   
 %%    math_trig.doc   cosine function
clear all
close all
clc

tMin = -360;
tMax= 360;
N = 2500;

t = linspace(tMin,tMax,N);
y = cosd(t);
t1 = tMin:30:tMax;
y1 = cosd(t1);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
    xP = t1;   yP = y1;
   h_plot = plot(xP,yP,'o');
   set(h_plot,'MarkerFaceColor','b','MarkerEdgeColor','b');
   col = 'k';LW = 1.5;
   xP = [tMin tMax];   yP = [0 0];
   plot(xP,yP,col,'lineWidth',LW);
   
   col = 'k';LW = 0.5;
   xP = [0 0];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   
   col = 'k';LW = 0.5;
   xP = [90 90];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [180 180];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [270 270];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW);
   xP = [360 360];   yP = [-1.2 1.2];
   plot(xP,yP,col,'lineWidth',LW); 
   
   xP = [-360 360];   yP = [0.5 0.5];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [1/sqrt(2) 1/sqrt(2)];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [sqrt(3)/2 sqrt(3)/2];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-0.5 -0.5];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-1/sqrt(2) -1/sqrt(2)];
   plot(xP,yP,col,'lineWidth',LW); 
   xP = [-360 360];   yP = [-sqrt(3)/2 -sqrt(3)/2];
   plot(xP,yP,col,'lineWidth',LW); 
   
   
   xlabel(tx); ylabel(ty);
   grid on; box on;
   
   text(15,1.2,'1st');
   text(15,1.1,'quadrant');
   text(15+90,1.2,'2nd');
   text(15+90,1.1,'quadrant');
   text(15+180,1.2,'3nd');
   text(15+180,1.1,'quadrant');
   text(15+270,1.2,'4th');
   text(15+270,1.1,'quadrant');
   
   set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',tMin:45:tMax);
   set(gca,'YLim',[-1.2 1.3]);  
   
%%    math_trig.doc   y = y0+ymax sin(2pit/T + phi) function
clear all
close all
clc

tMin = 0;
tMax= 100;
T = 25;
y0 = 6;
phi = pi/4;
ymax = 4;
N = 2500;
Ylimit = 12;
t = linspace(tMin,tMax,N);

y = y0 + ymax .* sin(2*pi*t/T + phi);

figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t;   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('time  t  [a.u.]');
   ylabel('displacement  y  [a.u.]');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   %set(gca,'Xtick',tMin:45:tMax);
   set(gca,'YLim',[0 Ylimit]);
   
  %%    math_trig.doc   y = sin(theta) function
clear all
close all
clc

tMin = -1*pi;
tMax= 1*pi;

N = 2500;
Ylimit = 1.2;
t = linspace(tMin,tMax,N);

y = sin(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t/(pi);   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('phase angle   \theta  / \pi');
   ylabel('y');
    title('y = sin(\theta)');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   
   hold on
   yP = [-0.5 0.5]; xP = [-1/6 1/6];
   plot(xP,yP,'ro','MarkerFaceColor','r');
    yP = [-0.5 0.5]; xP = [-1+1/6 1-1/6];
   plot(xP,yP,'ro','MarkerFaceColor','r');
  
   %%    math_trig.doc   y = cos(theta) function
clear all
close all
clc

tMin = -1*pi;
tMax= 1*pi;

N = 2500;
Ylimit = 1.2;
t = linspace(tMin,tMax,N);

y = cos(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t/(pi);   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('phase angle   \theta  / \pi');
   ylabel('y');
    title('y = cos(\theta)');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   
   hold on
   yP = [-0.5 0.5]; xP = [1-1/3 1/3];
   plot(xP,yP,'ro','MarkerFaceColor','r');
    yP = [-0.5 0.5]; xP = [-1+1/3 -1/3];
   plot(xP,yP,'ro','MarkerFaceColor','r');
   
     %%    math_trig.doc   y = cosec(theta) function
clear all
close all
clc

tMin = -2*pi;
tMax= 2*pi;

N = 12521;
Ylimit = 20;
t = linspace(tMin,tMax,N);

y = 1 ./sin(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t/(pi);   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('phase angle   \theta  / \pi');
   ylabel('y');
   title('y = cosec(\theta)');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   
   %%    math_trig.doc   y = sec(theta) function
clear all
close all
clc

tMin = -2*pi;
tMax= 2*pi;

N = 12521;
Ylimit = 20;
t = linspace(tMin,tMax,N);

y = 1 ./cos(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t/(pi);   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('phase angle   \theta  / \pi');
   ylabel('y');
   title('y = sec(\theta)');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   
     %%    math_trig.doc   y = cot(theta) function
clear all
close all
clc

tMin = -2*pi;
tMax= 2*pi;

N = 12521;
Ylimit = 20;
t = linspace(tMin,tMax,N);

y = 1 ./tan(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t/(pi);   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('phase angle   \theta  / \pi');
   ylabel('y');
   title('y = cot(\theta)');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
 
  
   
    %%    math_trig.doc   theta = asin(y) function
clear all
close all
clc

tMin = -1*pi;
tMax= 1*pi;

N = 2500;
Ylimit = 1.2;
t = linspace(tMin,tMax,N);

y = sin(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   yP = t/(pi);   xP = y;
   plot(xP,yP,'lineWidth',LW);
  
   ylabel('\theta  / \pi');
   xlabel('y');
    title('\theta = asin(y)');
   grid on; box on;
   set(gca,'Ytick',[-1:0.5:1]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   
   hold on
   xP = [-0.5 0.5]; yP = [-1/6 1/6];
   plot(xP,yP,'ro','MarkerFaceColor','r');
    xP = [-0.5 0.5]; yP = [-1+1/6 1-1/6];
   plot(xP,yP,'ro','MarkerFaceColor','r');
   
   
    %%    math_trig.doc   theta = acos(y) function
clear all
close all
clc

tMin = -1*pi;
tMax= 1*pi;

N = 2500;
Ylimit = 1.2;
t = linspace(tMin,tMax,N);

y = cos(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.25]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   yP = t/(pi);   xP = y;
   plot(xP,yP,'lineWidth',LW);
  
   ylabel('\theta  / \pi');
   xlabel('y');
    title('\theta = acos(y)');
   grid on; box on;
   set(gca,'Ytick',[-1:0.5:1]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   
   hold on
   xP = [-0.5 0.5]; yP = [1-1/3 1/3];
   plot(xP,yP,'ro','MarkerFaceColor','r');
    xP = [-0.5 0.5]; yP = [-(1-1/3) -1/3];
   plot(xP,yP,'ro','MarkerFaceColor','r');
   
   
   
      %%    math_trig.doc   theta = atan(y) function
clear all
close all
clc

tMin = -1*pi;
tMax= 1*pi;

N = 2500;
Ylimit = 1.2;
t = linspace(tMin,tMax,N);

y = tan(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.28]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   yP = t/(pi);   xP = y;
   plot(xP,yP,'lineWidth',LW);
  
   ylabel('\theta  / \pi');
   xlabel('y');
    title('\theta = atan(y)');
   grid on; box on;
   set(gca,'Ytick',[-1:0.5:1]);
  % set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   set(gca,'xLim',[-20 20]); 
%    hold on
%    xP = [-0.5 0.5]; yP = [1-1/3 1/3];
%    plot(xP,yP,'ro','MarkerFaceColor','r');
%     xP = [-0.5 0.5]; yP = [-(1-1/3) -1/3];
%    plot(xP,yP,'ro','MarkerFaceColor','r');
%    
%    

     %%    math_trig.doc   y = tan(theta) function
clear all
close all
clc

tMin = -1*pi;
tMax= 1*pi;

N = 12521;
Ylimit = 25;
t = linspace(tMin,tMax,N);

y = tan(t);
figure(1)
   fs = 10; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.28]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'angle \theta  [degrees]';  ty = 'y'; 
   
   xP = t/(pi);   yP = y;
   plot(xP,yP,'lineWidth',LW);
  
   xlabel('phase angle   \theta  / \pi');
   ylabel('y');
   title('y = cosec(\theta)');
   grid on; box on;
   %set(gca,'XLim',[-360 360]);
   set(gca,'Xtick',-2:0.5:2);
   set(gca,'YLim',[-Ylimit Ylimit]); 
   