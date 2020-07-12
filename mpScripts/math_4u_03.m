% math_4u_03.m


%%
clear all
close all
clc

xMin = 0.001;
xMax= 5;
yMin = -5;
yMax = 5;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = x;
y2 = log(x);
y = y1 .* y2;
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   yP = y; col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
   
   xlabel(tx); ylabel(ty);
   
   grid on; box on;
  % axis equal
   set(gca,'yTick',yMin:1:yMax);
   set(gca,'xLim',[xMin xMax]);
   set(gca,'yLim',[yMin yMax]);
   legend('x','log_ex','x log_ex');
   
   col = 'r';
   plot(2,2,'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
   col = 'm';
   plot(2,log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
   col = 'b';
   plot(2,2*log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
   
       
  %%   y = x + 2sin(3x)
clear all
close all
clc

xMin = -pi;
xMax= pi;
yMin = -5;
yMax = 5;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = x;
y2 = 2.*sin(3*x);
y = y1 + y2;
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x  [rad]';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   yP = y; col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
   
   xlabel(tx); ylabel(ty);
   
   grid on; box on;
  % axis equal
 ax = gca;
  set(gca,'yTick',yMin:1:yMax);
  %set(gca,'xLim',[xMini xMax]);
 % set(gca,'xTick',[xMin/pi :1/6:xMax/pi]);
 
 ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
 legend('x','2 sin(x)','x + 2 sin(x)');
   
%    col = 'r';
%    plot(2,2,'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    col = 'm';
%    plot(2,log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    col = 'b';
%    plot(2,2*log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    
              
  %%   y = | 2sin(3x-3) + x % ===============================================================================================================
clear all
close all
clc

xMin = -pi;
xMax= pi;
yMin = -5;
yMax = 5;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = x;
y2 = 2.*sin(3.*x-3);
y = y1 + y2;
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x  [rad]';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   yP = abs(y); col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
   
   xlabel(tx); ylabel(ty);
   
   grid on; box on;
  % axis equal
 ax = gca;
  set(gca,'yTick',yMin:1:yMax);
  %set(gca,'xLim',[xMini xMax]);
 % set(gca,'xTick',[xMin/pi :1/6:xMax/pi]);
 
 %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
 legend('x','2 sin(x)','|x + 2 sin(x)|');
   
%    col = 'r';
%    plot(2,2,'o','MarkerFaceColor',col,'MarkerEdgeColor',col);1+pi/2

%    col = 'm';
%    plot(2,log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    col = 'b';
%    plot(2,2*log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    



    %%   y = x exp(-x) % ===============================================================================================================
clear all
close all
clc

xMin = -1;
xMax= 5;
yMin = -20;
yMax = 30;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = 3.*x;
y2 = exp(3-x);
y = y1 .* y2;
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   yP = y; col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
   
   xP = [xMin xMax]; yP = [0 0];
   plot(xP,yP,'k', 'lineWidth',0.5);
   
   yP = [yMin yMax]; xP = [0 0];
   plot(xP,yP,'k', 'lineWidth',0.5);
   
   xlabel(tx); ylabel(ty);
   
   
   grid on; box on;
  % axis equal
 ax = gca;
 % set(gca,'yTick',yMin:1:yMax);
 set(gca,'yLim',[yMin yMax]);
 set(gca,'YMinorGrid','on');
 set(gca,'MinorGridAlpha', 0.5);
 set(gca,'MinorGridLineStyle','-');
 set(gca,'MinorGridColor',[0.8 0.8 0.8]);
 set(gca,'GridColor','k');
 %set(gca,'YMinorGrid','on');
 % set(gca,'xTick',[xMin/pi :1/6:xMax/pi]);
 
 %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
 legend('3x','e^{3-x}','3x * e^{3-x}');
   
%    col = 'r';
%    plot(2,2,'o','MarkerFaceColor',col,'MarkerEdgeColor',col);1+pi/2

%    col = 'm';
%    plot(2,log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    col = 'b';
%    plot(2,2*log(2),'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
%    



    %%   y = x*(x+1)/(x-2)) % ===============================================================================================================
clear all
close all
clc

xMin = 0;
xMax= 20;
yMin = -50;
yMax = 50;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = x.^2 + 10.*x;
y2 = x-2;
y = sqrt(y1 ./ y2);
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   yP = y; col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
   
   xP = [xMin xMax]; yP = [0 0];
   plot(xP,yP,'k', 'lineWidth',0.5);
   
   yP = [yMin yMax]; xP = [0 0];
   plot(xP,yP,'k', 'lineWidth',0.5);
   
   xlabel(tx); ylabel(ty);
   
   
   grid on; box on;
  % axis equal
 ax = gca;
 % set(gca,'yTick',yMin:1:yMax);
 set(gca,'yLim',[yMin yMax]);
 set(gca,'YMinorGrid','on');
 set(gca,'MinorGridAlpha', 0.5);
 set(gca,'MinorGridLineStyle','-');
 set(gca,'MinorGridColor',[0.8 0.8 0.8]);
 set(gca,'GridColor','k');
 %set(gca,'YMinorGrid','on');
% set(gca,'xTick',xMin:1:xMax);
 
 %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
 legend('x^2 + 10 x','x - 2','(x^2 + 10 x) / (x - 2)');
   
    %%   y = x*(x+1)/(x-2)) % ===============================================================================================================
clear all
close all
clc

xMin = 0;
xMax= 20;
yMin = -50;
yMax = 50;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = x.^2 + 10.*x;
y2 = x-2;
y = sqrt(y1 ./ y2);
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   yP = y; col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
   
   xP = [xMin xMax]; yP = [0 0];
   plot(xP,yP,'k', 'lineWidth',0.5);
   
   yP = [yMin yMax]; xP = [0 0];
   plot(xP,yP,'k', 'lineWidth',0.5);
   
   xlabel(tx); ylabel(ty);
   
   
   grid on; box on;
  % axis equal
 ax = gca;
 % set(gca,'yTick',yMin:1:yMax);
 set(gca,'yLim',[yMin yMax]);
 set(gca,'YMinorGrid','on');
 set(gca,'MinorGridAlpha', 0.5);
 set(gca,'MinorGridLineStyle','-');
 set(gca,'MinorGridColor',[0.8 0.8 0.8]);
 set(gca,'GridColor','k');
 %set(gca,'YMinorGrid','on');
% set(gca,'xTick',xMin:1:xMax);
 
 %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
 legend('x^2 + 10 x','x - 2','(x^2 + 10 x) / (x - 2)');
   
    %%   |x-3| + |x-5| = 10 % ===============================================================================================================
clear all
close all
clc

xMin = -5;
xMax= 15;
yMin = -5;
yMax = 20;
N = 5500;

x = linspace(xMin,xMax,N);
y1 = abs(x-3);
y2 = abs(x-5);
y12 = y1+y2;

figure(1)
   fs = 12; LW = 1;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1; col = 'r';
   plot(xP,yP,col,'lineWidth',LW);
   hold on
   
   yP = y2; col = 'm';
   plot(xP,yP,col,'lineWidth',LW);
     
   LW = 2;
   yP = y12; col = 'b';
   plot(xP,yP,col, 'lineWidth',LW);
%   
    LW = 1;
    xP = [xMin xMax]; yP = [10 10];
    plot(xP,yP,'k', 'lineWidth',LW);
%    
%    yP = [yMin yMax]; xP = [0 0];
%    plot(xP,yP,'k', 'lineWidth',0.5);
%    
    xlabel(tx); ylabel(ty);
%    
   
   grid on; box on;
  % axis equal
% ax = gca;
 % set(gca,'yTick',yMin:1:yMax);
  set(gca,'yLim',[yMin yMax]);
%  set(gca,'YMinorGrid','on');
%  set(gca,'MinorGridAlpha', 0.5);
%  set(gca,'MinorGridLineStyle','-');
%  set(gca,'MinorGridColor',[0.8 0.8 0.8]);
%  set(gca,'GridColor','k');
 %set(gca,'YMinorGrid','on');
set(gca,'xTick',xMin:5:xMax);
 
 %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
legend('|x-3|','|x - 5|','|x-3|+|x-5|', 'y = 10');
   

    %%   ellipse ===============================================================================================================
clear all
close all
clc

a = 5;
b = 3;
tMin = 0;
tMax= 2*pi;

xMax = 10; xMin = - xMax;
yMax = 10; yMin = -yMax;

N = 5500;

t = linspace(tMin,tMax,N);
x = a .* cos(t);
y = b.* sin(t);

xP = 4;
yP = sqrt(b^2 - (b/a)^2 .* xP^2);

x1 = linspace(xMin,xMax,N);
M1 = -(b/a)^2 * (xP/yP);
B1 = yP + (b/a)^2 * (xP^2/yP);

y1 = M1 .* x1 + B1;

M2 = -1/M1;
B2 = yP - M2 * xP;

y2 = M2.* x1 + B2;

figure(1)
   fs = 12; LW = 1;
   set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   X = x;   Y = y; col = 'b';
   plot(X,Y,col,'lineWidth',LW);
   hold on
   
   X = x1; Y = y1; col = 'r';
   plot(X,Y,col,'lineWidth',LW);
   
   X = x1; Y = y2; col = 'm';
   plot(X,Y,col,'lineWidth',LW);
%      
%    LW = 2;
%    yP = y12; col = 'b';
%    plot(xP,yP,col, 'lineWidth',LW);
% %   
%     LW = 1;
%     xP = [xMin xMax]; yP = [10 10];
%     plot(xP,yP,'k', 'lineWidth',LW);
%    
%    yP = [yMin yMax]; xP = [0 0];
%    plot(xP,yP,'k', 'lineWidth',0.5);
%    
  xlabel(tx); ylabel(ty);
%    
   
   grid on; box on;
  % axis equal
% ax = gca;
  set(gca,'xLim',[xMin xMax]);
  set(gca,'yLim',[yMin yMax]);
  axis equal
%  set(gca,'YMinorGrid','on');
%  set(gca,'MinorGridAlpha', 0.5);
%  set(gca,'MinorGridLineStyle','-');
%  set(gca,'MinorGridColor',[0.8 0.8 0.8]);
%  set(gca,'GridColor','k');
 %set(gca,'YMinorGrid','on');
%set(gca,'xTick',xMin:5:xMax);
 
 %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
 %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
 %ax.xLim = [-pi, pi];
 %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
%legend('|x-3|','|x - 5|','|x-3|+|x-5|', 'y = 10');
   
set(gca,'xLim',[xMin xMax]);
set(gca,'yLim',[yMin yMax]);

%%
clear all
close all
clc

p = 3/10; q = 7/10; n = 10;


for c = 1:10
  k = c-1;
  x(c) = k;
  P(c) = nchoosek(n,k) .* p.^k .* q.^(n-k);
end;

figure(1)
bar(x,P)
xlabel('# black balls' );
ylabel('Prob');

%%
clear all
close all
clc


n1 = 8
k1 = 3
n2 = 15
k2= 0

c1 = nchoosek(n1,k1)
c2 = nchoosek(n2,k2)

c1
c2

c = c1*c2






      