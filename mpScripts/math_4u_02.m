% math_4u_02.m


%%
clear all
close all
clc

xMin = -15;
xMax= 15;
N = 500;

x = linspace(xMin,xMax,N);
y = -2.*x + 10;
y1 = 3.*x - 5;
y2 = 0.5.*x + 20;
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.3 0.5]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
  % xP = x1;   yP = y2;
  % h_area = area(xP,yP,'lineWidth',LW);
  % set(h_area,'FaceColor',[0.5 0.5 0.8],'EdgeColor',[0.5 0.5 0.8]);
  % hold on
   
   xP = x;   yP = y;
   plot(xP,yP,'lineWidth',LW);
   hold on
   col = 'r';
   yP = y1;
   plot(xP,yP,'lineWidth',LW);
   col = 'k';
   yP = y2;
   plot(xP,yP,'lineWidth',LW);
   xlabel(tx); ylabel(ty);
   grid on; box on;
   axis equal
   set(gca,'Xtick',[-20:5:20]);
   set(gca,'yLim',[-10 30]);
   
%%   plotting polynomial functions
clear all
close all
clc

a0 = 48;   
a1 = 32;   
a2 = -19;   
a3 = -2;   
a4 = 1;

xMin = -5;
xMax= 5;
yMax = 100;
yMin = -yMax;
N = 500;

x = linspace(xMin,xMax,N);

y = a0 + a1.*x + a2.*x.^2 + a3.*x.^3 +a4.*x.^4;

%y1 = a0 - a1.*x + a2.*x.^2 - a3.*x.^3 +a4.*x.^4; 
y1 = -y;
figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.3 0.4]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y;
   plot(xP,yP,'lineWidth',LW); 
   hold on
   xP = x;   yP = y1;
   plot(xP,yP,'r','lineWidth',LW); 
   
   
   set(gca,'Xtick',[-5:1:5]);
   set(gca,'yLim',[yMin yMax]);
   xlabel('x','fontsize',fs);
   ylabel('y','fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
   legend('y = f(x)','y = -f(x)')
   
  % xP = [0 0];   yP = [yMin yMax];
    yP = [0 0];   xP = [xMin xMax];
   plot(xP,yP,'k','lineWidth',1); 
   grid on;
   box on;
   
 
%%   plotting rectangular hyperbola   graphs_03.docx
clear all
close all
clc

k1 = 1;
k2 = 10;

xMin = -2.1;
xMax= 2.1;
yMax = 60;
yMin = -yMax;
N = 889;

x = linspace(xMin,xMax,N);

y1 = k1 ./ x;
y2 = k2 ./ x;


figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1;
   plot(xP,yP,'lineWidth',LW); 
   hold on
   xP = x;   yP = y2;
   plot(xP,yP,'r','lineWidth',LW); 
   
   set(gca,'xtick',[-2:0.5:2]);
   set(gca,'ytick',[-yMax:20:yMax]);
   set(gca,'xLim',[(xMin+0.1) (xMax-0.1)]);
   set(gca,'yLim',[yMin yMax]);
   xlabel('x','fontsize',fs);
   ylabel('y','fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
   legend('y = 1 / x','y = 10 / x')
   
  % xP = [0 0];   yP = [yMin yMax];
    yP = [0 0];   xP = [xMin xMax];
   plot(xP,yP,'k','lineWidth',1); 
   grid on;
   box on;
   
  
%%   plotting expontial functions   graphs_03.docx
clear all
close all
clc

k1 = exp(1);   % 0.2
k2 = 10;   % 0.5

xMin = -2.1;
xMax= 2.1;
yMax = 25;
yMin = -5;
N = 889;

x = linspace(xMin,xMax,N);

y1 = k1 .^x;
y2 = k2 .^x;


figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1;
   plot(xP,yP,'lineWidth',LW); 
   hold on
   xP = x;   yP = y2;
   plot(xP,yP,'r','lineWidth',LW); 
   
   set(gca,'xtick',[-2:0.5:2]);
   set(gca,'ytick',[-yMax:5:yMax]);
   set(gca,'xLim',[(xMin+0.1) (xMax-0.1)]);
   set(gca,'yLim',[yMin yMax]);
   xlabel('x','fontsize',fs);
   ylabel('y','fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
   legend('y = e^x','y = 10^x')
   
  % xP = [0 0];   yP = [yMin yMax];
    yP = [0 0];   xP = [xMin xMax];
   plot(xP,yP,'k','lineWidth',1); 
   grid on;
   box on;

   
%%   plotting logarithmic   graphs_03.docx
clear all
close all
clc

k1 = 1;   % 0.2
k2 = 1;   % 0.5

xMin = 0.01;
xMax = 100;
yMax = 5;
yMin = 0;
N = 889;

x = linspace(xMin,xMax,N);

y1 = log10(x);
y2 = log(x);


figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1;
   plot(xP,yP,'lineWidth',LW); 
   hold on
   xP = x;   yP = y2;
   plot(xP,yP,'r','lineWidth',LW); 
   
   set(gca,'xtick',[0:1:10]);
   %set(gca,'ytick',[-yMax:5:yMax]);
   set(gca,'yLim',[-2 5]);
   set(gca,'xLim',[0 11]);
   xlabel('x','fontsize',fs);
   ylabel('y','fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
   legend('y = log_{10}x','y = log_{e}x')
   
  % xP = [0 0];   yP = [yMin yMax];
    yP = [0 0];   xP = [xMin xMax];
   plot(xP,yP,'k','lineWidth',1); 
   grid on;
   box on;


   %%   plotting logarithmic   graphs_03.docx
clear all
close all
clc

k1 = 1;   % 0.2
k2 = 1;   % 0.5

xMin = 0;
xMax = 100;
yMax = 5;
yMin = 0;
N = 889;

x = linspace(xMin,xMax,N);

y1 = x;
y2 = x.^(1/2);
y3 = x.^(1/3);


figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y1;
   plot(xP,yP,'lineWidth',LW); 
   hold on
   xP = x;   yP = y2;
   plot(xP,yP,'r','lineWidth',LW); 
   
   %set(gca,'xtick',[0:1:10]);
   %set(gca,'ytick',[-yMax:5:yMax]);
   set(gca,'yLim',[0 20]);
   %set(gca,'xLim',[0 11]);
   xlabel('x','fontsize',fs);
   ylabel('y','fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
  
   
    xP = x;   yP = y3;
    plot(xP,yP,'m','lineWidth',2); 
    legend('y = x','y = x^{1/2}','y = x^{1/3}')
    grid on;
    box on;
    
    %%   plotting radioactive decay   graphs_03.docx
clear all
close all
clc

t12 = 20;
N0 = 100;
xMin = 0;
xMax = 100;
yMax = 100;
yMin = 0;
N = 1889;

x = linspace(xMin,xMax,N);

y = N0 .* exp(-x.*log(2)./t12);



figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.35]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'x';  ty = 'y'; 
   
   xP = x;   yP = y;
   plot(xP,yP,'lineWidth',LW); 
   hold on
  
   
   %set(gca,'xtick',[0:1:10]);
   set(gca,'ytick',[0:25/2:yMax]);
  % set(gca,'yLim',[0 20]);
   %set(gca,'xLim',[0 11]);
   ylabel('No. unstable nuclei  N','fontsize',fs);
   xlabel('time  t','fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
  
 
    grid on;
    box on;
       
    
  %%   plotting turning pints
clear all
close all
clc


xMin = -2;
xMax = 2;
yMax = 50;
yMin = -50;
N = 1889;

x = linspace(xMin,xMax,N);

%y = x.^3 - 5.*x.^2 + 1.*x + 10; 
%y = (x+1).*(x-2).*(x-4);
y = x.^3 - 5.*x.^2 + 2.*x + 8;
y1 = 3.*x.^2 - 10.*x + 2;
y2 = 6.*x - 10;

figure(1)
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.65]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   
   tx = 'x';  ty = 'y'; 
  
   subplot(3,1,1)
   xP = x;   yP = y;
   plot(xP,yP,'lineWidth',LW); 
   grid on
   set(gca,'xtick',-2:1:6);
   ylabel(ty,'fontsize',fs);
   xlabel(tx,'fontsize',fs);
   
   subplot(3,1,2)
   tx = 'x';  ty = 'dy/dx'; 
   xP = x;   yP = y1;
   plot(xP,yP,'lineWidth',LW)
   grid on
   set(gca,'xtick',-2:1:6);
   ylabel(ty,'fontsize',fs);
   xlabel(tx,'fontsize',fs);
   
   subplot(3,1,3)
   xP = x;   yP = y2;
   tx = 'x';  ty = 'd^2y/dx^2'; 
   plot(xP,yP,'lineWidth',LW)
   grid on
   set(gca,'xtick',-2:1:6);
   %set(gca,'ytick',[0:25/2:yMax]);
  % set(gca,'yLim',[0 20]);
   %set(gca,'xLim',[0 11]);
   ylabel(ty,'fontsize',fs);
   xlabel(tx,'fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
  
 
    grid on;
    box on;
       
    
    %%   plotting curve sketching    ap/maths/graph07.docx
clear all
close all
clc


xMin = -3;
xMax = 3;
yMax = 100;
yMin = -100;
N = 88889;

x = linspace(xMin,xMax,N);

y1 = x.^2 - 3.*x + 2;
y2 = x.^2 + 3.*x + 2;
y = y1./y2;

xx  = -1.9

yy  = (xx.^2 - 3.*xx + 2) / (xx.^2 + 3.*xx + 2)


figure(1)  % ------------------------------------------------------------
   fs = 14; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.25 0.65]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   
   tx = 'x';  ty = 'numerator'; 
  
   subplot(3,1,1)
   xP = x;   yP = y1;
   plot(xP,yP,'lineWidth',LW); 
   grid on
   %set(gca,'xtick',-2:1:6);
   ylabel(ty,'fontsize',fs);
   xlabel(tx,'fontsize',fs);
   
   subplot(3,1,2)
   tx = 'x';  ty = 'demoninator'; 
   xP = x;   yP = y2;
   plot(xP,yP,'lineWidth',LW)
   grid on
   %set(gca,'xtick',-2:1:6);
   ylabel(ty,'fontsize',fs);
   xlabel(tx,'fontsize',fs);
   
   subplot(3,1,3)
   xP = x;   yP = y;
   tx = 'x';  ty = 'y'; 
   plot(xP,yP,'lineWidth',LW)
   grid on
   %set(gca,'xtick',-3:1:6);
   %set(gca,'ytick',[0:25/2:yMax]);
   set(gca,'yLim',[yMin yMax]);
   %set(gca,'xLim',[0 11]);
   ylabel(ty,'fontsize',fs);
   xlabel(tx,'fontsize',fs);
   %legend('y = f(x)','y = f(-x)');
  
 
    grid on;
    box on;
       
         
