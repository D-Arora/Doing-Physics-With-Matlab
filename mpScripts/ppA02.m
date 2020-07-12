% ppA01.m


clear 
close all
clc

tic

hold on

N = 50000;
x = zeros(N,1);
y = zeros(N,1);
L = 20;
x(1) = -5;
y(1) = 15;

tMax = 1.5;
t = linspace(0,tMax,N);
dt = t(2)-t(1);


x(2) = x(1) + dt*y(1)*(x(1)+1); 
y(2) = y(1) + dt*x(1)*(2*y(1)+3);

flagS = 0; c = 3;
while flagS == 0
   x(c) = x(c-2) + 2*dt*y(c-1)*(x(c-1)+1); 
   y(c) = y(c-2) + 2*dt*x(c-1)*(y(c-1)+3);
   if abs(x(c)) > L; flagS = 1; end
   if abs(x(c)) > L; flagS = 1; end
   if c > N-10; flagS = 1; end
   c = c+1;
end

XY = 1:c-1;


% Phase space setup: dimensions [-10 to + 10] / number of vectors nX [16] 
   x1Min = -10; 
   x1Max = 10;
   x2Min = -10; 
   x2Max = 10;
   nX = 16;

  x1 = linspace(-L,L,nX);
  x2 = x1;
  [xx, yy] = meshgrid(x1,x2);
  f = yy.*(xx+1);
  g = xx.*(yy+3);
  
  fs = f./sqrt(f.^2 + g.^2);    % unit vectors
  gs = g./sqrt(f.^2  +g.^2);

FS = 14;  
figure(1)
pos = [0.05 0.05 0.25 0.29];
set(gcf,'Units','normalized');
set(gcf,'Position',pos);
set(gcf,'color','w');
plot(t(XY),x(XY),'b')
hold on
plot(t(XY),y(XY),'r')
axis square
grid on
set(gca,'fontsize',FS)

figure(2)
  pos = [0.35 0.05 0.29 0.39];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
   
   hq = quiver(xx,yy,fs,gs);
   xlim([x1Min x1Max])
   ylim([x2Min x2Max])
   set(hq,'color',[0.2 0.2 0.2],'AutoScaleFactor',0.4);
   set(gca,'fontsize',FS)
   
   hold on
   
   plot(x(XY),y(XY),'b','linewidth',2.5)
   xlim([-L L])
   ylim([-L L])
   grid on
   
   
 
     



% figure(2)
% pos = [0.42 0.05 0.29 0.29];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
% plot(x(XY),y(XY))
% xlim([-L L])
% ylim([-L L])
% grid on

toc

