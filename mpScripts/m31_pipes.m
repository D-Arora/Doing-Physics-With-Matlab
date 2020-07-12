% m31_pipes.m

close all
clear all
clc



L = 1;
x = linspace(0,L,5500);

n = 4;
m = 7;
wL = 2*L / n;

k = 2*pi/wL;

A = 1;

sP = A .* sin(k*x);
sS = A .* sin(k*x+pi/2);

figure(1)

fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.25]);
set(gcf,'color',[1 1 1]);

subplot(2,1,1)
xP = x; yP = sP;
plot(xP,yP,'b','lineWidth',2);
hold on
plot(xP,-yP,'r','lineWidth',2);
xP = [0 1]; yP = [1.1 1.1];
plot(xP,-yP,'k','lineWidth',3);
plot(xP,yP,'k','lineWidth',3);
set(gca,'fontsize',fs);
axis off
title('displacement');
set(gca,'fontsize',fs);
axis off
title('pressure','FontWeight','normal');

subplot(2,1,2)
xP = x; yP = sS;
plot(xP,yP,'b','lineWidth',2);
hold on
plot(xP,-yP,'r','lineWidth',2);
xP = [0 1]; yP = [1.1 1.1];
plot(xP,-yP,'k','lineWidth',3);
plot(xP,yP,'k','lineWidth',3);
set(gca,'fontsize',fs);
axis off
title('displacement','FontWeight','normal')

figure(99)   % 99999999999999999999999999999999999999999999999999999999999
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.2 0.4 0.2 0.15]);
set(gcf,'color',[1 1 1]);
xP = x; yP = sP;
plot(xP,yP,'b','lineWidth',2);
hold on
plot(xP,-yP,'b','lineWidth',2);
xP = [0 1]; yP = [1.1 1.1];
plot(xP,-yP,'k','lineWidth',3);
plot(xP,yP,'k','lineWidth',3);
set(gca,'fontsize',fs);
axis off
set(gca,'fontsize',fs);
axis off
t1 = 'pressure  n =  ';
t2 = num2str(n,'%2.0f \n');
tm = [t1 t2];
title(tm,'FontWeight','normal');


figure(2)   % 222222222222222222222222222222222222222222222222222222222
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.5 0.2 0.28 0.15]);
set(gcf,'color',[1 1 1]);

xP = x; yP = sP.*sP;
plot(xP,yP,'m','lineWidth',2);
%axis off
%title('average pressure','FontWeight','normal')
set(gca,'fontsize',12);
grid on
set(gca,'xTick',0:0.25:1);
xlabel('x  [ m ]');
ylabel('P_{avg}  [ a.u.]')


figure(3)  % ======================================================
wL = 4*L / m;

k = 2*pi/wL;

sP = A .* sin(k*x);
sS = A .* sin(k*x+pi/2);
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.3 0.4 0.2 0.25]);
set(gcf,'color',[1 1 1]);

subplot(2,1,1)
xP = x; yP = sP;
plot(xP,yP,'b','lineWidth',2);
hold on
plot(xP,-yP,'b','lineWidth',2);
xP = [0 1]; yP = [1.1 1.1];
plot(xP,-yP,'k','lineWidth',3);
plot(xP,yP,'k','lineWidth',3);
xP = [x(end) x(end)]; yP = [-1.15 1.15];
plot(xP,yP,'k','lineWidth',5);
set(gca,'fontsize',fs);
axis off
%title('displacement');
set(gca,'fontsize',fs);
axis off
title('pressure','FontWeight','normal');

subplot(2,1,2)
xP = x; yP = sS;
plot(xP,yP,'b','lineWidth',2);
hold on
plot(xP,-yP,'b','lineWidth',2);
xP = [0 1]; yP = [1.1 1.1];
plot(xP,-yP,'k','lineWidth',3);
plot(xP,yP,'k','lineWidth',3);
xP = [x(end) x(end)]; yP = [-1.15 1.15];
plot(xP,yP,'k','lineWidth',5);
set(gca,'fontsize',fs);
axis off
title('displacement','FontWeight','normal')


figure(88)   % 99999999999999999999999999999999999999999999999999999999999
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.6 0.4 0.2 0.15]);
set(gcf,'color',[1 1 1]);
xP = x; yP = sP;
plot(xP,yP,'b','lineWidth',2);
hold on
plot(xP,-yP,'b','lineWidth',2);
xP = [0 1]; yP = [1.1 1.1];
plot(xP,-yP,'k','lineWidth',3);
plot(xP,yP,'k','lineWidth',3);
xP = [x(end) x(end)]; yP = [-1.15 1.15];
plot(xP,yP,'k','lineWidth',5);
set(gca,'fontsize',fs);
axis off
set(gca,'fontsize',fs);
axis off
t1 = 'pressure  m =  ';
t2 = num2str(m,'%2.0f \n');
tm = [t1 t2];
title(tm,'FontWeight','normal');
