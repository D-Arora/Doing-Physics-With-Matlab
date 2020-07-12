%m5_simulations.m
clear all
close all
clc


v1 = 3.24;
theta = 60;
ay = -9.81;

tMin = 0; tMax = 1;
Nt = 2000;

t = linspace(tMin,tMax, Nt);

v1x = v1 * cosd(theta);
v1y = v1 * sind(theta);

x = v1x .* t;
y = v1y .* t + 0.5*ay.*t.^2;
vx = v1x .* ones(1,length(t));
vy = v1y + ay .* t;

N = find(y<0,1);


% calculation dt
s2 = 0.1281 - 0.06711
n2 = 1
s3 = 0.2695 - 0.06711
n3 = 17
s_dt =sqrt( 2*(s3-s2*n3/n2) / (ay*(n3^2-n2*n3)) )

sf = 1/11.65;

d2 = sf*0.74 
d3 = sf*2.39 
d_dt =sqrt( 2*(s3-s2*n3/n2) / (ay*(n3^2-n2*n3)) )

% graphical data
xP = [0 0.57 0.97 1.54 2.02 2.50 2.94 3.41 3.90 4.36 4.88 5.34 5.80 6.27 6.73 7.21 7.59 8.09];
yP = [0 0.74 1.42 1.95 2.47 2.80 3.26 3.49 3.65 3.79 3.85 3.85 3.75 3.66 3.41 3.11 2.73 2.39];
xP = xP.*sf+0.040525;
yP = yP.*sf+0.06711;
ttP = t(2:50:N);
tP = ttP(1:18);

figure(1)
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.1 0.1 0.33 0.3]);
plot(x(1:N),y(1:N),'linewidth',2);
hold on
plot(xP,yP,'o');

set(gca,'fontSize',fs);
xlabel('s_x   [ m ]');
ylabel('s_y   [ m ]');
grid on
%set(gca,'yLim',[0, 1.1*max(y)]);
axis equal
axis([0 1 0 0.5])
%axis off

figure(2)
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.1 0.5 0.33 0.3]);
plot(x(1:50:N),y(1:50:N),'bo');
set(gca,'fontSize',fs);
xlabel('s_x   [ m ]');
ylabel('s_y   [ m ]');
grid on
%set(gca,'yLim',[0, 1.1*max(y)]);
axis equal
axis([0 1 0 0.5])
%axis off

figure(12)
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.1 0.5 0.33 0.3]);
plot(tP,yP,'bo');
hold on
plot(tP,xP,'ro');
set(gca,'fontSize',fs);
xlabel('t   [ m ]');
ylabel('s_y   [ m ]');
grid on
%set(gca,'yLim',[0, 1.1*max(y)]);
%axis equal
axis([0 0.5 0 1])
%axis off


figure(4)
fs = 14;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.5 0.1 0.33 0.8]);

subplot(3,1,1)
plot(t(1:N),y(1:N),'b','linewidth',2);
hold on
plot(t(1:N),x(1:N),'r','linewidth',2);
set(gca,'fontSize',fs);
xlabel('t   [ s ]');
ylabel('s_x   s_y   [ m ]');
grid on
legend('Y','X');

subplot(3,1,2)
plot(t(1:N),vy(1:N),'b','linewidth',2);
hold on
plot(t(1:N),vx(1:N),'r','linewidth',2);
set(gca,'fontSize',fs);
xlabel('t   [ s ]');
ylabel('v_x   v_y   [ m/s ]');
grid on

subplot(3,1,3)
plot(t(1:N),ay.*ones(N,1),'b','linewidth',2);
hold on
plot(t(1:N),zeros(N,1),'r','linewidth',2);
set(gca,'fontSize',fs);
xlabel('t   [ s ]');
ylabel('a_x   a_y   [ m/s^2 ]');
grid on
%set(gca,'yLim',[0, 1.1*max(y)]);
%axis equal
%axis([0 1 0 0.5])
%axis off
