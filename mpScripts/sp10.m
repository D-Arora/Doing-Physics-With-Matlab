% sp10A.m

clear all
close all
clc

g = 9.81;
d = 1;
dt = 0.02;
a = g*sind(45);
ax = a*cosd(45);
ay = -a*sind(45);

tMax = 2.0;

N = 600;

sx = zeros(N,1);
sy = zeros(N,1);
vx = zeros(N,1);
vy = zeros(N,1);
t = zeros(N,1);
sy(1) = d;
flagA = 0;
c = 1;

while flagA < d
  
  t(c+1) = t(c) + dt;
  
  vx(c+1) = vx(1) + ax*t(c+1);
  sx(c+1) = sx(1) + vx(1)*t(c+1) + 0.5*ax*t(c+1)^2;
  
  vy(c+1) = vy(1) + ay*t(c+1);
  sy(c+1) = sy(1) + vy(1)*t(c+1) + 0.5*ay*t(c+1)^2;
  
  flagA = sx(c);
  c = c+1;
end
 G = 6;
c = c-1; 
flagA = 0;
vxx = sqrt(2*g*d);
vyy = 0;
sxx = d;
syy = 0;

  
while flagA < 3*d
  t(c+1) = t(c)+dt;  
  vx(c+1) = vxx;
  vy(c+1) = vyy;
  sx(c+1) = sx(c) + vxx*dt;
  sy(c+1) = syy;
  c = c + 1;
  flagA = sx(c);
end

vxx = sqrt(2*g*d)*cosd(45);
vyy = sqrt(2*g*d)*sind(45);
vy(c) = vyy;
sxx = sx(c);
syy = 0;
txy = t(c-1);
ax = -a*cosd(45);
ay = -a*sind(45);

while flagA > 0
  t(c+1) = t(c) + dt;  
  vx(c+1) = vxx +  ax*(t(c)-txy);
  sx(c+1) = sxx + vxx*(t(c)-txy)  + 0.5*ax*(t(c)-txy)^2;
  vy(c+1)   = vyy +  ay*(t(c)-txy);
  sy(c+1)   = syy +  vyy*(t(c)-txy) + 0.5*ay*(t(c)-txy)^2;
  flagA = vx(c+1);
  c = c + 1;
end

K = find(sx == 0, 2);
K = K(end)-1;

clear ax ay;
ax = zeros(K,1);
ay = zeros(K,1);
ax(1) = a*cosd(45);
ay(1) = -a*sind(45);
ax(end) = -a*cosd(45);
ay(end) = -a*sind(45);

for c = 2:K-1
   ax(c) = (vx(c+1) - vx(c-1))/(2*dt);
   ay(c) = (vy(c+1) - vy(c-1))/(2*dt);
end

% ========================================================================
%  GRAPHICS
% ========================================================================

figure(1)
set(gcf,'units','normalized','position',[0.01 0.2 0.23 0.62]);
subplot(3,1,1)
plot(t(1:K),ax(1:K),'b','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('a_x  [m.s^{-2}]');
set(gca,'fontsize',14);
% set(gca,'xLim',[0 tMax]);
% set(gca,'yTick',-10:5:10);

subplot(3,1,2)
plot(t(1:K),vx(1:K),'r','linewidth',2)
grid on
xlabel('t  [s]');
ylabel('v_x  [m.s^{-1}]');
set(gca,'fontsize',14);
% set(gca,'xLim',[0 tMax]);
% set(gca,'yLim',[0 5]);
% set(gca,'yTick',0:1:4);

subplot(3,1,3);
plot(t(1:K),sx(1:K),'k','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('s_x  [m]');
set(gca,'fontsize',14);
% set(gca,'xLim',[0 tMax]);
% set(gca,'yTick',0:1:4);


figure(2)
set(gcf,'units','normalized','position',[0.41 0.2 0.23 0.62]);
subplot(3,1,1)
plot(t(1:K),ay(1:K),'b','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('a_y  [m.s^{-2}]');
set(gca,'fontsize',14);
% set(gca,'xLim',[0 1.5]);
% set(gca,'yTick',-10:5:10);

subplot(3,1,2)
plot(t(1:K),vy(1:K),'r','linewidth',2)
grid on
xlabel('t  [s]');
ylabel('v_y  [m.s^{-1}]');
set(gca,'fontsize',14);
% set(gca,'xLim',[0 1.5]);
%set(gca,'yTick',0:1:4);

subplot(3,1,3);
plot(t(1:K),sy(1:K),'k','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('s_y  [m]');
set(gca,'fontsize',14);
% set(gca,'xLim',[0 1.5]);
% set(gca,'yTick',0:1:4);