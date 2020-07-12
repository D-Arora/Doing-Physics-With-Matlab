% sp10A.m

clear all
close all
clc

g = 9.81;
d = 0.5*g*sind(45)*0.5^2;

t = 0;
dt = 0.005;
sx = 0; sy = d;
vx = 0; vy =0;

ax(1) =  g*sind(45);
ay(1) = -g*sind(45);

flagA = 0;
flagB = 0;
flagC = 0;
c = 1;
while vx >= 0
  vx(c+1) = vx(c) + ax(c)*dt;
  sx(c+1) = sx(c) + vx(c)*dt + 0.5*ax(c)*dt^2;
  
  vy(c+1) = vy(c) + ay(c)*dt;
  sy(c+1) = sy(c) + vy(c)*dt + 0.5*ay(c)*dt^2;
  t(c+1) = t(c) + dt;
  c = c + 1;
  
  ax(c) = g*sind(45);
  ay(c) = -g*sind(45);
  
  if sx(c) > d
    ax(c) = 0;
    ay(c) = 0;
     if flagA == 0
       vx(c) = sqrt(vx(c)^2 + vy(c)^2);
       vy(c)= 0;
     end
  end

  if sx(c) > 3*d
      flagA = 1;
      ax(c) = -g*sind(45);
      ay(c) = -g*sind(45);
      if flagC == 0
         vy(c) = sqrt(2*vx(c-1));
         vx(c) = sqrt(2*vx(c-1));
         flagC = 1; 
      end
  end    
end

figure(1)
set(gcf,'units','normalized','position',[0.01 0.2 0.23 0.62]);
subplot(3,1,1)
plot(t,ax,'b','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('a_x  [m.s^{-2}]');
set(gca,'fontsize',14);
set(gca,'xLim',[0 1.5]);
set(gca,'yTick',-10:5:10);

subplot(3,1,2)
plot(t,vx,'r','linewidth',2)
grid on
xlabel('t  [s]');
ylabel('v_x  [m.s^{-1}]');
set(gca,'fontsize',14);
set(gca,'xLim',[0 1.5]);
set(gca,'yLim',[0 5]);
set(gca,'yTick',0:1:4);

subplot(3,1,3);
plot(t,sx,'k','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('s_x  [m]');
set(gca,'fontsize',14);
set(gca,'xLim',[0 1.5]);
set(gca,'yTick',0:1:4);


figure(2)
set(gcf,'units','normalized','position',[0.41 0.2 0.23 0.62]);
subplot(3,1,1)
plot(t,ay,'b','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('a_y  [m.s^{-2}]');
set(gca,'fontsize',14);
set(gca,'xLim',[0 1.5]);
% set(gca,'yTick',-10:5:10);

subplot(3,1,2)
plot(t,vy,'r','linewidth',2)
grid on
xlabel('t  [s]');
ylabel('v_y  [m.s^{-1}]');
set(gca,'fontsize',14);
set(gca,'xLim',[0 1.5]);
%set(gca,'yTick',0:1:4);

subplot(3,1,3);
plot(t,sy,'k','linewidth',2);
grid on
xlabel('t  [s]');
ylabel('s_y  [m]');
set(gca,'fontsize',14);
set(gca,'xLim',[0 1.5]);
% set(gca,'yTick',0:1:4);