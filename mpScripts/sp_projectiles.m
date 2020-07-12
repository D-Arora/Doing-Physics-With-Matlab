% sp_projectiles.m

clear all
close all
clc


% INPUTS ================================================================
   g = 9.81;     
   u = [8 10]; 
   A = [30 60];
   h = 4;
   
   ax = 0; ay = - g;
   
% Initial velocites
   ux = u .* cosd(A)
   uy = u .* sind(A)
   
% When both balls are their highest positions
%  times, velocities and displacements 
   tH = -uy ./ ay
   sHy = uy.*tH + 0.5 * ay * tH.^2
   sHx = ux .* tH
   
% When ball A is at its highest position: time tH(1)   
%    Position of ball B 
   sAx = sHx(1);
   sAy = sHy(1);
   sBx = ux(2) * tH(1)
   sBy = uy(2) * tH(1) + 0.5*ay*tH(1)^2
   
   sBAx = sBx - sAx
   sBAy = sBy - sAy
   
   sBA = sqrt(sBAx^2 + sBAy^2)
   angleBA = atan2d(sBAy,sBAx)
   
%  When ball A is at its highest position: time tH(1)   
%    velocity of ball B
   vAx = ux(1)
   vAy = 0
   vBx = ux(2)
   vBy = uy(2) + ay * tH(1)
   
   vBAx = vBx - vAx
   vBAy = vBy - vAy
   
   vBA = sqrt(vBAx^2 + vBAy^2)
   anglevBA = atan2d(vBAy,vBAx)
   
% As the balls enter the water
  sy = -h
  
  vAWx = ux(1)
  vAWy = -sqrt(uy(1)^2+2*ay*sy)
  vAW = sqrt(vAWy^2+vAWy^2)
  angleAW = atan2d(vAWy,vAWx)
  
  vBWx = ux(2)
  vBWy = -sqrt(uy(2)^2+2*ay*sy)
  vBW = sqrt(vBWy^2+vBWy^2)
  angleBW = atan2d(vBWy,vBWx)
  
  tAW = (vAWy - uy(1))/ay
  tBW = (vBWy - uy(2))/ay
  
  sAWx = ux(1) * tAW
  sBWx = ux(2) * tBW
  
  
  
%% GRAPHICS --------------------------------------------------------------

close all

N = 2000;
tMin = 0;
tMax = 2.2;
t = linspace(tMin,tMax,N)';

uPx = ux .* ones(N,2);
uPy = uy .* ones(N,2);

vPx = uPx;
vPy = uPy + ay .* t;

sPx = uPx .* t;
sPy = uPy .*t + (0.5*ay).*t.^2;


figure(1)
   fs = 14;
   set(gcf,'units','normalized','position',[0.1 0.1 0.5 0.8]);
    stopA = 1:find(sPy(:,1)<-h,1);
    stopB = 1:find(sPy(:,2)<-h,1);

% Plot trajectory   
   pos1 = [0.1 0.78 0.8 0.2];
   subplot('Position',pos1);
   xP = sPx; yP = sPy;
   plot(xP,yP,'lineWidth',2);
   legend('A','B');
   grid on
   xlabel('s_x   [ m ] ');
   ylabel('s_y   [ m ] ');
   set(gca,'YLim',[-4 4]);
   set(gca,'fontsize',fs);
 
 % PLot sx / t  
   pos1 = [0.1 0.5 0.35 0.15];
   subplot('Position',pos1);
   yP = sPx(stopA,1); xP = t(stopA);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   yP = sPx(stopB,2); xP = t(stopB);
   plot(xP,yP,'r','lineWidth',2);
   legend('A','B','Location','northwest');
   grid on
   xlabel('t   [ s ] ');
   ylabel('s_x   [ m ] ');
   % set(gca,'YLim',[-4 4]);
   set(gca,'fontsize',fs);  

 % Plot sy / t
   pos1 = [0.56 0.5 0.35 0.15];
   subplot('Position',pos1);
   yP = sPy(:,1); xP = t;
   plot(xP,yP,'b','lineWidth',2);
   hold on
   yP = sPy(:,2); 
   plot(xP,yP,'r','lineWidth',2);
   legend('A','B','Location','northeast');
   grid on
   xlabel('t   [ s ] ');
   ylabel('s_y   [ m ] ');
   set(gca,'YLim',[-4 4]);
   set(gca,'fontsize',fs);    
   
 % Plot vx / t
   pos1 = [0.1 0.25 0.35 0.15];
   subplot('Position',pos1);
  
   yP = vPx(stopA,1); xP = t(stopA);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   
   yP = vPx(stopB,2); xP = t(stopB); 
   plot(xP,yP,'r','lineWidth',2);
   legend('A','B','Location','east');
   grid on
   xlabel('t   [ s ] ');
   ylabel('v_x   [ m/s ] ');
   set(gca,'fontsize',fs);    
   
% Plot vy / t
   pos1 = [0.56 0.25 0.35 0.15];
   subplot('Position',pos1);
   yP = vPy(stopA,1); xP = t(stopA);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   yP = vPy(stopB,2); xP  = t(stopB);
   plot(xP,yP,'r','lineWidth',2);
   legend('A','B','Location','east');
   grid on
   xlabel('t   [ s ] ');
   ylabel('v_y   [ m/s ] ');
   set(gca,'fontsize',fs);       