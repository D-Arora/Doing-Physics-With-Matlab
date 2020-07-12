% mec_car_01.m
% 10 jun 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% ../mphome.htm
% SI units used for all quantities
clear all
close all
clc

% INPUTS =================================================================
% mass of car
   m = 1324;
% coeff of friction
   u = 0.70;
% max time
   tMax = 18.1;
% Number of time intervals
   nT = 15101;
% coeff of rolling resistance
  tau = 0.015;
% cross-sectional area
   A = 2.4;
% drag coefficient
  Cd = 0.33;
% air density
   rho = 1.2;
% engine power
   P = 165e3;
% power loss factor
   f = 0.90;
% max time for display
  tLim = tMax-0.1;
%   VW GOLF GTI 2015 acceleration data
   tData = [2.3 3.5 4.6 5.9 7.8 9.7 11.9 15];
   vData = 1.609344 .* [30 40 50 60 70 80 90 100];
   
   
% SETUP ==================================================================
tMin = 0;                    % initial time
g = 9.8;                     % acc due to gravity
FG = (m*g).*ones(1,nT);                    % weight of car

t = linspace(tMin,tMax,nT);

a = zeros(1,nT);            % acceleration
v = zeros(1,nT);            % velocity
x = zeros(1,nT);            % displacement
Pke = zeros(1,nT);
dt = t(2)-t(1);             % time increment
Fnet = zeros(1,nT);         % net force acting on car
Fdrive = u .* FG;       % driving force on car due to friction
Frolling = tau .* FG;   % rolling resistance
Fdrag = zeros(1,nT);                  % air drag
Fnet(1) = Fdrive(1)-Frolling(1);
a(1) = Fnet(1) / m;

% Numerical finite difference procedure ==================================
v(2) = v(1) + a(1) * dt;
%v(2) = 0.5 * (v(1) + v(2));
x(2) = x(1) + v(1)*dt + 0.5*a(1)*dt^2;
Fdrag(2) = 0.5 * Cd * A * rho * v(2)^2;

for c = 3:nT
   Fdrag(c-1) = 0.5 * Cd * A * rho * v(c-1)^2;
   if Fdrive(c-1) > f*P / v(c-1); Fdrive(c-1) = f*P / v(c-1); end;
   Fnet(c-1) = Fdrive(c-1) - Frolling(c-1) - Fdrag(c-1);
   a(c-1)= Fnet(c-1)./m;
   v(c) = v(c-2) + (2*dt) .* a(c-1);
end

for c = 3:nT
   x(c) = x(c-2) + (2*dt) .* v(c-1);
end
   a(nT) = (v(nT)-v(nT-1))/dt;
   Fdrive(nT) = f*P / v(nT);
   Fdrag(nT) = 0.5 * Cd * A * rho * v(nT)^2;
   Fnet(nT) = Fdrive(nT) - Frolling(nT) - Fdrag(nT);

   KE = (0.5*m) .* v.^2;
   
   for c = 2 : nT-1
      Pke(c) = (KE(c+1) - KE(c-1))/(2*dt);
   end
      %Pke(1) = (KE(2) - KE(1))/dt; Pke(nT) = (KE(end) - KE(end-1))/dt;
       Pke(1) = Pke(2); Pke(end) = Pke(end-1);
       
% GRAPHICS ===============================================================
figure(1)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.1 0.1 0.3 0.6]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'time  t  [s]';  ty = 'displacement s [m]'; 
   xP = t;   yP = x;
   subplot(3,1,1);
   plot(xP,yP,'lineWidth',LW);
   xlabel(tx); ylabel(ty);
   grid on; box on;
   set(gca,'Xlim',[0 tLim]); 
   
   col = 'm';
   tx = 'time t [s]';  ty = 'velocity v [km.h^{-1}]'; 
   xP = t;   yP = 3.6.*v;
   subplot(3,1,2);
   plot(xP,yP,'lineWidth',LW);
   xlabel(tx); ylabel(ty);
   grid on; box on;
   set(gca,'Xlim',[0 tLim]);
   hold on 
   xP = tData; yP = vData;
   plot(xP,yP,'o');
     
   col = 'r';
   tx = 'time t [s]';  ty = 'acceleration a [m.s^{-2}]'; 
   xP = t;   yP = a;
   subplot(3,1,3);
   plot(xP,yP,'lineWidth',LW);
   xlabel(tx); ylabel(ty);
   grid on; box on;
   set(gca,'Ylim',[0 10]);
   set(gca,'Xlim',[0 tLim]);  
   
figure(2)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.5 0.5 0.3 0.3]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'time  t  [s]';  ty = 'force  F [N]'; 
   %subplot(2,1,1);
   
   hold on
   
   xP = t;   yP = -Frolling;
   plot(xP,yP,'lineWidth',LW);
   xlabel(tx); ylabel(ty);
     
   xP = t;   yP = -Fdrag;
   plot(xP,yP,'lineWidth',LW);
   
   xP = t;   yP = Fdrive;
   plot(xP,yP,'lineWidth',LW);
   
   xP = t;   yP = Fnet;
   plot(xP,yP,'lineWidth',LW);
   
   xlabel(tx); ylabel(ty);
   %set(gca,'Ylim',[0 1e4]);   
   set(gca,'Xlim',[0 tLim]); 
   legend('Frolling','Fdrag','Fdrive','Fnet');
   grid on; box on;

figure(3)
   fs = 12; LW = 2;
   set(gcf,'units','normalized','position',[0.5 0.1 0.3 0.3]);
   set(gca,'fontsize',fs);
   
   col = 'b';
   tx = 'time  t  [s]';  ty = 'kinetic energy KE [J]'; 
   %subplot(2,1,1);
   
   hold on
   
   xP = t;   yP = KE; yP2 = Pke;
   plotyy(xP,yP,xP,yP2);
   xlabel(tx); ylabel(ty);
%  
hold on
    xP = t;   yP2 = Fnet.*v;
     %plotyy(xP,yP,xP,yP2);
%    
%    xP = t;   yP = Fdrive;
%    plot(xP,yP,'lineWidth',LW);
%    
%    xP = t;   yP = Fnet;
%    plot(xP,yP,'lineWidth',LW);
   
   xlabel(tx); ylabel(ty);
   %set(gca,'Ylim',[0 1e4]);   
   set(gca,'Xlim',[0 tLim]); 
%   legend('Frolling','Fdrag','Fdrive','Fnet');
   grid on; box on;   
  