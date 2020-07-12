% mec_stokes.m

% 21 apr 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% ../mphome.htm
% Motion of an object falling through a fluid
% Stokes Law - numerical calculation of a, v and x
% S.I. units used for all physical quantities

clear all
close all
clc
hold on
% INPUTS -----------------------------------------------------------------
  R = 0.5;%2e-4;                 % radius of object
  CD = 0.5; %1e3;                 % drag coefficient
  rho = 7800;                      % density of object
  rhoF = 1.225 ;%1.2;               % density of fluid
  visc = 1.789e-5 ;%1.8e-5;         % viscosity of fluid
  tMin = 0;
  tMax = 20;
  N = 10000;
% SETUP ------------------------------------------------------------------
  g = 9.8;                  % acceleration due to gravity
  Vol = (4/3)*pi*R^3;       % volume of object
  m = Vol*rho;              % mass of object
  m =  70; % 2.7e-3;              % mass instead of calc from density
  FG = m*g;                 % weigth of object
  A = pi*R^2;               % Cross-sectional area of object
  v0 = 0;                   % initial velocity of object
 
  k1 = 0.5 * CD * A * rhoF / m;     % Drag force constant
  
  vT = sqrt(2 * m * g / (CD * A * rhoF));    % terminal velocity

  t = linspace(tMin,tMax,N);          % time
  dt = t(2)-t(1);                     % time step
  x = zeros(N,1);                     % displacement
  v = zeros(N,1);                     % velocity
  a = zeros(N,1);                     % acceleration

  % Starting Values  
  x(1) = 0;
  v(1) = v0;
  if v(1) == 0 
     a(1) = g;
  else
        a(1) = g - k1 * v(1)^3 / abs(v(1));
  end
  
  v(2) = v(1) + a(1) * dt;
  a(2) = g - k1 * v(2)^3 / abs(v(2));
  v(2) = v(1) + 0.5 *(a(1)+a(2)) * dt;
  

% Finite Difference procedure
for k = 1 : N-2
    v(k+2) = v(k) + 2*dt * (g - k1 * v(k+1)^3 / abs(v(k+1)));
end

 for k = 1 : N-2
     x(k+2) = x(k) + 2*dt * v(k+1);
 end
 
  for k = 3 : N
     a(k) = g - k1 * v(k).^3 ./ abs(v(k));
  end
  
  NR = (2*rhoF*R/visc) .* abs(v);           % Reynolds number
 
% forces
F_G = FG .* ones(1,N);
F_R = m.*a - F_G';
  
% GRAPHICS ================================================================
%xValues = 0:0.1:0.4;
figure(1)    % velocity / time
col = 'b';
fs = 9;
set(gcf,'units','normalized','position',[0.2 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'velocity  v  [m.s^{-1}]';
xP = t; yP = v;

%Xrange = [0 1.0 * max(xP)];
Xrange = [0 1.0*tMax];
%Yrange = [0 1.1 * max(yP)];
  plot(xP,yP,col,'lineWidth',2)
xlabel(tx); ylabel(ty);
set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
grid on

hold on
 plot([0 tMax],[vT vT],'r');
box on
figure(2)   % acceleration / time  aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
fs = 9;
set(gcf,'units','normalized','position',[0.4 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'acceleration  a  [m.s^{-2}]';
xP = t; yP = a;

% Xrange = [0 1.0 * max(xP)];
% Yrange = [-1.1*max(yP) 1.1 * max(yP)];
% Xrange = [0 1.0 * 2];
% Yrange = [-1.1*10 1.1 * 10];
   plot(xP,yP,col','lineWidth',2)
grid on
xlabel(tx); ylabel(ty);
set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
% 
figure(3)   % displacement / time xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
fs = 9;
set(gcf,'units','normalized','position',[0.6 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'displacement  x  [m]';
xP = t; yP = x;

% Xrange = [0 1.0 * max(xP)];
% Yrange = [0 1.1 * max(yP)];
% Xrange = [0 1.0 * 2];
% Yrange = [-20 1.1 * 40];
   plot(xP,yP,col','lineWidth',2)
grid on
xlabel(tx); ylabel(ty);
set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
hold on
box on

figure(4)   % displacement  /  velocity  xvxvxvxvxvxvxvxvxvxvxvxvxvxvxvxv
fs = 9;
set(gcf,'units','normalized','position',[0.8 0.2 0.18 0.22]);
tx = 'velocity  v   [m.s^{-1}]';
ty = 'displacement  x  [m]';
xP = v; yP = x;

% Xrange = [0 1.0 * max(xP)];
% Yrange = [0 1.1 * max(yP)];
% Xrange = [0 1.0 * 2];
% Yrange = [0 1.1 * 10];
   plot(xP,yP,col','lineWidth',2)
grid on
xlabel(tx); ylabel(ty);
%set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);

figure(5)   % forces / time  ffffffffffffffffffffffffffffffffffffffffffff
fs = 9;
set(gcf,'units','normalized','position',[0.1 0.5 0.18 0.22]);
tx = 'time  t  [s]';
ty = 'Force  [N]';
xP = t; yP = F_G;

% Xrange = [0 1.0 * max(xP)];
% Yrange = [0 1.1 * max(yP)];
% Xrange = [0 1.0 * 2];
% Yrange = [0 1.1 * 10];
   plot(xP,yP,col','lineWidth',2)
grid on
xlabel(tx); ylabel(ty);

hold on
   col = 'r';
   yP = F_R;
   plot(xP,yP,col','lineWidth',2)
   
   yP = F_G + F_R' ;
  
   col = 'k';
   plot(xP,yP,col','lineWidth',1)
   
set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
%legend('F_G','F_R','F_{net}','Orientation','horizontal');
title('    ');

figure(6)   % Reynolds Number / time RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
fs = 9;
col = 'b';
set(gcf,'units','normalized','position',[0.4 0.5 0.18 0.22]);
tx = 'time  t  [s]';
ty = 'Reynolds No.  N_R';

xP = t; yP = NR;

% Xrange = [0 1.0 * max(xP)];
% Yrange = [0 1.1 * max(yP)];
% Xrange = [0 1.0 * 2];
% Yrange = [0 1.1 * 10];
   plot(xP,yP,col','lineWidth',2)
grid on
xlabel(tx); ylabel(ty);

set(gca,'Xlim',Xrange);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);


