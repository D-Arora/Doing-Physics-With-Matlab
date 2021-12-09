% atmPendulum.m
%211119


clear
close all
clc

% SI units

% SETUP  ==============================================================

% Time 
  tMin = 0; tMax = 20001;  N = 99999;
% acceleration due to gravity
  g = 9.8;
% Length of pendulum
  L = 9.8;
% damping coefficient
  b = 0.5;%0.2;
% Driving amplitude, period and frequency
  AD = 0.8;
  
  wD = 2/3;  %TD = 5;  wD = 2*pi/TD; %wD = 2/3;

% Initial angular  X == theta  [rad]
  X0 = 1;
% Initial angular velocity  w == omega  
  w0 = 3;
% Coefficient vector
  K = [g L b AD wD];

% Simulation time
  tSpan = linspace(tMin, tMax, N);
%  tSpan = [tMin tMax];
% Initial conditions  
  s0 = [X0 w0];


opts = odeset('RelTol',1e-6);  

[t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0,opts); 


% Angular displacement and angular velocity
  X = sol(:,1); w = sol(:,2);
% Driving input
  F = AD.*cos(wD*t);

 % XR = atan2(tan(X),1);



% Natural frequency and period
  wN = sqrt(g/L)
  TN = 2*pi/wN
  TD = 2*pi/wD


  [pks, loc] = findpeaks(X,t);

  TP = (loc(end) - loc(end-2))/2


% % POINCARE SECTION ==================================================
  dt = t(2) - t(1);
  NP = round(TD/dt);
  P = 1:NP:N;
  tP = t(P);
  XP = X(P);
%  XP = atan2(tan(XP),1);
  wP = w(P);

% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.30 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
%  subplot(2,1,1) 
%  % yyaxis left
%     xP = t; yP = X;   
%     plot(xP,yP,'b','linewidth',2)
%     grid on
%    % xlim([0 140])
%    % xticks(0:20:140)
%     ylabel('\theta  [ rad / \pi ]')
%     xlabel('t  [ s ]')
% 
%     txt = sprintf('L = %2.2f   b = %2.2f   A_D = %2.2f   \\omega_D = %2.2f \n',L,b,AD,wD);
%     title(txt,'FontWeight','normal')
%     set(gca,'fontsize',12)
%     hold on
 %   set(gca,'YColor','b')

%   yyaxis right
%     xP = t; yP = F;   
%     plot(xP,yP,'r-','linewidth',1)
%     grid on
%     ylabel('F_D  [ rad.s^{-2} ]')
%     set(gca,'YColor','r')
% 

 %subplot(2,1,2)   
     xP = XP./pi; yP = wP;   
     HP = plot(xP,yP,'bo');
     set(HP,'markerfacecolor','b','markersize',2)
     grid on
     ylabel('\omega  [ rad.s^{-1} ]')
     xlabel('\theta  [ rad / \pi ]')
%     set(gca,'fontsize',12)
% 
% 
% % POINCARE SECTION 
% 
%  subplot(3,1,3)   
%    % xP = L.*sin(X); yP = -L.*cos(X); 
%     xP = XP./pi; yP = w(P); 
%     plot(xP,yP,'b','linewidth',2)
% 
%     hold on
%     HP = plot(xP(1),yP(1),'go');
%     set(HP,'markerfacecolor','g','markersize',10)
% 
%  %  xlim([-30 1])
%     grid on
%     xlabel('\theta  [ rad/\pi ]')
%     ylabel('\omega  [ rad.s^{-1} ]')
%     set(gca,'fontsize',12)

%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  g = K(1); L = K(2); b = K(3); AD = K(4); wD = K(5);
  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -(g/L)*sin(X) - b*w + AD*cos(wD*t);

  sDot = sDot';


end

