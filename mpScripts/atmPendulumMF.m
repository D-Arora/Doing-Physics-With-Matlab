% atmPendulum.m
%211119


clear
close all
clc


tic
% SI units

% INPUTS & SETUP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Time 
  tMax = 10;
  N = 9999;
  tMin = 0; 

  tSpan = [tMin tMax];
  %tSpan = linspace(tMin,tMax,N);
 
% Driving signal: strength / period / angular frequency
  gamma = 1.01;
  TD = 1;
  wD = 2*pi/TD;

% Natutal frequency
  wN = 1.5*wD;
  TN = 2*pi/wN;

% Damping   
  beta = wN/4;
%  beta = 0.5;

% Initial conditions: displacement / velocity
% X = theta   x = omega  
  X0 = 0; -pi/2;
  w0 = 0;
  s0 = [X0 w0];

%  SETUP  ==================================================== 
% Coefficient vector
  K(1) = wN^2;
  K(2) = 2*beta;
  K(3) = gamma*wN^2;
  K(4) = wD;

% Solve ODE
  opts = odeset('RelTol',1e-6);  
  [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 

% Angular displacement and angular velocity
  X = sol(:,1); w = sol(:,2);
% Restrict angular displacment -pi to + pi
  % X = atan2(tan(X),1);
% Driving signal and frequencies
  F = K(3).*cos(K(4)*t);
  fN = 1/TN;
  fD = 1/TD;

 [pks, loc] = findpeaks(X,t);

  TP = (loc(end) - loc(end-2))/2
 
% OUPUTS ===========================================================
   flagP = 1;  

   txt = sprintf('theta(0) = %2.3f  \n',X0);
    disp(txt) 
  txt = sprintf('omega(0) = %2.3f  \n',w0);
    disp(txt) 
  txt = sprintf('Driving period TD = %2.3f  \n',TD);
    disp(txt);
  txt = sprintf('Driving freq wD = %2.3f  \n',wD);
    disp(txt)
  txt = sprintf('Natural period TN = %2.3f  \n',TN);
    disp(txt)
  txt = sprintf('Natural freq  wN = %2.3f  \n',wN);
    disp(txt)
  if flagP == 1
  txt = sprintf('period from peaks  TP = %2.3f  \n',TP);
    disp(txt)
  end
  txt = sprintf('beta = %2.3f  \n',beta);
    disp(txt) 
  txt = sprintf('gamma = %2.3f  \n',gamma);
    disp(txt)   


% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.25 0.50];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
 subplot(3,1,1) 
  yyaxis left
    xP = t;
   yP = X./pi;
    yP = X;
    plot(xP,yP,'b','linewidth',2)
    grid on
   ylabel('\theta  [ rad / \pi ]')
    ylabel('\theta  [ rad ]')
    xlabel('t  [ s ]')
   hold on
%     txt = sprintf('L = %2.2f   b = %2.2f   A_D = %2.2f   \\omega_D = %2.2f \n',L,b,AD,wD);
%     title(txt,'FontWeight','normal')
    set(gca,'fontsize',12)
    hold on
    set(gca,'YColor','b')

  yyaxis right
    xP = t; yP = F;   
    plot(xP,yP,'r-','linewidth',1)
    grid on
    ylabel('F_D  [ rad.s^{-2} ]')
    set(gca,'YColor','r')


subplot(3,1,2)   
    xP = t; yP = w;   
    plot(xP,yP,'b','linewidth',2)
    grid on
    ylabel('\omega  [ rad.s^{-1} ]')
    xlabel('t  [ s ]')
    set(gca,'fontsize',12)
    hold on

    
%  subplot(3,1,3)   
%    xP = L.*sin(X); yP = -L.*cos(X); 
%     xP = X./pi; yP = w; 
%     plot(xP,yP,'b','linewidth',2)
% 
%     hold on
%     HP = plot(xP(1),yP(1),'go');
%     set(HP,'markerfacecolor','g','markersize',10)
% 
%   xlim([-30 1])
%     grid on
%     xlabel('\theta  [ rad/\pi ]')
%     ylabel('\omega  [ rad.s^{-1} ]')
%     set(gca,'fontsize',12)


 figure(9)
   pos = [0.5 0.1 0.25 0.250];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
  
   xP = t; yP = X;   
   plot(xP,yP,'b-','linewidth',2)
   ylabel('\theta  [ rad ]')
   xlabel('t  [ s ]')
   grid on
    txt = sprintf('\\gamma = %2.4f  \n',gamma);
   title(txt)
   set(gca,'fontsize',14)
   hold on

figure(10)
   pos = [0.5 0.5 0.28 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
 
subplot(1,2,1)
   xP = t; yP = X./pi;   
   plot(xP,yP,'b-','linewidth',2)
   ylabel('\theta  [ rad / \pi ]')
   xlabel('t  [ s ]')
   grid on
   xlim([0 30])
   txt = sprintf('\\gamma = %2.4f  \n',gamma);
   title(txt)
   set(gca,'fontsize',14)
   hold on

   subplot(1,2,2)
   xP = t; yP = X./pi;   
   plot(xP,yP,'b-','linewidth',2)
   ylabel('\theta  [ rad / \pi ]')
   xlabel('t  [ s ]')
   xlim([20 30])
   xticks(20:2:30)
  ylim([0.98 1.0])
   grid on
   txt = sprintf('T = %2.2f s \n',1.0);
  title(txt)
   set(gca,'fontsize',14)
   hold on


   toc

%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -K(1)*sin(X) - K(2)*w + K(3)*cos(K(4)*t);

  sDot = sDot';


end

