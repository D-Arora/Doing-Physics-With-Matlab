% atmPendulum.m
%211119


clear
close all
clc

% SI units

% SETUP  ==============================================================

% Time 
  tMin = 0; tMax = 30;
% acceleration due to gravity
  g = 9.8;
% Length of pendulum
  L = 0.25*9.8;
% damping coefficient
  b = 0.2;
% Driving amplitude, period and frequency
  AD = 2.5;%
  wD = 1.26;   % TD = 5; wD = 2*pi/TD;
% Initial angular  X == theta  [rad]
  X0 = 0.2; X01 = 0.201;
% Initial angular velocity  w == omega  
  w0 = 0;
% Coefficient vector
  K = [g L b AD wD];

% Simulation time
%  tSpan = [tMin tMax];
  tSpan = linspace(tMin,tMax,9999);
% Initial conditions  
  s0 = [X0 w0];
  

opts = odeset('RelTol',1e-6);  

[t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 


% Angular displacement and angular velocity
  X = sol(:,1); w = sol(:,2);
% Driving input
  F = AD.*cos(wD*t);

 % X = atan2(tan(X),1);

 % SOLUTION 2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 s0 = [X01,w0];
 [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 
 X1  = sol(:,1); w1 = sol(:,2);

 
 % Natural frequency and period
  wN = sqrt(g/L)
  TN = 2*pi/wN
  TD = 2*pi/wD


  [pks, loc] = findpeaks(X,t);

  TP = (loc(end) - loc(end-2))/2



% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.25 0.50];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
 subplot(2,1,1) 
 
    xP = t; yP = X./pi;   
    plot(xP,yP,'b','linewidth',3)
    hold on
    yP = X1./pi;   
    plot(xP,yP,'r','linewidth',1)
    grid on
    ylabel('\theta  [ rad / \pi ]')
    xlabel('t  [ s ]')

    txt1 = sprintf('\\theta_1(0) = %2.3f \n',X0);
    txt2 = sprintf('\\theta_2(0) = %2.3f  \n',X01);
    legend(txt1, txt2,'location','south','Orientation','horizontal')

    txt = sprintf('L = %2.2f   b = %2.2f   A_D = %2.2f   \\omega_D = %2.2f \n',L,b,AD,wD);
    title(txt,'FontWeight','normal')
    set(gca,'fontsize',12)
   
   

  

subplot(2,1,2)   
    X = atan2(tan(X),1);
    X1 = atan2(tan(X1),1);
    dX = abs(X - X1);
 
    xP = t; yP = dX;   
    %plot(xP,yP,'b','linewidth',2)
    semilogy(xP,yP,'b','linewidth',2)
    grid on
    yticks(10.^(-5:0))
    ylabel('\Delta\theta  [ rad ]')
    xlabel('t  [ s ]')
    set(gca,'fontsize',12)



%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  g = K(1); L = K(2); b = K(3); AD = K(4); wD = K(5);
  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -(g/L)*sin(X) - b*w + AD*cos(wD*t);

  sDot = sDot';


end

