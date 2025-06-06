% atmPendulumPS.m

% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% IAN COOPER
%    matlabvisualphysics@gmail.com
% LINK: DOWNLOAD MATLAB SCRIPTS 
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb
% MATLAB VERSIOB R2021b
% 211204

% POINCARE SECTION for a damped driven pendulum
% SI units are for all variables
%  theta - angular displacement / omega - angular frequency / gamma - drive strength
%  ode45 solves the equation of motion

%  Ouput:  Phase space plots
%          
% Supporting documentation
%  https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/atmDDP.htm 



clear
close all
clc


tic
% SI units

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Time 
  tMax = 60000;
  N = 99999;
  tMin = 0; 

% POINCARE SECTION - starting time
  tStart = 10;

% Drive strength
  gamma = 1.5;

% Initial conditions: displacement X = theta / velocity w = omega
  X0 = -pi/2;
  w0 = 0;  


% SETUP  ===========================================================
  tSpan = [tMin tMax];
%  tSpan = linspace(tMin,tMax,N);
 
% Driving signal: strength / period / angular frequency
  TD = 1;
  wD = 2*pi/TD;

% Natutal frequency
  wN = 1.5*wD;
  TN = 2*pi/wN;

% Damping   
  beta = wN/8;

% Initial conditions vector
  s0 = [X0 w0];

% Coefficient vector
  K(1) = wN^2;
  K(2) = 2*beta;
  K(3) = gamma*wN^2;
  K(4) = wD;

% Solve ODE
  opts = odeset('RelTol',1e-10,'AbsTol',1e-6);  
  [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 

% Angular displacement and angular velocity
  X = sol(:,1); w = sol(:,2);
% Restrict angular displacment -pi to + pi
  XR = atan2(tan(X),1);

% Driving signal and frequencies
%   F = K(3).*cos(K(4)*t);
%   fN = 1/TN;
%   fD = 1/TD;


% POINCARE SECTION ==================================================
  tPS = tStart:1:t(end)-1;
  tPSLen = length(tPS);
  n = zeros(tPSLen,1);
  for c = 1 : tPSLen
      n(c) = find(t >= tPS(c),1);
  end
 
   XPS = X(n);
   wPS = w(n);
  
 %  XPSR = XPS; %atan2(tan(XPS),1); 
   XPSR = atan2(tan(XPS),1)./pi; 

% OUPUTS ===========================================================
  
%    txt = sprintf('theta(0) = %2.3f  \n',X0);
%     disp(txt) 
%   txt = sprintf('omega(0) = %2.3f  \n',w0);
%     disp(txt) 
%   txt = sprintf('Driving period TD = %2.3f  \n',TD);
%     disp(txt);
%   txt = sprintf('Driving freq wD = %2.3f  \n',wD);
%     disp(txt)
%   txt = sprintf('Natural period TN = %2.3f  \n',TN);
%     disp(txt)
%   txt = sprintf('Natural freq  wN = %2.3f  \n',wN);
%     disp(txt)
%   txt = sprintf('beta = %2.3f  \n',beta);
%     disp(txt) 
%   txt = sprintf('gamma = %2.3f  \n',gamma);
%     disp(txt)   


% GRAPHICS  ===========================================================

 figure(9)
   pos = [0.05 0.5 0.2 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

subplot(2,1,1)
   xP = t; yP = X;
   plot(xP,yP,'b','linewidth',2)
   grid on
   hold on

   xP = tPS; yP = XPS;
   HP = plot(xP,yP,'ro');
   set(HP,'markerfacecolor','r','markersize',2)

   ylabel('\theta  [ rad / \pi ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',14)

 subplot(2,1,2)
   xP = t; yP = w;
   plot(xP,yP,'b','linewidth',2)
   grid on
   hold on

   xP = tPS; yP = wPS;
   HP = plot(xP,yP,'ro');
   set(HP,'markerfacecolor','r','markersize',2)

   ylabel('\omega  [ rad / s ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',14)  



figure(1)
   pos = [0.05 0.05 0.2 0.2];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

   Xlen = length(X);
   R = Xlen-100:Xlen;

    xP = X./pi;
    yP = w;
    plot(xP,yP,'b','linewidth',2)
    grid on
    hold on
    HP = plot(xP(1),yP(1),'go');
    set(HP,'markerfacecolor','g','markersize',8)
   
    plot(xP(R),yP(R),'r','linewidth',2)

    HP = plot(XPS,wPS,'mo');
    set(HP,'markerfacecolor','m','markersize',6)


    ylabel('\omega  [ rad / s ]')
    xlabel('\theta  [ rad / \pi ]')
    txt = sprintf('\\gamma =  %2.5f \n',gamma);
    title(txt,'FontWeight','normal')
    set(gca,'fontsize',12)
    limX = [min(xticks), max(xticks)];
    limY = [min(yticks), max(yticks)];



figure(2)
   pos = [0.5 0.05 0.3 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
       
    HP = plot(XPSR,wPS,'bo');
    set(HP,'markerfacecolor','b','markersize',1)
    hold on
  %  HP = plot(xP,yP,'r-','linewidth',0.5);

    grid on
    xlim([-1.1 1.1])  
%    xlim([limX(1) limX(2)])   
    ylim([limY(1) limY(2)])    
    ylabel('\omega  [ rad / s ]')
    xlabel('\theta  [ rad / \pi ]')
    set(gca,'fontsize',14)





   toc

%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -K(1)*sin(X) - K(2)*w + K(3)*cos(K(4)*t);

  sDot = sDot';


end

