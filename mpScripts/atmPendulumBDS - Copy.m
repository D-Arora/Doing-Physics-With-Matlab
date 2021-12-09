% atmPendulumBDS.m
% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% IAN COOPER
%    matlabvisualphysics@gmail.com
% LINK: DOWNLOAD MATLAB SCRIPTS 
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb
% MATLAB VERSIOB R2021b
% 211204

% Bifurcation diagran for a damped driven pendulum
% SI units are for all variables
%  theta - angular displacement / omega - angular frequency / gamma - drive strength
%  ode45 solves the equation of motion

%  Ouput:  Bifurcation diagrams:  theta vs gamma / omega vs gamma
%          
% Supporting documentation
%  https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/atmDDP.htm 


clear
close all
clc

tic


% INPUTS & SETUP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Time 
  tMax = 330;  %600
  N = 9999;
  tMin = 0; 

  tSpan = [tMin tMax];
% tSpan = 0: 0.01: 600;  %linspace(tMin,tMax,N);
 
% Driving signal: strength / period / angular frequency
 
  gammaBD = 1.0:0.0001:1.2;
  
  Lgamma = length(gammaBD);
 
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
  X0 = -pi/2; % -pi/2;
  w0 = 0;
  s0 = [X0 w0];

%  SETUP  ==================================================== 
% Coefficient vector
  K(1) = wN^2;
  K(2) = 2*beta;
  K(4) = wD;
  

% Number of cycle steps;
  nT = 100;

% Solve ODE
  opts = odeset('RelTol',1e-10); 
  tBD = zeros(nT,1);
  xBD = zeros(nT,1); wBD = zeros(nT,1);
  zBD = zeros(nT,1); tBD = zeros(nT,1);


% Solutions X and w stored as a [2D] array

xSol = zeros(Lgamma,nT);
wSol = zeros(Lgamma,nT);
zSol = zeros(Lgamma,nT);
tSol = zeros(Lgamma,nT);

 for k = 1 : Lgamma
     K(3) = gammaBD(k)*wN^2;
     [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 
    

% Angular displacement and Velocity
     w = sol(:,2);
     X = sol(:,1);
   %   w = atan2(tan(w),1);
     t1 = 201+.5;  %501
    % tBD = zeros(100,1); wBD = zeros(100,1); zBD = zeros(100,1);
 
    for c = 1 : nT
      z = find(t >= t1,1);
      xBD(c) = X(z);
      wBD(c) = w(z);
      tBD(c) = t(z);
      zBD(c) = z;
      t1 = t1 + 1;
    end  

    xSol(k,:) = atan2(tan(xBD),1);
    wSol(k,:) = wBD;
    zSol(k,:) = zBD;
    tSol(k,:) = tBD;

end

% Restrict angular displacment -pi to + pi
  % X = atan2(tan(X),1);


% GRAPHICS  ===========================================================
   figure(1)
   pos = [0.05 0.05 0.30 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   hold on

    for k = 1:Lgamma
      HP = plot(gammaBD(k),xSol(k,:),'bo');
      set(HP,'markerfacecolor','b','markersize',1)
      hold on
    end

 %  xlim([1.06 1.087])
 %  ylim([-2 0])

   ylabel('\theta  [ rad ]') 
   xlabel('\gamma')
   grid on
   box on
   xtickformat('%2.3f')
   ytickformat('%2.2f')
   set(gca,'fontsize',12)
   

figure(2)
   pos = [0.55 0.05 0.30 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   hold on

    for k = 1:Lgamma
      HP = plot(gammaBD(k),wSol(k,:),'bo');
      set(HP,'markerfacecolor','b','markersize',1)
      hold on
    end

 %  xlim([1.06 1.087])
 %  ylim([-2 0])

   ylabel('\omega  [ rad/s ]') 
   xlabel('\gamma')
   grid on
   box on
   xtickformat('%2.3f')
   ytickformat('%2.2f')
   set(gca,'fontsize',12)
   

   toc



%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -K(1)*sin(X) - K(2)*w + K(3)*cos(K(4)*t);

  sDot = sDot';


end

