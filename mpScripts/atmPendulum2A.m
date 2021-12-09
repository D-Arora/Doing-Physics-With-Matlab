% atmPendulum.m
%
% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% IAN COOPER
%    matlabvisualphysics@gmail.com
% LINK: DOWNLOAD MATLAB SCRIPTS 
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb
% MATLAB VERSIOB R2021b
% 211204

% Damped Driven pendulum:  Lyapunov Exponents
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

% SI units

% INPUTS & SETUP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Time 
  tMax = 20;
  N = 99999;
  tMin = 0; 

 % tSpan = [tMin tMax];
  tSpan = linspace(tMin,tMax,N);
 
% Driving signal: strength / period / angular frequency
  gamma = 1.075;
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
  X0 = 0;
  w0 = 0;
  s0 = [X0 w0];

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
  
 % SOLUTION 2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  X01 = 0.001;
  s0 = [X01,w0];
 [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 
 X1  = sol(:,1); w1 = sol(:,2);
 t1 = t;
 
 % Natural frequency and period
  wN = 1.5*wD;
  TN = 2*pi/wN
  TD = 2*pi/wD


 


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

   txt = sprintf('\\gamma = %2.4f  \n',gamma);
   title(txt)
%  title(txt,'FontWeight','normal')
    set(gca,'fontsize',12)
   
   

  

subplot(2,1,2)   
 %   X = atan2(tan(X),1);
 %   X1 = atan2(tan(X1),1);
    dX = abs(X - X1);
 
    xP = t; yP = dX;   
    %plot(xP,yP,'b','linewidth',2)
    semilogy(xP,yP,'b','linewidth',2)
    grid on
    yticks(10.^(-10:2:0))
    ylabel('| \Delta\theta |  [ rad ]')
    xlabel('t  [ s ]')
    set(gca,'fontsize',12)



%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -K(1)*sin(X) - K(2)*w + K(3)*cos(K(4)*t);

  sDot = sDot';

end


