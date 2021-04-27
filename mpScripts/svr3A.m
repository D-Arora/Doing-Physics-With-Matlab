% svr3.m

% Numerical Model for the time evolution of a disease spreading through a
% population

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/sird1.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b
% 200304

close all
clear
clc

% *********************************************************************  
% INPUTS: Model Parameters
% ********************************************************************* 

% Model: number of time steps 
  num = 5000;
% Initial number of infections   [0.001] 
  I = zeros(num,1);
  I(1) = 0.001;
% Maximum number of days for model   [100]
  tMax = 500;
% X axis tick marks 
  dXT = 100;
  XT = 0:dXT:tMax;
% Adjustible model parameters  
% S --> I  [0.35] 
  a = 0.6;
% I --> R  [0.03]
  b = 0.2;
% R --> S   [0] 
  c = 0.004;
% I --> D  [0.01]
  d = 0.01;

% Oscillations: tMax = 500 a = 0.6 b = 0.2 c = 0.004  d = 0.01   

% *********************************************************************  
% SETUP
% ********************************************************************* 
  
% Initialize arrays
  N = zeros(num,1);
  S = zeros(num,1);
  D = zeros(num,1);
  R = zeros(num,1);
   Itot = zeros(num,1);
 
  
% Initial values
  N(1) = 100;
  S(1) = N(1) - I(1); 
  D(1) = 0;
  R(1) = 0;

  t = linspace(0,tMax,num);
  dt = t(2) - t(1);


% *********************************************************************  
% MODEL: Time Evolution
% ********************************************************************* 
   
for n = 2:num
  S(n) = S(n-1) - dt * (a*S(n-1)*I(n-1)/N(n-1) - c*R(n-1));
  I(n) = I(n-1) + dt * (a*S(n-1)*I(n-1)/N(n-1) - b*I(n-1));
  R(n) = R(n-1) + dt * (b*I(n-1) - c*R(n-1));
  Itot(n) = Itot(n-1) + dt * a*S(n-1)*I(n-1)/N(n-1);
  
  if S(n) < 0; S(n) = 0; end
  
  deaths = dt * d*I(n-1);
  D(n) = D(n-1) + deaths;
  I(n) = I(n) - deaths;
  N(n) = S(n) + I(n) + R(n);
end


% *********************************************************************  
% GRAPHICS
% ********************************************************************* 

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.4 0.6]);
  set(gcf,'color','w');

% Y axis tick marks
  YT = 0:20:100;
subplot(3,2,1)  
  xP = t; yP = N;
    plot(xP,yP,'k','linewidth',2)
    grid on
    set(gca,'xtick',XT)
    xlim([0, tMax])
   % xlabel('time elapsed  [a.u.]')
    ylabel('N %');
    set(gca,'fontsize',12)
    
     txt = sprintf('I(0) = %3.4f   a = %2.3f   b = %2.3f',I(1),a,b);
     Htitle = title(txt);
     set(Htitle,'FontWeight','normal')    
  
subplot(3,2,2)  
  xP = t; yP = S;
    plot(xP,yP,'k','linewidth',2)
    grid on
    set(gca,'xtick',XT)
    set(gca,'ytick',YT)
    xlim([0, tMax])
 %   xlabel('time elapsed  [a.u.]')
    ylabel('S  %');
    set(gca,'fontsize',12)
    
    txt = sprintf('c = %2.3f   d = %2.3f',c, d);
    Htitle = title(txt);
    set(Htitle,'FontWeight','normal')   
    
subplot(3,2,3)  
  xP = t; yP = I;
    plot(xP,yP,'m','linewidth',2)
    grid on
    set(gca,'xtick',XT)
    xlim([0, tMax])
 %   xlabel('time elapsed  [a.u.]')
    ylabel('I  %');
    set(gca,'fontsize',12)   
    
subplot(3,2,4)  
  xP = t; yP = Itot;
    plot(xP,yP,'m','linewidth',2)
    grid on
    set(gca,'xtick',XT)
    set(gca,'ytick',0:25:max(Itot))
    xlim([0, tMax])
%    xlabel('time elapsed  [a.u.]')
    ylabel('I_{tot}');
    set(gca,'fontsize',12)       
    
 subplot(3,2,5)  
  xP = t; yP = R;
    plot(xP,yP,'b','linewidth',2)
    grid on
    set(gca,'xtick',XT)
    set(gca,'ytick',YT)
    xlim([0, tMax])
    xlabel('time elapsed  [a.u.]')
    ylabel('R  %');
    set(gca,'fontsize',12)   
    
subplot(3,2,6)  
  xP = t; yP = D;
    plot(xP,yP,'r','linewidth',2)
    grid on
    set(gca,'xtick',XT)
    xlim([0, tMax])
    xlabel('time elapsed  [a.u.]')
    ylabel('D  %');
    set(gca,'fontsize',12)       
    
