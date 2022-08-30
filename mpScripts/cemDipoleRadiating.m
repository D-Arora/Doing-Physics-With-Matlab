% cemDipoleRadiating

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/


clear; %close all
clc

% Constant C
  CMin = 0.25; CMax = 1.35; dC = 0.30;
  C = CMin:dC:CMax;
% time
  wtMin =  0; wtMax = 2*pi; dwt = pi/64;
  wt = wtMin:dwt:wtMax;
% Displacement
  krMin = 0.05; krMax = 12; dkr = 0.05;
  kr = krMin:dkr:krMax; fr = krMax/11;

  for c = 1:length(C)
  s2 = abs(C(c)./( (1./kr).*(sin(kr-wt(1))) - cos(kr - wt(1) )  ));

 figure(1)
 set(gcf,'color',[1 1 1])
 for n = 1:length(s2)
% if s2(n) < 1 && s2(n) > 0
 x = real(fr.*kr.*sqrt(s2));
 y = real(fr.*kr.*sqrt(1-s2));
 
 plot(x,y,'b-','linewidth',1)
 hold on
 plot(x, -y,'b-','linewidth',1);
 plot(-x,-y,'b-','linewidth',1);
 plot(-x, y,'b-','linewidth',1);
 plot([-12 12], [0 0],'color',[1 1 1],'linewidth',1.0)
 xlim([-12 12])
 ylim([-12 12])
 axis square
 axis off
  end
  end