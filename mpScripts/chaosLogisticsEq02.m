% chaosLogisticsEq01.m
% Logistic Equation
% Bifurcation Diagram

% DOING PHYSICS WITH MATLAB
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/???
% Scripts
%   https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%   https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% IAN COOPER
% email    matlabvisualphysics@gmail.com

% Matlab 2021b     220526

clear
close all
clc

tic

% Control parameter 0 <= r <= 1
  r = 0.77;

% Initial condition   0 <= x0 <= 1
  x0 = 0.5001;

% Number of iterations
  N = 100;

% Logistic Equation

  x = zeros(N,1);
  x(1) = x0;
  
 for c = 1:N-1
     x(c+1) = 4*r*x(c)*(1-x(c));
 end
     xEND = x(c+1);
     xSTABLE = 1-1/(4*r);

% GRAPHICS  =========================================================== 

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'Position', [0.05 0.05 0.25 0.30])
  set(gcf,'color','w');
  FS = 14;
   hold on
   xP = 1:N; yP = x;
   plot(xP,yP,'b','linewidth',2)
   txt = sprintf('x(0) = %2.3f    r = %2.3f   x_{end} = %2.3f   \n',x0, r, xEND);
   title(txt)
   xlabel('n')
   ylabel('x(n)')
   ylim([0 1.1])
   grid on
   box on
   set(gca,'fontsize',FS)
 
   


 toc
