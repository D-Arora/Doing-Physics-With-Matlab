% QMG1G.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230415   Matlab R2021b

% GRIFFITH  problem 1.10
clear; close all;clc
% SETUP
 PI = '3141592653589793238462643';
 L = length(PI);
 N = zeros(10,1);
 PIs = sort(PI);
% Frequency of each digit 
 for c = 1:10
   N(c) = count(PI,num2str(c-1));
 end
% Probabilities 
 prob = N./25;
% Expectation values
  c = 1:9;
  r = N(2:end)';
  cavg = sum(c.*r)/L;
  c2avg = sum(c.^2.*r)/L;
  sigma = sqrt(c2avg - cavg^2);

% OUTPUT 
  disp(PI)
  disp(PIs)
  fprintf('Median value = %2.0f  \n', str2double(PIs(12)))
  fprintf('Most common value = %2.0f  \n', c(N == max(N))-1)
  fprintf('Average value = %2.4f  \n', cavg)
  fprintf('standard deviation = %2.4f  \n', sigma)
  ct = 0:9; ct = ct';
  table(ct,N,prob)

% GRAPHICS
 figure(1)
   figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.2 0.2]);
  set(gcf,'color','w');
  FS = 14;
   x = 0:9;
   bar(x,N,'b')
   xlabel('c')
   ylabel('N(c)')
   set(gca,'fontsize',14)