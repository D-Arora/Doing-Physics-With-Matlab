% QMG1C.m
% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm

% IAN COOPER
% matlabvisualphysics@gmail.com
% 230508   Matlab R2021b

clear; close all; clc

% SETUP: Gaussian function - probability density,rho
% Peak value
  A = 1;
% Decay constant
  k = 2;
% Limits x
  L1 = -2; L2 = 4;
% Centre of Gaussian function 
  a = 1;
% x values
  N = 9999;
  x = linspace(L1,L2,N);
% Probability density, rho
  rho = A.*exp(-k.*(x-a).^2);
% Normalize probability density
  AN = simpson1d(rho,L1,L2);
  A = A/AN;
% Gaussian function: probability density, rho
  rho = A.*exp(-k.*(x-a).^2);
  check = simpson1d(rho,L1,L2);
% Expectation value <x>, xavg
  fn = x.*rho;
  xavg = simpson1d(fn,L1,L2);
% Expectation value <x^2>, x2avg
  fn = x.^2.*rho;
  x2avg = simpson1d(fn,L1,L2);
% Standard deviation, sigma
  sigma = sqrt(x2avg - xavg^2);

% Theoretical predictions
  AT = 1/(sigma*sqrt(2*pi));
  kT = 1/(2*sigma^2);
  xavgT = a;

% OUTPUT
  fprintf('Check normalization: prob = %2.4f  \n', check)
  fprintf('normalized amplitude, A = %2.4f  \n', A)
  fprintf('THEORY: normalized amplitude, A = %2.4f  \n', AT)
  fprintf('decay constant, k = %2.4f  \n', k)
  fprintf('THEORY: decay constant, k = %2.4f  \n', kT)
  fprintf('expectation value, <x> = %2.4f  \n', xavg)
  fprintf('THEORY: expectation value, <x> = %2.4f  \n', xavgT)
  fprintf('expectation value, <x^2> = %2.4f  \n', x2avg)
  fprintf('Standard deviation, sigma = %2.4f  \n', sigma)


% GRAPHICS
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.20 0.20]);
  set(gcf,'color','w');
  FS = 14;
  plot(x,rho,'b','LineWidth',2)
  xticks([L1:1:L2])
  grid on
  xlabel('x')
  ylabel('\rho')
  set(gca,'FontSize',FS)

 