% QMG1E.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230508   Matlab R2021b

clear; close all;clc

% SETUP
  A = 1;
  k = 1.5;
  x1 = -2.5; x2 = 2.5; N = 9999;
  x = linspace(x1,x2,N);
  psi = A.*exp(-k*abs(x));
  
% Normalize the wavefuntion
  AN = simpson1d(psi.^2,x1,x2);
  psi = psi./sqrt(AN);
  A = A/sqrt(AN);
  check = simpson1d(psi.^2,x1,x2);
% Most likely postion to find the particle
  xM = x(psi == max(psi));
% % Probability of x < a
%   probA = simpson1d(psi.^2,0,a);  
% % Probability of x > a
%   probB = simpson1d(psi.^2,a,b);
% Expectation <x>
  fn = psi.*x.*psi;
  xavg = simpson1d(fn,x1,x2);
% Expectation <x2>
  fn = psi.*x.^2.*psi;
  x2avg = simpson1d(fn,x1,x2);
  sigma = sqrt(x2avg - xavg^2);
% Probability of xavg-sigm < x < xavg+sigma
  ind1 = find(x>xavg-sigma,1);ind2 = find(x>xavg+sigma,1);
  R = ind1:ind2;
  % WARNING R must be an odd number
  L = length(R); if mod(L,2) == 0, ind2 = ind2+1; end
  R = ind1:ind2;    
  probS = simpson1d(psi(R).^2,x(ind1),x(ind2));

% OUTPUT
  fprintf('Check normalization: prob = %2.4f  \n', check)
  fprintf('Normalized amplitude, A = %2.4f  \n', A)
  fprintf('Most likely postion, xM = %2.4f  \n', xM)
  fprintf('Expectation value, <x> = %2.4f  \n', xavg)
  fprintf('Standard deviation, sigma = %2.4f  \n', sigma) 
  fprintf('Probability x-sigma < x < x+sigma = %2.4f  \n', probS)

% GRAPHICS
  figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.25 0.35]);
  set(gcf,'color','w');
  FS = 14;
subplot(2,1,1)  
  xP = x; yP = psi;
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlabel('x')
  ylabel('\psi(x,0)')
  set(gca,'FontSize',FS)
subplot(2,1,2)  
  xP = x; yP = psi.^2;
  plot(xP,yP,'b','LineWidth',2)
  hold on
  area(x(R),psi(R).*psi(R))
  grid on
  xlabel('x')
  ylabel('|\psi(x,0)|^2')
  set(gca,'FontSize',FS)


