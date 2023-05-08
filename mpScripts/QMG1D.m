% QMC1D.m
% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230508   Matlab R2021b

% RUN CELLS #! and CELL #2 separately

%% CELL 1
clear; close all;clc

% SETUP
  A = 1; a = 2; b = 5;
  N = 999;
  x = linspace(0,b,N);
  psi = (A/a).*x;
  psi(x>a) = A.*(b-x(x>a))./(b-a);
% Normalize the wavefuntion
  AN = simpson1d(psi.^2,0,b);
  psi = psi./sqrt(AN);
  A = A/sqrt(AN);
  check = simpson1d(psi.^2,0,b);
% Most likely postion to find the particle
  xM = x(psi == max(psi));
% Probability of x < a
  ind1 = 1;ind2 = find(x>a,1);
  R = ind1:ind2;
  % WARNING R must be an odd number
  L = length(R); if mod(L,2) == 0, ind2 = ind2+1; end
  R = ind1:ind2;    
  probA = simpson1d(psi(R).^2,x(ind1),x(ind2));
% Probability of x > a
  ind1 = find(x>a,1); ind2 = N;
  R = ind1:ind2;
  % WARNING R must be an odd number
  L = length(R); if mod(L,2) == 0, ind2 = ind2+1; end
  R = ind1:ind2;    
  probB = simpson1d(psi(R).^2,x(ind1),x(ind2));
% Expectation <x>
  fn = psi.*x.*psi;
  xavg = simpson1d(fn,0,b);
% Expectation <x2>
  fn = psi.*x.^2.*psi;
  x2avg = simpson1d(fn,0,b);
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
  fprintf('Probability x < a = %2.4f  \n', probA)
  fprintf('Probability x > a = %2.4f  \n', probB)
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

%%  CELL #2
clear; close all; clc

syms x a b A
fn1 = x/a
S1 = int(fn1^2, [0 a])
fn2 = ((b-x)/(b-a))
S2 = int(fn2^2,[a b])
S = S1 + S2
A = 1/sqrt(S)
psi1 = A*fn1
psi2 = A*fn2
probA1 = int(psi1^2,[0 a])
probA2 = int(psi2^2,[a b])

xavg1 = int(x*psi1^2,[0 a])
xavg2 = int(x*psi2^2,[a b])
xavg = xavg1 + xavg2

