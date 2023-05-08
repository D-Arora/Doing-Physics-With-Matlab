% QMG1E.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230415   Matlab R2021b

% GRIFFITH Problem 1.16
clear; close all;clc
% SETUP
  hbar = 1.05457182e-34;
  a = 1;
  x1 = -a; x2 = a; N = 9999;
  x = linspace(x1,x2,N);
  dx = x(2)-x(1);
  psi = (a^2 - x.^2);
% Normalize the wavefuntion
  AN = simpson1d(psi.^2,x1,x2);
  psi = psi./sqrt(AN);
  A = 1/sqrt(AN);
  check = simpson1d(psi.^2,x1,x2);
% Most likely postion to find the particle
  xM = x(psi == max(psi));
% Expectation <x>
   fn = psi.*x.*psi;
   xavg = simpson1d(fn,x1,x2);
% Expectation <x2>
  fn = psi.*x.^2.*psi;
  x2avg = simpson1d(fn,x1,x2);
  sigmaX = sqrt(x2avg - xavg^2);
% Momentum operator
  Gpsi = gradient(psi,dx);
  fn = psi.*Gpsi;
  pavg = -1i*hbar*simpson1d(fn,x1,x2);
  G2psi = gradient(Gpsi,dx);
  fn = psi.*G2psi;
  p2avg = (-1i*hbar)^2*simpson1d(fn,x1,x2);
  sigmaP = sqrt(p2avg - conj(pavg).*pavg);
 
% OUTPUT
   fprintf('Check normalization: prob = %2.4f  \n', check)
   fprintf('Normalized amplitude, A = %2.4f  \n', A)
   fprintf('Most likely postion, xM = %2.4f  \n', xM)
   fprintf('Expectation value, <x> = %2.4f   \n', xavg)
   fprintf('Standard deviation, sigmaX = %2.4f   \n', sigmaX) 
   fprintf('Expectation value <p>  real = %2.4f   imag = %2.4f  \n', real(pavg), imag(pavg))
   fprintf('Standard deviation, sigmaP/hbar = %2.4f   \n', sigmaP/hbar)
   fprintf('Uncertainty Principle: sigmaX*sigmaP/hbar = %2.4f > 0.5  \n', sigmaX*sigmaP/hbar)
  
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
  ylabel('\psi')
  title('Wavefunctiopn','FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(2,1,2)  
  xP = x; yP = psi.^2;
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlabel('x')
  ylabel('\psi^2')
  title('Probability density','FontWeight','normal')
  set(gca,'FontSize',FS)


