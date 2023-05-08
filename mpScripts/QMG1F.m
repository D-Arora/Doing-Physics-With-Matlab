% QMG1E.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG01.pdf
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230415   Matlab R2021b

clear; close all;clc
% SETUP
  n = 1;
  L = 0.5e-9;
  k0 = pi/(L);
  x1 = -L/2; x2 = L/2; N = 9999;
  x = linspace(x1,x2,N);
  dx = x(2)-x(1);
  psi = sin(n*k0*(L/2-x));
  h = 6.62607015e-34;
  hbar = 1.05457182e-34;
  m = 9.1093837e-31;
  e = 1.60217663e-19;
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
  % Energies
  Tavg = (p2avg/(2*m))/e;
  Eavg = Tavg;
  % THEORY
    wL = 2*L/n;
    p = h/wL;
    T = (p^2/(2*m))/e;
    V = 0;
    E = T + V;

% OUTPUT
   fprintf('Quantum number, n = %2.0f  \n',n)
   fprintf('Check normalization: prob = %2.4f  \n', check)
   fprintf('Normalized amplitude, A = %2.4e  \n', A)
   fprintf('Most likely postion, xM = %2.4e  m  \n', xM)
   fprintf('Expectation value, <x> = %2.4e  m  \n', xavg)
   fprintf('Standard deviation, sigmaX = %2.4e  m  \n', sigmaX) 
   fprintf('Expectation value <p>  real = %2.4e  N.s  imag = %2.4e  N.s  \n', real(pavg), imag(pavg))
   fprintf('Expectation value <p*p>  real = %2.4e  N^2.s^2  imag = %2.4e  N^2.s^2  \n', real(p2avg), imag(p2avg))
   fprintf('Uncertainty Principle: sigmaX*sigmaP  = %2.4e  kg.m^2.s^-1  \n', sigmaX*sigmaP)
   fprintf('hbar / 2 = %2.4e  kg.m^2.s^-1 \n',hbar/2')
   fprintf('Wavelength, lambda = %2.4e  m  \n', wL)
   fprintf('Expectation value, <T> = %2.4f  eV  \n', Tavg)
   fprintf('Expectation value, <E> = %2.4f  eV  \n', Eavg)
   fprintf('momentum, p = %2.4e  N.s  \n', p)
   fprintf('kinetic energy, T = %2.4e eV  \n', T)
   fprintf('total energy, E = %2.4e eV \n', E)

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
  ylabel('\psi^2)')
  title('Probability density','FontWeight','normal')
  set(gca,'FontSize',FS)


