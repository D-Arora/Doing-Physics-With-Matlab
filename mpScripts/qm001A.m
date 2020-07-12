% qm001.m



clear
close all
clc


global E k
tic

% INPUTS ==============================================================

% Energy of particles  [eV];
  E = 200;
% X axis  [nm]
  xMin = 0;
  xMax = 5;
% Grid points 
  N = 500;

  % CONSTANTS ===========================================================

  e = 1.6021766208e-19;           % Elementary charge [C]
  h = 6.626070040e-34;            % Planck constant [ J.s]
  hbar = 1.054571800e-34;         % Reduced Planck's constant
  me = 9.10938356e-31;            % Electron mass [kg]

  s1 = 1e-9;                       % m --> nm
  s2 = e;                          % J --> eV
 
   
% SETUP ===============================================================
  
  k = 2*me*s1^2*s2/hbar^2;
  
  k = 1*pi^2;
  
% ode45
  xSpan = linspace(xMin,xMax,N);
  u0 = [1,1i];
  options = odeset('RelTol',1e-6);
 
  [x,u] = ode45(@FNode, xSpan, u0, options);
  
  psi = u(:,1);
  toc
 
  
  % GRAPHICS ==========================================================
  
  figure(1)
  
  plot(x,real(psi))
  
  
 
 % FUNCTIONS ==========================================================
 
 function du = FNode(x,u)
 global E k
 
  y = u(1);
  ydot = u(2);
  du = zeros(2,1);
  du(1) = ydot;
  du(2) = -k*y;
 end

 