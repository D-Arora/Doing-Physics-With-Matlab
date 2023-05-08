% QMG1H.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230415   Matlab R2021b

% GRIFFITH  problems 1.11 and 1.12
clear; close all;clc
% SETUP
  A = 1;                    % amplitude
  k = 2;                    % spring constant
  m = 2;                    % mass
  E = 0.5*k*A^2;            % total energy
  N = 99999;                % grid points
  x1 = -A; x2 = A;          % limits a and b 
  x = linspace(x1,x2,N);    % x position
  V = 0.5*k*x.^2;           % potential energy
  v = sqrt(2*(E-V)/m);      % velocity
  Tab = pi*sqrt(m/k);       % time interval a to b
  rho = 1./(v.*Tab);        % probability density
% POSITION: Probabilities / expectation values / standard deviation
  Xrho = x.*rho;
  R = 2:N-1;                % v --> then rho --> infinity
  fn = x(R).^2.*rho(R);     % range does not include end points at a and b
  check = simpson1d(rho(R),x(R(1)),x(R(end)));
  xavg  = simpson1d(Xrho(R),x(R(1)),x(R(end)));
  x2avg = simpson1d(fn,x(R(1)),x(R(end)));
  sigmaX = sqrt(x2avg - xavg^2);
% MOMENTUM   p > 0  --> / p < 0 <--
% To calculate pavg need to consider both the motion to the left and right  
  p = m*v;
  fn = p.*rho;
  pavg = simpson1d(fn(R),x(R(1)),x(R(end)));
  pavg = (pavg - simpson1d(fn(R),x(R(1)),x(R(end))))/2;
  fn = p.^2.*rho;
  p2avg = simpson1d(fn(R),x(R(1)),x(R(end)));
  sigmaP = sqrt(p2avg^2 - pavg^2);
%  **** WARNING the standard deviation sigmaP is wrong: ans should be
%  sqrt(2) not 2

% OUTPUT 
  fprintf('Check normalization = %2.4f  \n', check)
  fprintf('Average position, <x> = %2.4f  \n', xavg)
  fprintf('Standard deviation position, sigmaX = %2.4f  \n', sigmaX)
  fprintf('THEORY: standard deviation position, sigmaX = %2.4f  \n', A/sqrt(2))
  fprintf('Average momentum, <p> = %2.4f  \n', pavg)
  fprintf('**** Standard deviation momentum,  sigmaP = %2.4f  \n', sigmaP)
  fprintf('THEORY: standard deviation momentum, sigmaP = %2.4f  \n', sqrt(m*E))
  disp('**** WARNING the standard deviation sigmaP is wrong: ans should be sqrt(2) not 2')

% GRAPHICS
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.2 0.50]);
  set(gcf,'color','w');
  FS = 14;
subplot(3,1,1)
  xP = x; yP = V;
  plot(xP,yP,'b','LineWidth',2)
  %xticks([L1:1:L2])
  grid on
  xlabel('x')
  ylabel('V')
  title('potential energy V(x)','FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(3,1,2)
  xP = x; yP = rho;
  plot(xP,yP,'b','LineWidth',2)
  %xticks([L1:1:L2])
  grid on
  xlabel('x')
  ylabel('\rho')
  set(gca,'FontSize',FS)
  title('probability density \rho(x)','FontWeight','normal')
subplot(3,1,3)
  xP = x; yP = p;
  plot(xP,yP,'b','LineWidth',2)
  %xticks([L1:1:L2])
  grid on
  xlabel('x')
  ylabel('p')
  title('momentum p(x)','FontWeight','normal')
  set(gca,'FontSize',FS)  
