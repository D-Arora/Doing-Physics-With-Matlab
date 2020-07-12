% qmTwo.m

% QUANTUM MECHANICS
% Two Electrons in an Infinite Potential Well
% Symmetrical and asymmetrical states
% Fermions and Bosons

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 200429 / Matlab version R2020a

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/qmTwo.htm
% Sripts: Download Location
%    http://www.physics.usyd.edu.au/teach_res/mp/mscripts/



clear
close all
clc

% >>>>> Input principle quantum numbers
  m = 4;
  n = 2;

% SETUP  
  L = 1;
  N = 500;
  k = pi/L;
  x1 = linspace(0,L,N);
  x2 = x1;
  [xx, yy] = meshgrid(x1,x2);

% WAVEFUNCTIONS
  psi1 = sqrt(2/L).*sin(m*k*xx).*sin(n*k*yy);
  psi2 = sqrt(2/L).*sin(m*k*yy).*sin(n*k*xx);

  psiS = psi1 + psi2;   % symmetric state
  psiA = psi1 - psi2;   % antisymmetic state 

% Probability Density  
  probS = psiS.^2;
  probA = psiA.^2;


% GRAPHICS
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.1 0.2 0.5]);
  set(gcf,'color','w');

subplot(2,1,1)
  pcolor(xx,yy,probS)
  shading interp
  set(gca,'xtick',0:0.2:1);
  axis square
  txt = sprintf('Symmetric  |\\psi|^2  m = %2.0f  n = %2.0f', m,n);
  title(txt)
  xlabel('particle 1')
  ylabel('particle 2')
  set(gca,'fontsize',12)
  
subplot(2,1,2)
  pcolor(xx,yy,probA)
  shading interp
  axis square
  set(gca,'xtick',0:0.2:1);
  txt = sprintf('Antisymmetric  |\\psi|^2  m = %2.0f  n = %2.0f', m,n);
  title(txt)
  xlabel('particle 1')
  ylabel('particle 2')
  set(gca,'fontsize',12)


