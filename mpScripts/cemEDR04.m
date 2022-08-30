% cemEDR04.m

% ELECTRIC DIPOLE RADIATION
%   GENERATION OF ELECTROMAGNETIC WAVES
%   Dipole in Z direction

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Notes
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemEDR.htm
% SCRIPTS
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb


close all
clc
clear

% Speed of light
 c = 2.99792458e8;
% wavelength wL / propagation constant k /
% freq f / period / ang freq  w / w*t
  wL =  1;
  k = 2*pi/wL;
  f = c/wL;
  T = 1/f;
  w = 2*pi*f;

% ===================================================================
% E & B FIELD CALCULATIONS  as function of radial displacement
 t = 0e-9;
 th = pi/2;

 N = 599;
 r1 = 1; r2 = 10;
 r = linspace(r1,r2,N);

 kr = k.*r;
 wt = w*t;
 ph = 1i.*(kr - wt);


  Er   = (1 + 1i./kr).*exp(ph)./r.^2;
  Eth  = sin(th).*(1i./kr.^2 + 1./kr - 1i) .* exp(ph)./r;
  Bphi = sin(th).*(1./kr - 1i).*exp(ph)./r;

  
% GRAPHICS
figure(1)
  set(gcf,'Units','normalized');
  set(gcf,'Position',[0.1 0.1 0.20,0.50])
  
subplot(3,1,1)
  xP = r; yP = real(Eth);
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = 0.99./r;
  plot(xP,yP,'r')
  grid on; box on
  ylabel('E_\theta  [ a.u. ]')
  xticks(0:2:10)
  set(gca,'FontSize',14)

subplot(3,1,2)
  xP = r; yP = real(Er);
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = 0.99./r.^2;
  plot(xP,yP,'r')
  grid on; box on
  ylabel('E_r  [ a.u. ]')
  xticks(0:2:10)
  set(gca,'FontSize',14)

subplot(3,1,3)
  xP = r; yP = real(Bphi);
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = 0.99./r;
  plot(xP,yP,'r')
  grid on; box on
  xlabel(' r  [ m ]'); ylabel('B_\phi  [ a.u. ]')
  xticks(0:2:10)
  set(gca,'FontSize',14)  


% ===================================================================
% E & B FIELD CALCULATIONS  as function of polar angle theta
 t = 0e-9;
 r = 8.15; %4.15; %8.15;

 N = 599; th1 = 0; th2 = pi;
 th = linspace(th1,th2,N);

 kr = k.*r;
 wt = w*t;
 ph = 1i.*(kr - wt);

 Er   = ones(1,N).*(1 + 1i./kr).*exp(ph)./r.^2;
 Eth  = sin(th).*(1i./kr.^2 + 1./kr - 1i) .* exp(ph)./r;
 Bphi = sin(th).*(1./kr - 1i).*exp(ph)./r;

% GRAPHICS
figure(2)
  set(gcf,'Units','normalized');
  set(gcf,'Position',[0.50 0.1 0.20,0.50])
  
subplot(3,1,1)
  hold on
  xP = th./pi; yP = real(Eth);
  plot(xP,yP,'b','linewidth',2)
  grid on; box on
  ylabel('E_\theta  [ a.u. ]')
  xticks(0:0.25:1)
  set(gca,'FontSize',14)

subplot(3,1,2)
  hold on
  xP = th./pi; yP = real(Er);
  plot(xP,yP,'b','linewidth',2)
  grid on; box on
  ylabel('E_r  [ a.u. ]')
  xticks(0:0.25:1)
  ylim([0 0.2])
  set(gca,'FontSize',14)

subplot(3,1,3)
  hold on
  xP = th./pi; yP = real(Bphi);
  plot(xP,yP,'b','linewidth',2)
  grid on; box on
  xlabel(' \theta  [ rad / \pi ]'); ylabel('B_\phi  [ a.u. ]')
  xticks(0:0.25:1)
  set(gca,'FontSize',14) 

