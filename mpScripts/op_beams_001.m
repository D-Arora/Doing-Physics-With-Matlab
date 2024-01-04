% op_beams_001.m

% Modelling scalar GAUSSIAN BEAMS
%   Assume a Gaussian profile with spherical wavefronts in the paraxial
%   regeme.

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op_beams_001.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  190131

close all
clear
clc

% INPUTS ==============================================================

% Max power transmitted in beam  [W]
P0 = 1e-3;

% beam waist  [m]
w0 = 0.5e-3;

% wavelength  [m]    [HENe laser  wL = 632.8e-9 m]
wL = 632.8e-9;

% Length of Z domain in multiples of zR :  z = nZ * zR
nZ = 5;

% Number of grid points in calculations  (must be an odd number)
N = 501;

% Z position of XY plane for radial irradiance plots:  zPR = zP / zR
%   [e,g,  zPR = 0, 1, 2, ... ]
zPR = 0;

% Max radial distance for radial irradiance plot: Rmax = nR * w0
nR = 5;


% CALCULATIONS  =======================================================

% Wave number
   k = 2*pi/wL;
   
% Speed of light
   c = 299792458;

% Permittivity of free space
   eps0 = 8.85418782e-12;

% Electric field amplitude  [V/m]
   E0 = sqrt(4*P0/(pi*c*eps0*w0^2));

% Max beam irradiance  [W/m^2]
  Smax = 2*P0/(pi*w0^2);
   
% Z Domain  z_zR = z / zR
   z_zR = linspace(0,nZ,N);
   
% Rayleigh range  [m]
   zR = pi*w0^2/wL;
   
% Z domain  z  [m]   
    z = zR .* z_zR;
   
% Beam divergence angle [rad]  [deg]
   thetaR = w0/zR;
   theta = rad2deg(thetaR);
  
% Beam spot  [m]
   w = w0 .* sqrt(1+z_zR.^2);
   
% Gouy Phase [rad/pi]
   phi = atan(z_zR)./pi;
    
% Radius of curvature of wavefront  [m]:   Rmax = 20
   R = (z + eps) + zR^2./(z + eps);
   R(R>20) = 20;
  
% Axial Irradiance  [W/m^2]
   Sz = (c*eps0*E0^2/2) ./ (1+z_zR.^2);

% Radial Irradiance  [W/m^2]
   rP = nR * w0;
   wP = w0 * sqrt(1+(zPR)^2);
   r = linspace(0,rP,N);
   Sr = ((c*eps0*E0^2/2)./(1+(zPR)^2)) .* exp(-2.*(r./wP).^2);
    
% Irradiance in XY plane  [W/m^2]
   x = linspace(-rP,rP,N);
   [xx, yy] = meshgrid(x,x);
   r2 = xx.^2 + yy.^2;
   K = -2.* r2 ./ wP^2;
   Sxy = (c*eps0*E0^2/2).*(w0/wP) .* exp(K);
   
% Irradiance in R-Z plane  [W/m^2]
   [zzR, xxR] = meshgrid(z,x);
   r2 = xxR.^2;
   wR = w0.*sqrt(1 + (zzR/zR).^2);
   KR = -2.* r2 ./ wR.^2;
   Szx = (c*eps0*E0^2/2).*(w0./wR) .* exp(KR);   
   

% Power through a circular aperture placed in beam at the waist  z = 0
   rA = linspace(0,rP,N);
   P = P0 .* (1 - exp(-2.*(rA./w0)).^2);
   

% phase of the wave along the Z axis
   E_chi = E0 .* (1./sqrt(1+z_zR.^2)) .* exp(1i*(k.*zR.*z_zR - atan(z_zR)));
   chi = angle(E_chi);

   
% GRAPHICS  =========================================================== 
  FS = 12; 
  col = ColorCode(wL);
  myMap = zeros(64,3);
     for c1 = 1:64
      myMap(c1,:) = col;
     end
     
figure(1)
  pos = [0.1 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = z_zR; yP = w.*1e3;
    plot(xP,yP,'color',col,'linewidth',2)

  hold on

  xP = [1 1]; yP = [0 3.5];
    plot(xP,yP,'k','linewidth',1)
    
  xP = [0 5]; yP = [0, 0.5 + 1e3*(5*wL/(pi*w0))];
    plot(xP,yP,'m','linewidth',1)

tm1 = num2str(zR,'z_R = %2.4f  m   ');
tm2 = '    \theta = ';
tm3 = num2str(theta,' %2.4f  deg');
tm = [tm1 tm2 tm3];
title(tm)

grid on
set(gca,'fontsize',12)
xlabel('z / z_R')
ylabel('beam spot     w   [ mm ]')
ylim([0 3.5])


figure(2)
  pos = [0.4 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = z_zR; yP = Sz;
    plot(xP,yP,'color',col,'linewidth',2)
    
   hold on

  xP = [1 1]; yP = [0 max(Sz)];
    plot(xP,yP,'k','linewidth',1)  

  grid on
  set(gca,'fontsize',12)
  xlabel('z / z_R')
  ylabel('axial irradiance  S_z   [ W / m^2 ]')  
  
  tm1 = num2str(zR,'z_R = %2.4f  m   ');
  tm2 = '    S_{max} = ';
  tm3 = num2str(Smax,' %2.2e  W.m^{-2}');
  tm = [tm1 tm2 tm3];
  title(tm)
  

  figure(3)
  pos = [0.7 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = r./w0; yP = Sr;
    plot(xP,yP,'color',col,'linewidth',2)

  grid on
  set(gca,'fontsize',12)
  xlabel('r / w_0')
  ylabel('radial irradiance  S_r    [ W / m^2 ]') 
  
  tm1 = num2str(zPR,'z / z_R = %2.2f');
  tm2 = '    w_0 = ';
  tm3 = num2str(1e3*w0,' %2.2f  mm');
  tm = [tm1 tm2 tm3];
  title(tm)
  
figure(4)
  pos = [0.1 0.5 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = z_zR; yP = phi;
    plot(xP,yP,'color',col,'linewidth',2)
    
 hold on
 
 xP = -z_zR; yP = -phi;
    plot(xP,yP,'color',col,'linewidth',2)
    
 xP = [1 1]; yP = [-0.5 0.5];
    plot(xP,yP,'k','linewidth',1)   
  
 grid on
 set(gca,'fontsize',12)
 set(gca,'xtick',-nZ:nZ)
 set(gca,'ytick',-0.5:0.25:0.5);
 xlabel('z / z_R')
 ylabel('\phi / \pi')   
 title('Gouy Phase Angle')
 
 
 figure(5)
  pos = [0.4 0.5 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = z_zR; yP = R;
    plot(xP,yP,'color',col,'linewidth',2)
    
   hold on

  xP = [1 1]; yP = [0 20];
    plot(xP,yP,'k','linewidth',1)  
    
  grid on
  set(gca,'fontsize',12)
  xlabel('z / z_R')
  ylabel('radius of curvature  R   [ m ]')   
  tm = num2str(zR,'z_R = %2.4f  m');
  title(tm)
  

 figure(6)
  pos = [0.25 0.2 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = rA/w0; yP = 100*P/P0;
    plot(xP,yP,'color',col,'linewidth',2)
  
  hold on
  
  xP = [1 1]; yP = [0 100];
    plot(xP,yP,'k','linewidth',1) 
    
  grid on
  set(gca,'fontsize',12)
  xlabel('radius / w_0')
  ylabel('% power through aperture') 
  
  tm = num2str(1e3*w0,' %2.2f  mm');
  title(tm)  
  
  
figure(7)
  pos = [0.45 0.2 0.25 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  pcolor(zzR,xxR.*1e3,Szx)
  shading interp
  shadingMap = gray(64);
  colormap(shadingMap(:,1)*col);
  
  hold on
  xP = z; yP = w0 .* sqrt(1+(z./zR).^2) .* 1e3;
     plot(xP,yP,'y','linewidth',1.5);
     plot(xP,-yP,'y','linewidth',1.5);     
     
  set(gca,'fontsize',12)
  xlabel('z  [m]')
  ylabel('x  [mm]') 
  colorbar
  

figure(8)
  pos = [0.25 0.2 0.25 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  subplot(2,1,1) 
  surf(xx.*1e3,yy.*1e3,Sxy)
  shading interp
  view(-36,58)
  shadingMap = gray(64);
  colormap(shadingMap(:,1)*col);
  set(gca,'fontsize',12)
  xlabel('x  [mm]')
  ylabel('y  [mm]') 
  zlabel('S_{XY}  [a.u.]') 

  tm = num2str(zPR,'z / z_R = %2.2f');
  title(tm)
  
  subplot(2,1,2) 
  pcolor(xx.*1e3,yy.*1e3,Sxy)
  shading interp
  shadingMap = gray(64);
  colormap(shadingMap(:,1)*col);
  set(gca,'fontsize',12)
  xlabel('x  [mm]')
  ylabel('y  [mm]') 
  colorbar
  axis square
  

figure(9)
  pos = [0.1 0.3 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = z_zR; yP = chi./pi;
    plot(xP,yP,'color',col,'linewidth',2)

  tm = num2str(zR,'z_R = %2.4f  m   ');
  title(tm)

  grid on
  set(gca,'fontsize',12)
  xlabel('z / z_R')
  ylabel('phase  \chi / \pi [ rad ]')

  

