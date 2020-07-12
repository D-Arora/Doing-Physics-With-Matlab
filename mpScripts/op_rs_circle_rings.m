close all;
clear all;
clc;

% Circular apertures - irradiance in XY observation plane
% Irradiance pattern has to be to be symmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]
% Calculation of energy enclosed in circles
% Uses functions
%     simpson1d.m  fn_distancePQ.m   turningPoints.m

% 28 oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
n1 = 100;             % Aperture: grid points for inner ring
n2 = 500;            % Aperture: grid points for outer ring
nR = 501;             % Aperture: number of rings    must be ODD
nP = 509;             % Screen (Observation plane XY): must be odd
wL = 632.8e-9;        % wavelength [m]
a = 1*1e-4;             % radius of circular aperture  [m]

uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]

zP = 1;               % Aperture to Screen distance  [m]
%zP = 6.5*a;
yP = 0;               % Observation point P
xPmax = 2e-2;         % Max radial distance from Z axis
%xPmax = 1.8*a;
% Default values
% n1 = 100;             % Aperture: grid points for inner ring
% n2 = 1000;            % Aperture: grid points for outer ring
% nR = 1000;             % Aperture: number of rings
% nP = 509;             % Screen (Observation plane XY): must be odd
% wL = 632.8e-9;        % wavelength [m]
% a = 1e-4;             % radius of circular aperture  [m]
% uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]
% zP = 1;               % Aperture to Screen distance  [m] 
% yP = 0;               % Observation point P
% xPmax = 2e-2;         % Max radial distance from Z axis


% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

EQmax = sqrt(2*uQmax/(cL*nRI*eps0));   % Aperture: Electric Field {V/m]

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

d_RL = 4*a^2/wL;          % Rayleigh distance  

% Aperture Space -------------------------------------------------------
zQ = 0;
A = zeros(nR,1);     % intgeral for each ring in aperture space
n = zeros(nR,1);     % number of points for ring in aperture space

UQ = uQmax * pi * a^2;   % Theoretical energy emitted from aperture  [J/s]


% Ring structure
%   radius of ring r  [m]     no. of data points around a ring n
%   Greater the circumference of a ring -->  more grid points
%   Width of each ring  dr
%   Total no. grid points for Aperture  nQ
rMax = a;
rMin = eps;
r = linspace(rMin, rMax, nR);
dr = r(2)-r(1);
m = (n2-n1) / (nR-1);
b = n2 - m * nR;

for c = 1 : nR
   n(c) = 2*round(0.5*(m * c + b))+1;
end
nQ = sum(n);

% energy emitted from aperture [J/s]
  
  for c = 1 : nR
     f = uQmax .* ones(1,n(c)); 
     UQring(c) = r(c) * dr * simpson1d(f,0,2*pi);
  end
  UQtheory = sum(UQring);
  
  
% Observation Space -----------------------------------------------------
xP = linspace(0,xPmax,nP);
dxP = xP(2)-xP(1);
% optical coordinates
vP = (2*pi*a/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
EP = zeros(1,nP);


% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    uP energy density (irradiance) along the X axis  [W/m^2]
%    UP energy from aperture to screen  [J/s]

for cP = 1 : nP
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t);
      unit = ones(1,n(cQ));
      rPQ = fn_distancePQ(xP(cP),yP,zP,xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      MP2 = zP .* (ik .* rPQ - unit);
      
      f = EQmax .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cP) = dr * sum(A)/ (2*pi);
end

uP = (cL*nRI*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));
uPdB = 10 .* log10(uP./uPmax);

% Energy received - XY Observation plane  -------------------------------
%   AP      energy within a ring  [J/s]
%   UPsum   energy within XY Observation plane  [J/s]
%   UPsum   energy enclosed within rings of increasing radius  
AP    = zeros(1,nP); 
UPsum = zeros(1,nP);
for cP = 1 : nP
    AP(cP) = 2 * pi * xP(cP) * uP(cP) * dxP ;
end
UP = sum(AP);

for cP = 1 : nP-1
 UPsum(cP+1) = UPsum(cP) + AP(cP);
end


% ======================================================================
% Turning Point Indices
   xData = xP; yData = uP;
   [indexMin, indexMax] = turningPoints(xData, yData) ;
% Zeros
   xP_zeros = xP(indexMin);
% Maxima
   xP_maxs = xP(indexMax);
   uP_maxs = uP(indexMax);
   

% =======================================================================
% OUTPUT to Command Window
% =======================================================================
disp('Parameter summary  [SI units]');
fprintf('wavelength [m]  =  %3.5g \n',wL);
fprintf('nQ  =  %3.3d \n',nQ);
fprintf('nP  =  %3.3d \n',nP);
disp('  ')
disp('Aperture Space');
fprintf('radius of aperture [m]  =  %3.3e \n',a);
fprintf('energy density [W/m2] uQmax  =  %3.3e \n',uQmax);
fprintf('energy from aperture [J/s]   UQ(theory) = %3.3e \n',UQ);
disp('  ')
disp('Observation Space');
fprintf('max radius rP [m] =  %3.3e \n',xPmax);
fprintf('distance aperture to observation plane [m]   zP = %3.3e \n',zP);
fprintf('Rayleigh distance  [m]   d_RL = %3.3e \n',d_RL);
disp('  ');
fprintf('energy: aperture to screen  [J/s]   UP = %3.3e \n',UP);
fprintf('max energy density  [W./m2]   uPmax = %3.3e \n',uPmax);
disp('   ');
disp('   ');

disp('Radial coordinates - zero positions in energy density');
for c = indexMin
   fprintf('     %3.3f \n',vP(c));
end

disp('   ');
disp('Radial coordinates - max positions in energy density');
disp('Relative intensities of peaks      ')
for c = indexMax
   fprintf('     %3.3f      %3.4f \n',vP(c), uP(c)./uPmax);
end

% =======================================================================
% GRAPHICS
% =======================================================================
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 14;

% energy density   uP vs Xp & yP
   subplot(2,2,1);
   tx = 'radial coordinate r_P  [m]';
   ty = 'energy density  u  [W.m^{-2}]';
   x = xP; y = uP;
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);

% energy density   uPdB  vs xP & yP
   subplot(2,2,3);   
   tx = 'radial coordinate  r_P  [m]';
   ty = 'energy density  u  [dB]';
   x = xP;   y = uPdB;
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);

% energy density   uP vs vPx & vPy
   subplot(2,2,2);
   tx = 'optical coordinate  v_P';
   ty = 'energy density  u  [a.u.]';
   x = vP;   y = uP./uPmax;
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);
   grid on
   set(gca,'Xlim',[0,20]);
% energy density   uPdB  vs xP & yP
   subplot(2,2,4);   
   tx = 'optical coordinate   v_P';
   ty = 'energy density  u  [dB]';
   x = vP;   y = uPdB;
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);
   set(gca,'Xlim',[0,20]);
   grid on

% Energy enclosed --------------------------------------------------------
figure(2)
   tx = 'optical coordinate  v_P ';
   ty = '% energy enclosed within circle of radius v_P ';
   x = vP; y = 100 * UPsum./UQ;
   plot(x,y,'lineWidth',2);
   xlabel(tx);   ylabel(ty);
   grid on

figure(22)
   tx = 'radial coordinate  x_P ';
   ty = '% energy enclosed within circle of radius x_P ';
   x = xP; y = 100 * UPsum./UQ;
   plot(x,y,'lineWidth',2);
   xlabel(tx);   ylabel(ty);
   hold on
   x = [a, a]; y = [0,100];
   plot(x,y,'r','lineWidth',1);
   grid on
   

% [3D] plots --------------------------------------------------------
  r = xP;
  t = linspace(0,2*pi,nP);

  xx = zeros(nP,nP);
  yy = zeros(nP,nP);

  for c = 1: nP
    xx(c,:) = r .* cos(t(c));
    yy(c,:) = r .* sin(t(c));
  end

  uPxy = meshgrid(uP,uP);

  
figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.3]);
pcolor(xx,yy,10.*(uPxy).^0.3);
shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')

figure(4)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.4 0.2 0.2 0.3]);
surf(xx,yy,10.*(uPxy).^0.3);
shading interp
colormap(jet)
axis off
set(gcf,'color','b')


figure(5)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 16;

% energy density   uP vs Xp & yP
   tx = 'radial coordinate r_P  [m]';
   ty = 'energy density  u  [W.m^{-2}]';
   x = xP; y = uP;
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);

% =======================================================================




% figure(1)
% for c = 1 : nR
%    t = linspace(0,2*pi,n(c));  
%    x = r(c) .* cos(t);
%    y = r(c) .* sin(t); 
%    plot(x,y,'o');
%    axis square
%    hold on
% end

toc



