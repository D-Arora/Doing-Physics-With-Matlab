close all;
clear all;
clc;

% Circular apertures - irradiance in XY observation plane
% ZONE PLATE
% Irradiance pattern has to be to be symmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]
% Uses functions
%     simpson1d.m  fn_distancePQ.m
% 07 nov 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic

% =======================================================================
% INPUTS 
% =======================================================================

fzone = 1;            % focal length of zone plate [m]
Nzones = 16;          % Number of zones  opaque / transparent

n1 = 200;             % Aperture: grid points for inner ring
n2 = 1200;             % Aperture: grid points for outer ring
nR = 1000;             % Aperture: number of rings
nP = 1309;             % Screen (Observation plane XY): must be odd
wL = 632.8e-9;        % wavelength [m]
uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]
zP = 0.2*fzone;         % Aperture to Screen distance  [m]
yP = 0;               % Observation point P
xPmax = 0.1999e-3; 


% Default values
% fzone = 1;            % focal length of zone plate 
% Nzones = 16;          % Number of zones  opaque / transparent
% 
% n1 = 100;             % Aperture: grid points for inner ring
% n2 = 400;             % Aperture: grid points for outer ring
% nR = 100;             % Aperture: number of rings
% nP = 309;             % Screen (Observation plane XY): must be odd
% wL = 632.8e-9;        % wavelength [m]
% uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]
% zP = 1*fzone;         % Aperture to Screen distance  [m]
% yP = 0;               % Observation point P
% xPmax = 1e-3; 


% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

% Aperture Space -------------------------------------------------------
zQ = 0;
A = zeros(nR,1);       % intgeral for each ring in aperture space
n = zeros(nR,1);       % number of points for ring in aperture space
EQ = ones(nR,1);   % Aperture: Electric Field {V/m]

%   Ring structure
%   radius of ring r  [m]     no. of data points around a ring n
%   Greater the circumference of a ring -->  more grid points
%   Width of each ring  dr
%   Total no. grid points for Aperture  nQ

% ZONE PLATE SETUP -------------------------------------------------------
   rz = zeros(Nzones,1);
   for nc = 1 : Nzones
    rz(nc) = sqrt((fzone + nc * wL / 2)^2 - fzone^2); 
   end

   rMax = rz(end);
   rMin = eps;
   r = linspace(rMin, rMax, nR);
   dr = r(2)-r(1);
   m = (n2-n1) / (nR-1);
   b = n2 - m * nR;

   for c = 1 : nR
      n(c) = 2*round(0.5*(m * c + b))+1;
   end

   nQ = sum(n);
 
   flag = zeros(nR,1);  flags = 1;
   for nc = 2 : Nzones
      if flags == 1; flag(r > rz(nc)) = 1; end;
      if flags ~= 1; flag(r > rz(nc)) = 0; end;
      flags = (-1)^(nc-1);
   end
   
   EQ = sqrt(2*uQmax/(cL*nRI*eps0)) .* EQ;  

   EQ = EQ .* ~flag;               % bright centre
   %EQ = EQ .* flag;               % dark centre
% Observation Space -----------------------------------------------------
   xP = linspace(0,xPmax,nP);
   dxP = xP(2)-xP(1);
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
      
      f = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cP) = dr * sum(A)/ (2*pi);
end

uP = (cL*nRI*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));
uPdB = 10 .* log10(uP./uPmax);

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
    fprintf('Number of zones  =  %3.3d \n',Nzones);
    fprintf('radius of zone plate [m]  =  %3.3e \n',r(end));
    disp('  ');
    disp('Observation Space');
    fprintf('max radius rP [m] =  %3.3e \n',xPmax);
    fprintf('focal length of zone plate  [m]   f = %3.3e \n',fzone);
    fprintf('distance aperture to observation plane [m]   zP = %3.3e \n',zP);
    fprintf('Energy density at xP = 0 yP = 0 zP   [W/m^2m]   uPmax = %3.3e \n',uPmax);
    disp('   ')

    disp('Radial coordinates - zero positions in energy density');
       for c = indexMin
           fprintf('     %3.3e \n',xP(c));
       end

   disp('   ');
   disp('Radial coordinates - max positions in energy density');
   disp('Relative intensities of peaks      ')
       for c = indexMax
           fprintf('     %3.3e      %3.4e \n',xP(c), uP(c)./uPmax);
       end
    
% =======================================================================
% GRAPHICS
% =======================================================================

% [3D] plots --------------------------------------------------------
% Observation space
  rP = xP;
  t = linspace(0,2*pi,nP);

  xx = zeros(nP,nP);
  yy = zeros(nP,nP);

  for c = 1: nP
    xx(c,:) = rP .* cos(t(c));
    yy(c,:) = rP .* sin(t(c));
  end

  uPxy = meshgrid(uP,uP);

% Aperture sapce
  t = linspace(0,2*pi,nR);
  xxQ = zeros(nR,nR);
  yyQ = zeros(nR,nR);

  for c = 1: nR
    xxQ(c,:) = r .* cos(t(c));
    yyQ(c,:) = r .* sin(t(c));
  end
  
  uQ = EQ .* EQ;
  uQxy = meshgrid(uQ,uQ);  


% IRRADIANCE variation in radial direction ------------------------------   
figure(1)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.6 0.6]);
   fs = 16;
   tx = 'radial coordinate r_P  [m]';
   ty = 'energy density  u  [W.m^{-2}]';
   x = xP; y = uP;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
   set(gca,'fontsize',fs);
  % set(gca,'Ylim', [0, 0.2]);
   
% ZONE PLATE IMAGE ------------------------------------------------------
figure(2)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.2 0.2 0.3]);
   pcolor(xxQ,yyQ, uQxy);
   shading interp
   axis equal
   colormap(gray)
   axis off
   set(gcf,'color','k')      
% Diffraction patterns XY plane ------------------------------------------
figure(3)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.3 0.2 0.2 0.3]);
   pcolor(xx,yy,10.*(uPxy).^0.3);
   shading interp
   axis equal
   colormap(gray)
   axis off
   set(gcf,'color','k')

figure(4)
   set(gcf,'Units','normalized');
  set(gcf,'Position',[0.5 0.2 0.2 0.3]);
   surf(xx,yy,10.*(uPxy).^0.3);
   shading interp
   colormap(jet)
  axis off
  set(gcf,'color','b')

% ======================================================================

toc



