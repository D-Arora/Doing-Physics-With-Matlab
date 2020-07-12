close all;
clear all;
clc;

% op_rs_point_source_xz.m
% Circular apertures - irradiance in ZY MERIDONAL observation plane
% POINT SOURCE illumination of aperture
%    point source must be on optical  xS = 0  ys = 0
% a
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

% 13 nov oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
n1 = 50;             % Aperture: grid points for inner ring
n2 = 80;            % Aperture: grid points for outer ring
nR = 101;             % Aperture: number of rings    must be ODD
nP = 121;             % Screen (Observation plane XY): must be odd
wL = 632.8e-9;        % wavelength [m]
a = 10*wL;             % radius of circular aperture  [m]

%uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]

% Source
    xS = 0; yS = 0;
    %zS = -1;
    zS = -1;
    ES = 1;

% Observation Space  [m]    
zPmin = 1*wL;
zPmax = 200*wL; 
xPmin = 0;
xPmax = 10*wL; 
yP = 0;               


% Default values --------------------------------------------------------
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

%EQmax = sqrt(2*uQmax/(cL*nRI*eps0));   % Aperture: Electric Field {V/m]

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

d_RL = 4*a^2/wL;          % Rayleigh distance  

% Aperture Space -------------------------------------------------------
zQ = 0;
A = zeros(nR,1);     % intgeral for each ring in aperture space
n = zeros(nR,1);     % number of points for ring in aperture space

%UQ = uQmax * pi * a^2;   % Theoretical energy emitted from aperture  [J/s]


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

  
% Observation Space -----------------------------------------------------
    zP = linspace(zPmin,zPmax,nP);
    dzP = zP(2)-zP(1);
    xP = linspace(xPmin,xPmax,nP);
    dxP = xP(2)-xP(1);
    [zz,xx] = meshgrid(zP./wL,xP./wL);

    EP = zeros(nP,nP);


% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    uP energy density (irradiance) along the X axis  [W/m^2]
% aperture electric field EQ / energy density / energy transmitted
     UQring = zeros(1,nR);  EQ = zeros(1,nR);
     for cQ = 1 : nR
        rSQ = sqrt(zS^2 + r(cQ)^2);
        EQ(cQ) = ES * exp(ik .* rSQ) ./ rSQ;
     end


for cz = 1 : nP
for cx = 1 : nP    
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t); 
      rPQ = fn_distancePQ(xP(cx),yP,zP(cz),xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      unit = ones(1,n(cQ));
      MP2 = zP(cz) .* (ik .* rPQ - unit);
      
      f = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cx,cz) = dr * sum(A)/ (2*pi);
end
end

uP = (cL*nRI*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));


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

disp('  ')
disp('Observation Space');

disp('  ');
disp('   ');
disp('   ');



% =======================================================================
% GRAPHICS
% =======================================================================
 figure(1)
    fs = 12;
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 14;
    pcolor(zz,xx,uP.^0.8);
    shading interp
    tx = 'axial position  z_P / \lambda';
    ty = 'radial position r_P / \lambda';
    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
    set(gca,'fontsize',fs);
   
   figure(2)
    plot(zP./wL, uP(1,:));
    
toc



