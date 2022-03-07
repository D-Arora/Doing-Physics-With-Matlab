% FOP01.m

% FOURIER OPTICS: PLANE WAVE ANALYSIS


close all
clc
clear

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Vacuum wavelength  [m]
  wL0 = 1.55e-6;
% Refractive index of medium
  n = 1.5;

% Direction of polarization
%  phi   (polar angle: -pi/2 to p/2 rad)
%  theta (azimuthal angle: -pi to pi rad)
  phiDeg   = 45;   % [deg];
  thetaDeg = 30;   % [deg]

% Electric field components at Origin [a.u.]
  E0x = 1;
  E0y = 2;
  E0z = -1.86;

% Observation point P [m]  
  x = 1e-6;
  y = 2e-6;
  z = 3e-6;

% Speed of light  [m/s]
  c = 3.00e8;


% CALCULATIONS  ===================================================
%  Angles: deg to rad
   phi   = deg2rad(phiDeg);
   theta = deg2rad(thetaDeg);

% speed of plane wave in medium  [m/s]
  v = c/n;
% wavelength in medium  [m]
  wL = wL0/n;
% Period  [s]
  T = wL / v;
% Temporal frequency  [Hz]
  f = 1/T;
% Angular temporal frequency  [1/s]
  w = 2*pi*f;


% Unit vector e components
  ex = sin(phi)*cos(theta);
  ey = sin(phi)*sin(theta);
  ez = cos(phi);

% Propagation constant in medium  [1/m]
  k = 2*pi/wL;
  kx = k*ex;
  ky = k*ey;
  kz = k*ez;

% Spatial frequencies  [1/m]
  fs = 1/wL; 
  fsx = fs*ex;
  fsy = fs*ey;
  fsz = fs*ez;
% Wavelength components [m]
  wLx = 1/fsx;
  wLy = 1/fsy;
  wLz = 1/fsz;

% OBSERVATION POINT P(x,y,z)
  r = sqrt(x^2+y^2+z^2);
  kr = kx*x + ky*y + kz*z;

  Ex = E0x*exp(1i*kr);
  Ey = E0y*exp(1i*kr);
  Ez = E0z*exp(1i*kr);

  E = sqrt(Ex*conj(Ex) + Ey*conj(Ey) + Ez*conj(Ez));





