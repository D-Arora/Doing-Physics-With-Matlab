% cemVectorsS.m

% Ian Cooper
% Email      matlabvisualphysics@gmail.com
% Website    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation   https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cem0010.pdf

% Orthogonal components of a [3D] vector
% Cartsian   X R Z
% Polar / cylindrical / pherical
% R  RHO  PHI (azimuthal angle)  THETA (elevation)  
% Date  220708    Matlab version  2021b   

clear 
close all
clc

% Cartesian to polar and cylindrical
% >>>>>  Enter the Cartesian components of the vector V
  V = [3; 5; -6];

  X = V(1); Y = V(2); Z = V(3);
  
  [PHI, RHO] = cart2pol(X, Y);
    disp('Cartesian to polar')  
    rho = RHO
    phi = PHI
    phiD = rad2deg(phi)
 
  [PHI, RHO, Z] = cart2pol(X, Y, Z);
  disp('Cartesian to cylindrical')
    rho = RHO
    phi = PHI
    phiD = rad2deg(phi)
    Vz = Z

 %  Polar / Cylindrical to Cartesian
% >>>>>  Enter the cylindrical components of the vector V
  V = [1.0304, 5.831, -6];

  PHI = V(1); RHO = V(2); Z = V(3);
  
  [X, Y, Z] = pol2cart(PHI, RHO, Z);
   disp('cylindrical to Cartesian')
    Vx = X
    Vy = Y
    Vz = Z

% Cartesian to spherical
% >>>>>  Enter the Cartesian components of the vector V
  V = [3; 5; -6];

  X = V(1); Y = V(2); Z = V(3);
    
  [PHI, THETA, R] = cart2sph(X, Y, Z);
  disp('Cartesian to spherical')
    R
    phi = PHI
    phiD = rad2deg(phi)
    theta = pi/2 - THETA
    thetaD = rad2deg(theta)

 % Spherical to Cartesian
% >>>>>  Enter the spherical components of the vector V
  V = [8.3666, 2.3705, 1.0304];
  
  R = V(1);  THETA = pi/2 - V(2); PHI = V(3);
    
  disp('spherical to Cartesian')
  [X, Y, Z] = sph2cart(PHI, THETA, R);
    Vx = X
    Vy = Y
    Vz = Z
    
  


