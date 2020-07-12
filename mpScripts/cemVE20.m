% cemVE20.m
% 98 May 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Electostatics - point charge
%    partial derivatives
%    divergence and curl
% SI units used unless stated otherwise

clear all
close all
clc

tic

% ========================================================================
% INPUTS  
% ========================================================================

% Value of the point charge Q located at the origin (0,0,0)
   Q = 3e-7;
% Increment h = dx = dy = dz
   h = 0.0001;
% Cartesian coordinates (x,y,z) for the point P at which the
%     partial derivatives,divergence and curl are calculated
   p = [0.71, 0.71, 0.71];
  
% =======================================================================
% SETUP 
% =======================================================================
% Cartesian coordinates of the point P and increments around P   
   x = p(1); x2 = x + h/2; x1 = x - h/2;
   y = p(2); y2 = y + h/2; y1 = y - h/2;
   z = p(3); z2 = z + h/2; z1 = z - h/2;

% Constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);
   k = kC * Q;
   
% Matrix for the 9 partial derivates of the electric field
   gradE = zeros(3,3);

% Distances from the origin (0,0,0) to the point P and increments around P 
   R   = sqrt(x^2 + y^2 + z^2);
   Rx2 = sqrt(x2^2 + y^2 + z^2);
   Rx1 = sqrt(x1^2 + y^2 + z^2);
   Ry2 = sqrt(x^2 + y2^2 + z^2);
   Ry1 = sqrt(x^2 + y1^2 + z^2);
   Rz2 = sqrt(x^2 + y^2 + z2^2);
   Rz1 = sqrt(x^2 + y^2 + z1^2);
   
   rho = Q / h^3;
   rho_eps0 = rho / eps0;

% =======================================================================
% CALCULATIONS: partial derivates / divergence / curl 
% =======================================================================

% dEx/dx
   E2 = k * x2 / Rx2^3;
   E1 = k * x1 / Rx1^3;
   gradE(1,1)  = (E2 - E1) / h;
   gradEA(1,1) = k * (1/R^3 - 3*x^2/R^5);
   
% dEx/dy;
   E2 = k * x / Ry2^3;
   E1 = k * x / Ry1^3;
   gradE(1,2)  = (E2 - E1) / h;
   gradEA(1,2) = -k * (3*x*y/R^5);
   
% dEx/dz;
   E2 = k * x / Rz2^3;
   E1 = k * x / Rz1^3;
   gradE(1,3)  = (E2 - E1) / h;
   gradEA(1,3) = -k * (3*x*z/R^5);
   
% dEy/dx
   E2 = k * y / Rx2^3;
   E1 = k * y / Rx1^3;
   gradE(2,1)  = (E2 - E1) / h;
   gradEA(2,1) = -k * (3*x*y/R^5);
   
% dEy/dy;
   E2 = k * y2 / Ry2^3;
   E1 = k * y1 / Ry1^3;
   gradE(2,2)  = (E2 - E1) / h;
   gradEA(2,2) = k * (1/R^3 - 3*y^2/R^5);
   
% dEy/dz;
   E2 = k * y / Rz2^3;
   E1 = k * y / Rz1^3;
   gradE(2,3)  = (E2 - E1) / h;
   gradEA(2,3) = -k * (3*y*z/R^5);

% dEz/dx
   E2 = k * z / Rx2^3;
   E1 = k * z / Rx1^3;
   gradE(3,1)  = (E2 - E1) / h;
   gradEA(3,1) = -k * (3*x*z/R^5);
% dEz/dy;
   E2 = k * z / Ry2^3;
   E1 = k * z / Ry1^3;
   gradE(3,2)  = (E2 - E1) / h;
   gradEA(3,2) = -k * (3*y*z/R^5);
   
% dEz/dz;
   E2 = k * z2 / Rz2^3;
   E1 = k * z1 / Rz1^3;
   gradE(3,3)  = (E2 - E1) / h;
   gradEA(3,3) = k * (1/R^3 - 3*z^2/R^5);
  
% divergence of the electric field at the point P(x,y,z)
   divE  = gradE(1,1)  + gradE(2,2)+  gradE(3,3);
   divEA = gradEA(1,1) + gradEA(2,2)+ gradEA(3,3);

% curl of the electric field at the point P(x,y,z)
   curlEx = gradE(3,2) - gradE(2,3);
   curlEy = gradE(1,3) - gradE(3,1);
   curlEz = gradE(2,1) - gradE(1,2);

% =======================================================================
% COMMAND WINDOW OUPUT: partial derivates / divergence / curl 
% =======================================================================
disp('   ');
fprintf('Charge at origin  Q = %2.3e  C \n',Q);
disp('   ');
fprintf('Charge density  rho = %2.3e  C/m^3 \n',rho);
disp('   ');
fprintf('  rho / eps0 = %2.3e  V/m^2 \n',rho_eps0);
disp('   ');
disp('Observation point P(x,y,z)');
fprintf('   x = %2.3e  m ',p(1));
fprintf('   y = %2.3e  m ',p(2));
fprintf('   z = %2.3e  m \n ',p(3));
disp('   ');
disp('Displacement increment h = dx = dy = dz  ');
fprintf('   h = %2.3e  m \n',h);
disp('  ');
disp('Partial derivatives: numerical calculations');
gradE
disp('   ');
disp('Partial derivatives: analytical calculations');
gradEA
disp('Divergence of E    ');
fprintf('   divE = %2.3e  m \n',divE); 
disp('  ');
fprintf('   divEA = %2.3e  m \n',divEA); 
disp('  ');
disp('Curl of E    ');
fprintf('   curlEx = %2.3e  m \n',curlEx); 
fprintf('   curlEy = %2.3e  m \n',curlEy); 
fprintf('   curlEz = %2.3e  m \n',curlEz); 
disp('  ');

  toc
