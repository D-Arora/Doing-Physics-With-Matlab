% cemVE20.m
% 98 May 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Electostatics  
%    partial derivatives
%    divergence and curl of an electric field specified by a function
% SI units used unless stated otherwise

clear all
close all
clc

tic

% ========================================================================
% INPUTS  
% ========================================================================

% Increment h = dx = dy = dz
   h = 0.0001;
% Cartesian coordinates (x,y,z) for the point P at which the
%     partial derivatives,divergence and curl are calculated
   p = [2, 0, 0];
  
% =======================================================================
% SETUP 
% =======================================================================

% Define function for electric field
K = 1;
Ex = @(x,R)K .* x .* R.^2;

% Cartesian coordinates of the point P and increments around P   
   x = p(1); x2 = x + h/2; x1 = x - h/2;
   y = p(2); y2 = y + h/2; y1 = y - h/2;
   z = p(3); z2 = z + h/2; z1 = z - h/2;

% Constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);
      
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
   
% =======================================================================
% CALCULATIONS: partial derivates / divergence / curl 
% =======================================================================
   E = K * R^3;
   
% dEx/dx
   E2 = Ex(x2,Rx2);
   E1 = Ex(x1,Rx1);
   gradE(1,1)  = (E2 - E1) / h;
      
% dEx/dy;
    
   E2 = Ex(x2,Ry2); 
   E1 = Ex(x1,Ry1);
   gradE(1,2)  = (E2 - E1) / h;
     
% dEx/dz;
   E2 = Ex(x2,Rz2);
   E1 = Ex(x1,Rz1);
   gradE(1,3)  = (E2 - E1) / h;
    
% dEy/dx
   E2 = Ex(y2,Rx2);
   E1 = Ex(y1,Rx1);
   gradE(2,1)  = (E2 - E1) / h;
    
% dEy/dy;
   E2 = Ex(y2,Ry2);
   E1 = Ex(y1,Ry1);
   gradE(2,2)  = (E2 - E1) / h;
   
% dEy/dz;
   E2 = Ex(y2,Rz2);
   E1 = Ex(y1,Rz1);
   gradE(2,3)  = (E2 - E1) / h;
 
% dEz/dx
   E2 = Ex(z2,Rx2);
   E1 = Ex(z1,Rx1);
   gradE(3,1)  = (E2 - E1) / h;
   
% dEz/dy;
   E2 = Ex(z2,Ry2);
   E1 = Ex(z1,Ry1); 
   gradE(3,2)  = (E2 - E1) / h;
    
% dEz/dz;
   E2 = Ex(z2,Rz2);
   E1 = Ex(z1,Rz1);
   gradE(3,3)  = (E2 - E1) / h;
  
 % divergence of the electric field at the point P(x,y,z)
   divEP  = gradE(1,1)  + gradE(2,2)+  gradE(3,3);
   
 % charge density rho at the point P
   rho = divEP * eps0;
   
% curl of the electric field at the point P(x,y,z)
   curlEx = gradE(3,2) - gradE(2,3);
   curlEy = gradE(1,3) - gradE(3,1);
   curlEz = gradE(2,1) - gradE(1,2);
   
   
% Variation of divergence as a function of radial position R --> x  
   N = 201;     % must be an odd number
   x = linspace(0,10,N);
   ER = K .* x.^3;

   x2 = x + h/2; x1 = x - h/2;
   R   = sqrt(x.^2 + y^2 + z^2);
   Rx2 = sqrt(x2.^2 + y^2 + z^2);
   Rx1 = sqrt(x1.^2 + y^2 + z^2);
   Ry2 = sqrt(x.^2 + y2^2 + z^2);
   Ry1 = sqrt(x.^2 + y1^2 + z^2);
   Rz2 = sqrt(x.^2 + y^2 + z2^2);
   Rz1 = sqrt(x.^2 + y^2 + z1^2);
   
% dEx/dx
   E2 = Ex(x2,Rx2);
   E1 = Ex(x1,Rx1);
   divEx  = (E2 - E1) ./ h;
% dEy/dy;
   E2 = Ex(y2,Ry2);
   E1 = Ex(y1,Ry1);
   divEy  = (E2 - E1) ./ h;
% dEz/dz;
   E2 = Ex(z2,Rz2);
   E1 = Ex(z1,Rz1);
   divEz  = (E2 - E1) ./ h;
   
   divE = divEx + divEy + divEz;
   
   rhoR = eps0 .* divE;

% Charge Q enclosed by a sphere of radius R
  dx = x(2)-x(1); Q(1) = x(1)^2 * rhoR(1); Q = zeros(N,1);
  for n = 2 : N
  Q(n) = Q(n-1) + x(n)^2 * rhoR(n); 
  end
  Q = Q .* (4*pi*dx); 
   
% =======================================================================
% COMMAND WINDOW OUPUT: partial derivates / divergence / curl 
% =======================================================================
disp('   ');
disp('Observation point P(x,y,z)');
fprintf('   x = %2.3e  m ',p(1));
fprintf('   y = %2.3e  m ',p(2));
fprintf('   z = %2.3e  m \n ',p(3));
disp('   ');
disp('Displacement increment h = dx = dy = dz  ');
fprintf('   h = %2.3e  m \n',h);
disp('  ');
disp('Electric field strength at point P');
fprintf('   E = %2.3e  V/m \n',E);
disp(' ');
disp('Partial derivatives: numerical calculations');
gradE
disp('   ');
disp('Divergence of E    ');
fprintf('   divE = %2.3e  m \n',divEP); 
disp('   ');
fprintf('Average charge density at point P  rho = %2.3e  C/m^3 \n',rho);
disp('  ');
disp('Curl of E    ');
fprintf('   curlEx = %2.3e  m \n',curlEx); 
fprintf('   curlEy = %2.3e  m \n',curlEy); 
fprintf('   curlEz = %2.3e  m \n',curlEz); 
disp('  ');

% % ======================================================================= 
% % GRAPHICS 
% % ======================================================================= 
   fs = 12;

figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.03 0.50 0.43 0.35]);
   xP = x; yP = rhoR;
   subplot(1,2,1)
   plot(xP,yP,'linewidth',2);
   tx = 'R  (m)';  ty = '\rho  [C/m^3]';
   xlabel(tx,'fontsize',fs'); ylabel(ty,'fontsize',fs');
   set(gca,'fontsize',fs');
   
   xP = x.^2; yP = rhoR;
   subplot(1,2,2)
   plot(xP,yP,'linewidth',2);
   tx = 'R^2  (m^2)';  ty = '\rho  [C/m^3]';
   xlabel(tx,'fontsize',fs'); ylabel(ty,'fontsize',fs');
   
   set(gca,'fontsize',fs');
   box on

figure(2)   % 22222222222222222222222222222222222222222222222222222222222 
   set(gcf,'units','normalized','position',[0.50 0.50 0.43 0.35]);
   xP = x; yP = Q;
   subplot(1,2,1)   
   plot(xP,yP,'linewidth',2);
   tx = 'R  (m)';  ty = 'Q_{enclosed}  [C]';
   xlabel(tx,'fontsize',fs'); ylabel(ty,'fontsize',fs');
   set(gca,'fontsize',fs');
   
   xP = x.^5; yP = Q;
   subplot(1,2,2)
   plot(xP,yP,'linewidth',2);
   tx = 'R^5  (m^5)';  ty = 'Q_{enclosed}  [C]';
   xlabel(tx,'fontsize',fs'); ylabel(ty,'fontsize',fs');
   
   set(gca,'fontsize',fs');
   box on

figure(3)   % 333333333333333333333333333333333333333333333333333333333333 
   set(gcf,'units','normalized','position',[0.03 0.06 0.43 0.35]);
   xP = x; yP = ER;
   subplot(1,2,1)   
   plot(xP,yP,'linewidth',2);
   tx = 'R  (m)';  ty = 'E  [V/m]';
   xlabel(tx,'fontsize',fs'); ylabel(ty,'fontsize',fs');
   set(gca,'fontsize',fs');
   
   xP = x.^5; yP = Q;
   subplot(1,2,2)
   plot(xP,yP,'linewidth',2);
   tx = 'R^3  (m^3)';  ty = 'E  [V/m]';
   xlabel(tx,'fontsize',fs'); ylabel(ty,'fontsize',fs');
   
   set(gca,'fontsize',fs');
   box on
   
   
   toc
