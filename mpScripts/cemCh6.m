% cemCh6.m

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

% Load numerical constants
% amu amu_e c e epso h hbar ke me me_e me_u mn mn_e mn_u mp mp_e mp_u mu0
% Do not use these symbols as variables
% SI units for all quantities
  clear; close all; clc
  constantsEM

% ==================================================================
% Electric field surrounding a sphere of radius a and charge Q
  a = 0.25; Q = 1;
% Grid
  N = 999; rMin = 0.0001; rMax = 1;
  r1 = linspace(rMin,rMax, N);
% Calculate electric field 
   E1 = ke.*Q./r1.^2;
   E1(r1<a) = 0;

% ==================================================================
% Electric field surrounding a long straight wire - linear density lambda
%  of radius a and charge Q
  lambda = 1;
% Grid
  N = 999; rMin = 0.5*a; rMax = 1;
  r2 = linspace(rMin,rMax, N);
% Calculate electric field 
   E2 = 2*ke.*lambda./r2;

% Graphics
  figure(1)
    FS = 14;
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.30,0.25])
    xP = r1; yP = E1;
    plot(xP,yP,'b','linewidth',2)
    hold on
    xP = r2; yP = E2;
    plot(xP,yP,'r','linewidth',2)
    grid on; box on
    xlabel('r  [ m ]'); ylabel('E [ N.C^{-1} ]')
    legend('sphere','wire','Orientation','horizontal','Location','best')
    set(gca,'fontsize',FS)

% ELECTRIC FIELD SURROUNDING A SHORT CHARGED STRAIGHT WIRE
% XY Grid for detector points
    N = 5;  L = 2;
    xMin = -L; xMax = L; yMin = -L; yMax = L;
    x = linspace(xMin,xMax,N); y = linspace(yMin,yMax,N);
    [xx, yy] = meshgrid(x,y);
% Source: charged elemenets: charge / linear density / position
   Q = 1; LE = 1; lambdaE = Q/LE;
   NE = 9;
   xE = 0; yEmin = -LE; yEmax = LE;
   yE = linspace(yEmin,yEmax,NE);
   dyE = yE(2) - yE(1);
   QE = lambdaE* dyE;
% Displacements element to detector points
   X = xx;
   for n = 1: NE
       Y = yy - yE(n);
   end
   R = sqrt(X.^2 + Y.^2); R(R<0.1) = 0.1;
   R3 = R.^3;
  
   Ex = ke.* QE.*X./R3; Ey = ke.*QE.*Y./R3;
   E = sqrt(Ex.^2 + Ey.^2);

% GRaphics
figure(2)
quiver(X,Y,Ex,Ey)

figure(3)
contourf(X,Y, E)

%%
clear;  close all; clc
n = rand(100000,1);
n(n<=0.5) = 0;
n(n>0.5) = 1;
n';
sum(n)



