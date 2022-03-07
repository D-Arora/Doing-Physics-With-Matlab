% em_01.m

% Vector coordinate systems: Cartesian, Cylindrical, Spherical
% CELL 1:  Cartesian  ---> cylindrical & spherical
% CELL 2:  Cylindical ---> Cartesian % spherical
% CELL 3:  Spherical  ---> Cartesian & cylindrical

% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

% DOWNLOAD Scripts from
%   https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%   https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% 220226  Matlab Version R2021b

  
%   Convert Cartesian coordinates of a vector into its
%   cylindrical and spherical coordinates


close all
clear
clc

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Vector Cartesian components
  Ax = -1.5;
  Ay = 0.8;
  Az = 1.3;

% CALCULATIONS ---------------------------------------------------
% Vector A  
  A = [Ax;Ay;Az];
% Magnitude Amag
  Ar = norm(A);
% Polar magnitude Arho
  Arho = sqrt(Ax^2 + Ay^2);
% Polar angle phi  [rad]
   phi = atan2(Ay,Ax);
% Azimuthal angle  theta [rad]
   theta = acos(Az/Ar);

% OUTPUT RESULTS TO COMMAND WINDOW
disp('  ')
fprintf('Ax = %2.4f   Ay = %2.4f   Az = %2.4f  \n',Ax,Ay,Az)
disp('  ')
fprintf('Ar = %2.4f   Arho = %2.4f \n',Ar,Arho)
disp('  ')
fprintf('phi = %2.4f  rad   theta = %2.4f rad  \n',phi, theta)
disp('  ')


phiHx = -sin(phi); phiHy = cos(phi); phiHz = 0;

thetaHx = cos(phi)*cos(theta);
thetaHy = sin(phi)*cos(theta);
thetaHz = -sin(theta);

figure(1)

xP = [0 1]; yP = [0 0]; zP = [0 0];
plot3(xP,yP,zP,'m','linewidth',2)

hold on
xP = [0 0]; yP = [0 1]; zP = [0 0];
plot3(xP,yP,zP,'m','linewidth',2)

xP = [0 0]; yP = [0 0]; zP = [0 1];
plot3(xP,yP,zP,'m','linewidth',2)

xP = [0 Ax]; yP = [0 Ay]; zP = [0 Az];
plot3(xP,yP,zP,'r','linewidth',3)

xP = [0 Ax/Ar]; yP = [0 Ay/Ar]; zP = [0 Az/Ar];
plot3(xP,yP,zP,'g','linewidth',3)

xP = [0 phiHx]; yP = [0 phiHy]; zP = [0 phiHz];
plot3(xP,yP,zP,'b','linewidth',3)

xP = [0 thetaHx]; yP = [0 thetaHy]; zP = [0 thetaHz];
plot3(xP,yP,zP,'k','linewidth',3)
% xP = [0 1]; yP = [0,0]; zP = [1 1];
% plot3(xP,yP,zP,'k','linewidth',0.2)
% xP = [0 0]; yP = [0,0]; zP = [0 1];
% plot3(xP,yP,zP,'k','linewidth',0.2)


xticks(-1:1)
yticks(-1:1)
zticks(-1:1)

title('A_r (green)  \phi (blue)  \theta (black)')
set(gca,'FontSize',12)
view([-30 45])
%axis equal
box on
grid on
