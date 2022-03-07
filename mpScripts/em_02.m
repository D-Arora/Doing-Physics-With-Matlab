% em_02.m

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



%%  CELL 1  
%   Convert Cartesian coordinates of a vector into its
%   cylindrical and spherical coordinates


close all
clear
clc

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Vector Cartesian components
  Ax = 3;
  Ay = 4;
  Az = 6;

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


%%  CELL 2  
%   Convert cylindrical coordinates of a vector into its
%   Cartesian and spherical coordinates

close all
clear
clc

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Vector Cylindrical components
  Arho = sqrt(2);
  phi = pi/4;   %  [rad]
  Az = 1;

% CALCULATIONS ---------------------------------------------------
% Vector A  
  A = [Arho;phi;Az];
% Magnitude Amag
  Ar = sqrt(Arho^2 + Az^2);
% X component
  Ax = Arho*cos(phi);
% Y component
  Ay = Arho*sin(phi) ;
% Azimuthal angle  theta [rad]
  theta = acos(Az/Ar);

% OUTPUT RESULTS TO COMMAND WINDOW
 disp('  ')
 fprintf('Arho = %2.4f   phi = %2.4f  rad   Az = %2.4f  \n',Arho,phi,Az)
 disp('  ')
 fprintf('Ar = %2.4f   Ax = %2.4f   Ay = %2.4f   Az = %2.4f \n',Ar,Ax,Ay, Az)
 disp('  ')
 fprintf('phi = %2.4f  rad   theta = %2.4f rad  \n',phi, theta)
 disp('  ')


%%  CELL 3  
%   Convert spherical coordinates of a vector into its
%   Cartesian and cylindrical coordinates

close all
clear
clc

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Vector Cylindrical components
  Ar = sqrt(3);
  phi = pi/4;   %  [rad]
  theta = 0.9553;

% CALCULATIONS ---------------------------------------------------
% Vector A  
  A = [Ar;phi;theta];
  % Arho
  Arho = Ar*sin(theta);

% X component
  Ax = Ar*cos(phi)*sin(theta);
% Y component
  Ay = Ar*sin(phi)*sin(theta) ;

% Z component
  Az = Ar*cos(theta);

% OUTPUT RESULTS TO COMMAND WINDOW
 disp('  ')
 fprintf('Ar = %2.4f   phi = %2.4f  rad   theta = %2.4f  rad  \n',Ar,phi,theta)
 disp('  ')
 fprintf('Arho = %2.4f   Ax = %2.4f   Ay = %2.4f   Az = %2.4f \n',Arho,Ax,Ay,Az)
 disp('  ')
 fprintf('phi = %2.4f  rad   theta = %2.4f rad  \n',phi, theta)
 disp('  ')

