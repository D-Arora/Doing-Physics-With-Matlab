% mpNPCB.m

% Calculation: NUCLEAR FUSION AND THE COULOMB BARRIER
%   Two nuclei collide head-on and fuse 
%   Atomic number Z / mass number A
%   Nuclear radius R
%   Average KE of nuclie  Kavg
%   Height of Coulomb Barrier  UC [MeV]
%   Required temperature for fusion  T  [K]

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/mp001.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190709   Matlab 2018b

clc
close all
clear



% INPUT: Z and A for the two nuclei  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Atomic numbers Z = [Z(1), Z(2)]
  Z = [1,1];
% Mass numbers A = [A(1),A(2)] 
  A = [2,3];

  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% DATA ================================================================

% Fundamental charge [C]
    e = 1.60217662e-19;
  % Boltzmann constant
    kB = 1.38064852e-23;
  % Permittivity of free space
    eps0 = 8.85418782e-12;
 

% CALCULATIONS ========================================================    

% Nuclear radii R(1) and R(2)  [m]
  R = 1.2e-15 .* A.^(1/3);
% Couloumb Barrier Height  [MeV]
  UC = e^2*Z(1)*Z(2) / (e*1e6*4*pi*eps0*(R(1)+R(2))); 
% Kinetic energy of incident nuclei
  Kavg = UC/2;
% Temperature 
  T = 2*(Kavg*1e6*e)/(3*kB);

  
% OUTPUT ==============================================================  
  
  disp('  ')
  fprintf('Nuclear radii R = %3.2f %3.2f  fm \n',1e15.*R);
  disp('  ')
  fprintf('Coulomb Barrier height UC = %3.2f  MeV  \n',UC);
  disp('  ')
  fprintf('Average kinetic energy per nuclei  Kavg = %3.2f  MeV  \n',Kavg);
  disp('  ')
  fprintf('Temperature  T = %3.2e  K  \n',T);
