% mpNPBE.m

% Calculation of the binding energy EB and EB/A  for a nucleus
%   Nucleus  --> protons + neutrons
% Masses are given in amu
% Mass: me electron me / mp proton (hydrogen nucleus Z = 1  A = 1) / mn neutron

% Nuclear mass: atomic number Z / mass number A
%  calls the function  mass(Z,A) / function returns mass of nucleus
%  Function mass.m: stored mass values for isopotes (Z,A)
%    A new value for isotope mass can easily be added to the data base
%    Mass values can be found at
%     https://wwwndc.jaea.go.jp/NuC/


% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod72/mod72K1.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
%   Download the Script(function) mass.m

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190713   Matlab 2018b


clc
close all
clear


% INPUT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  % Input atomic number Z and mass number A of isotope of element
  
  Z = 2;  
  A = 4;
  
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% MASS ==============================================================

% Electron mass [amu] 
  me = 5.485799093287202e-04;
% Proton mass [amu]  mass of hydrogen nucleus Z = 1, A = 1
  mp = 1.007276467284985; 
% Neutron mass [amu] 
  mn = 1.008664915821551; 
 
% Conversion factor   u --> MeV    u --> kg
  u_MeV = 931.4940954;
  u = 1.66053904e-27;

  
% CALCULATION =======================================================

% Number of protons
  nP = Z;
% Number of neutrons
  nN = A-Z;
% Mass of nucleus of isotope  [u]
   mNuc = mass(Z,A);
% Mass of nucleons (protons + neutrons)  [u]
   mPN = nP*mp + nN*mn;
% Mass Defect  [u]
  dm = mNuc - mPN;
% Q-value  [MeV]   (amu --> MeV)
  Q = u_MeV * dm;
% Binding Energy [MeV]
  EB = -Q;
% Binding Energy / nucleon   [MeV / nucleon]
  EB_A = EB/A;

  
% OUTPUT ============================================================
  
  fprintf('Atomic number Z = %3.0f    \n',Z);
  fprintf('Mass number A = %3.0f    \n',A);
  fprintf('Number of neutrons A - Z = %3.0f    \n',A - Z);
  fprintf('Mass proton mP = %3.8f    \n',mp);
  fprintf('Mass neutronn mN = %3.8f    \n',mn);
  fprintf('Mass nucleus mNUC = %3.8f    \n',mNuc);
  fprintf('Mass defect  dm = %3.6f  u  \n',dm);
  fprintf('Q-value  Q = %3.4f  u  \n',Q);
  fprintf('Binding energy  EB = %3.4f  MeV  \n',EB);
  fprintf('Binding energy / nucleon   EB/A = %3.4f  MeV  \n',EB_A);
