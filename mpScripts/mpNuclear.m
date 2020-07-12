% mpNuclear.m

% https://wwwndc.jaea.go.jp/NuC/

% All masses are in amu  [u]
% Mass of nuclei m  element / isotope

clc
close all
% clear
format long

% CONSTANTS ===========================================================
% fundamental charge
   e = 1.60217662e-19;
% Boltzmann constant
   kB = 1.38064852e-23;
% Permittivity of free space
   eps0 = 8.85418782e-12;
% Speed of light  [m/s]
  c = 2.99792458e8;
% Conversion factor   u --> MeV    u --> kg
   M = 931.4940954;
   u = 1.66053904;
% Mass of electron  [u]
  me = 5.4857990907e-04;
% Mass of neutron
  mn = 1.00866491588;  
% Mass of proton [kg]
  mp = 1.007276466879;
  
% NUCLEAR MASSES (not atom) amu [u] =================================== 
%  Element: atomic number Z  / Isotope: mass number A
   m = zeros(100,250);
   
Z = 1;     % hydrogen  
   A = 1; m(Z,A) = 1.00782503223 - Z*me;
   A = 2; m(Z,A) = 2.01410177812 - Z*me;
   A = 3; m(Z,A) = 3.01604927790 - Z*me;

Z = 2;     % heliumn  
   A = 3; m(Z,A) = 3.01602932007 - Z*me;
   A = 4; m(Z,A) = 4.00260325413 - Z*me;

Z = 4;     % beryllium  
   A = 8; m(Z,A) = 8.005305102 - Z*me;
    
Z = 6;     % carbonn  
   A = 12; m(Z,A) = 12.0000000000  - Z*me;
   A = 13; m(Z,A) = 13.00335483507 - Z*me;
 
Z = 7;     % nitrogenn  
   A = 13; m(Z,A) = 13.005738609   - Z*me;
   A = 14; m(Z,A) = 14.00307400443 - Z*me;
   A = 15; m(Z,A) = 15.00010889888 - Z*me;
  
Z = 8;     % oxygen  
   A = 15; m(Z,A) =  15.003065618 - Z*me;  
  
  
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% NUCLEAR REACTIONS: reactants ---> products 
% INPUT for Reactants m(Z,A)
  mReactants = m(1,2);
% INPUT for Products  m(Z,A) 
  mProducts = mn+mp ;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Mass defect  [u]
  dM = mProducts - mReactants;

% Q-value for nuclear reaction  [MeV]
%   Q < 0  exothermic reaction: energy released
%   Q > 0  endothermic reaction: energy absorbed
  Q = dM * M;

  fprintf('mass defect  dM = %3.6f  u  \n',dM);
  disp('  ')
  fprintf('Q-value       Q = %3.6f  MeV  \n',Q);
  
 
%% CELL #2
%  COULOMB BARRIER ====================================================
% Assume a head-on collsion between incident nuclei
% INPUTS
% Atomic numbers Z(1)and Z(2)
 Z = [6,1];
% Mass numbers A(1)and A(2) 
 A = [12,1];
 
% Nuclear radii R(1) and R(2)  [m]
  R = 1.2e-15 .* A.^(1/3);
% Couloumb Barrier Height  [MeV]
  UC = e^2*Z(1)*Z(2) / (e*1e6*4*pi*eps0*(R(1)+R(2))); 
% Kinetic energy of incident nuclei
  Kavg = UC/2;
% Temperature 
  T = 2*(Kavg*1e6*e)/(3*kB);

  disp('  ')
  disp('*************************************')
  fprintf('Nuclear radii R = %3.2f %3.2f  fm \n',1e15.*R);
  disp('  ')
  fprintf('Coulomb Barrier height UC = %3.2f  MeV  \n',UC);
  disp('  ')
  fprintf('Average kinetic energy per nuclei  Kavg = %3.2f  MeV  \n',Kavg);
  disp('  ')
  fprintf('Temperature  T = %3.2e  K  \n',T);

  
%% ENERGY OUTPUT OF SUN =============================================

format short
% Mass of proton [kg]
    mp = 1.672623e-27;
% Average density of core  [kg/m^3]
    d_core = 1.5e5;
% Radius of Sun [m]
    R_sun = 6.96e8;
% Distance Sun - Earth [m]
    R_SE = 1.496e11;
% Q-value for proton-proton cycle [J]
    Q_pp = 25*e*1e6;
% Power out of Sun [W]
    P_sun = 3.846e26;
% Solar constant [W.m^2]
   S0 = 1.368e3;
% Proportion of core hot enough to sustain nuclear fusion
   eta = 0.10;
   
% Radius of core [m]
    R_core = 0.25 * R_sun;
% Volume of core [m^3]
    V_core = (4/3)*pi*R_core^3;
% Mass of core [kg]
    M_core = d_core * V_core;
% Number of protons in core
    N_core = M_core / mp;
% Available energy from core by fusion of 4 protons to give helium [J]
    Q_core = eta * (N_core/4) * Q_pp
% Energy output from Sun in one year [J}
    E_sun = P_sun *(3600*24*365)
% Power output of sun from Solar Constant [W]
    P_solar = S0 * (4*pi*R_SE^2)
% Energy output of Sun from Solar Constant S_sun [J] 
    E_solar = S0 * (4*pi*R_SE^2) * (3600*24*365)
% Lifetime of Sun [years]
    LifeTime = Q_core / E_sun
    
 
    
    
  
  

