% mpNuclear.m

% Script for doing nuclear physics calculations:
%   Cell #1: Binding energies and binding energes per nucleon
%   Cell #2: Q-values for nuclear reactions
%   Cell #3: Temperature required for nuclear fusion: Coulomb barrier
%   Cell #4: Lifetime of our Sun

% Function constants: stored numerical values for physics constants

% Function mass: stored values of the isotope mass of atoms
%   New values for the mass of an isotope can be easily added to the data
%   base. Mass values can be found at
%         https://wwwndc.jaea.go.jp/NuC/

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS ONLINE: 
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
% http://www.physics.usyd.edu.au/teach_res/mp/doc/mp001.htm

% 190507   Matlab 2018b


%% Cell #1   **********************************************************
% Mass Defect [u] / Binding energy {MeV]
% Masses are for nuclii in amu [u] (does not include mass of electrons)
% Energy in MeV

global e kB eps0 c u_MeV u me mp mn

  format short
  clc
  close all
  constants
 
% INPUTS ==============================================================
  % Input atomic number Z and mass number Z of isotope of element
  Z = 2;  
  A = 4;
% =====================================================================  

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

  fprintf('Atomic number Z = %3.0f    \n',Z);
  disp('  ')
  fprintf('Mass number A = %3.0f    \n',A);
  disp('  ')
  fprintf('Number of neutrons A - Z = %3.0f    \n',A - Z);
  disp('  ')
  fprintf('Mass defect  dm = %3.4f  u  \n',dm);
  disp('  ')
  fprintf('Q-value  Q = %3.4f  u  \n',Q);
  disp('  ')
  fprintf('Binding energy  EB = %3.4f  MeV  \n',EB);
  disp('  ')
  fprintf('Binding energy / nucleon   EB/A = %3.4f  MeV  \n',EB_A);

  
%% CELL #2   **********************************************************
% NUCLEAR REACTIONS: MASS DEFECT and Q-values 
%  M(reactants) * c^2 ---> M(products) * c^2 + Q

global e kB eps0 c u_MeV u me mp mn

format short
  clc
  close all
  constants

% INPUT  ===============================================================  
% Reactants (Initial state) 
  mReactants = mass(2,4) + mass(7,14) ;
% Products (Final State)
  mProducts = mass(8,17) + mass(1,1) ;
% ======================================================================

% Mass defect  [u]
  dM = mReactants - mProducts;

% Q-value for nuclear reaction  [MeV]
%   Q > 0  exoenergetic reaction: energy released
%   Q < 0  endoenergetic reaction: energy absorbed
  Q = dM * u_MeV;

  fprintf('mass defect  dM = %3.6f  u  \n',dM);
  disp('  ')
  fprintf('Q-value       Q = %3.6f  MeV  \n',Q);
  
 
%% CELL #3   **********************************************************
%  NUCLEAR FUSION and the COULOMB BARRIER 
%  Assume a head-on collsion between incident nuclei
  clear
  close all
  clc
  global e kB eps0 c u_MeV u me mp mn

% INPUT ===========================================================
% Atomic numbers Z(1)and Z(2)
  Z = [1,1];
% Mass numbers A(1)and A(2) 
  A = [1,2];
% =================================================================  
  
% Call function for numerical value of constants
  constants
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

  

%% CELL 4   ***********************************************************
%        ENERGY OUTPUT OF SUN 
clc

% Mass of proton [kg]
    mP = 1.672623e-27;
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
    N_core = M_core / mP;
% Available energy from core by fusion of 4 protons to give helium [J]
    Q_core = eta * (N_core/4) * Q_pp;
% Energy output from Sun in one year [J]
    E_sun = P_sun *(3600*24*365);
% Power output of sun from Solar Constant [W]
    P_solar = S0 * (4*pi*R_SE^2);
% Energy output of Sun from Solar Constant S_sun [J] 
    E_solar = S0 * (4*pi*R_SE^2) * (3600*24*365);
% Lifetime of Sun [years]
    LifeTime = Q_core / E_sun;

  disp('  ')
    fprintf('Energy output the core of the Sun  Q_core = %3.2e  J \n',Q_core);
  disp('  ')
  fprintf('Energy output of Sun in one year  E_sun = %3.2e  J/y  \n',E_sun);
  disp('  ')
  fprintf('Lifetime of Sun = %3.1e  years  \n',LifeTime);
  
    
    
%%   ******************************************************************
function constants

global e kB eps0 c u_MeV u me mn mp

  % fundamental charge
    e = 1.60217662e-19;
  % Boltzmann constant
    kB = 1.38064852e-23;
  % Permittivity of free space
    eps0 = 8.85418782e-12;
  % Speed of light  [m/s]
    c = 2.99792458e8;
  % Conversion factor   u --> MeV    u --> kg
     u_MeV = 931.4940954;
     u = 1.66053904e-27;
  % Mass of electron
    me = 5.4857990907e-04;
  % Mass of neutron
    mn = 1.00866491588; 
  % Mass of proton
    mp = 1.0072766;
 
end


% %   *******************************************************************
% function m = mass(Z,A)
% % NUCLEAR MASSES (not atom) [amu u] 
% %  Element: atomic number Z  / Isotope: mass number A
%   global me
%   constants
%   
% if Z == 1     % hydrogen  
%    if A == 1; m = 1.00782503223; end
%    if A == 2; m = 2.01410177812; end
%    if A == 3; m = 3.01604927790; end
% end   
% 
% if Z == 2     % heliumn  
%    if A == 3; m = 3.01602932007; end
%    if A == 4; m = 4.00260325413; end
% end
%    
% if Z == 3     % lithium 
%    if A == 6;  m = 6.01512288741;   end 
%    if A == 7;  m = 7.01600342665;   end
%    if A == 8;  m = 8.022486238;     end
%    if A == 9;  m = 9.026790199;     end
%    if A == 11; m = 11.043723754;    end
% end   
% 
% if Z == 4     % beryllium 
%    if A == 7;  m = 7.016928707;   end 
%    if A == 8;  m = 8.005305102;   end
%    if A == 9;  m = 9.012183050;   end
%    if A == 10; m = 10.013534679;  end
%    if A == 11; m = 11.021661078;  end
%    if A == 12; m = 12.026920848 ; end
%    if A == 14; m = 14.042892920;  end
% end
% 
% if Z == 6     % carbon
%    if A == 11; m = 11.011433611;   end
%    if A == 12; m = 12.0000000000;  end
%    if A == 13; m = 13.00335483507; end
%    if A == 14; m = 14.00324198842; end
% end
% 
% 
% if Z == 7     % nitrogenn  
%    if A == 13; m = 13.005738609; end
%    if A == 14; m = 14.00307400443; end
%    if A == 15; m = 15.00010889888; end
% end
%   
% if Z == 8     % oxygen  
%    if A == 15; m = 15.003065618; end  
%    if A == 16; m = 15.99491461957; end  
%    if A == 17; m = 16.99913175650; end  
% end 
% 
% if Z == 25     % iron  
%    if A == 55; m = 54.938043937; end  
%    
% end 
% 
% if Z == 26     % iron  
%    if A == 55; m = 54.938292023; end  
%    if A == 56; m = 55.934936355; end 
% end 
% 
% if Z == 36     % krypton  
%    if A == 83; m = 82.914127143;   end 
%    if A == 84; m = 83.91149772821; end
%    if A == 85; m = 84.912527262;   end
%    if A == 86; m = 85.910610627;   end
%    if A == 87; m = 86.913354760;   end
%    if A == 88; m = 87.914447881;   end
%    if A == 89; m = 88.917835451;   end
%    if A == 90; m = 90.923806311;   end
%    if A == 91; m = 90.923806311;   end
%    if A == 92; m = 91.926173095;   end
%    if A == 93; m = 92.931147175;   end
%    if A == 94; m = 93.934140455;   end
%    if A == 95; m = 94.939710924;   end
%    if A == 96; m = 95.943016618;   end
%    if A == 97; m = 96.949088785;   end
% end 
% 
% if Z == 56     % barium  
%    if A == 138; m = 137.905246899; end 
%    if A == 139 ; m = 138.908841004; end
%    if A == 140; m = 139.910605633; end
%    if A == 141; m = 140.914402957; end
%    if A == 142; m = 141.916429935; end
%    if A == 143; m = 142.920625195; end
%    if A == 144; m = 143.922954824; end
%    if A == 145; m = 144.927518400; end
%    if A == 146; m = 145.930283286; end
%    if A == 147; m = 146.935303900; end
%    if A == 148; m = 147.938170578; end
%    if A == 149; m = 148.942920; end
%    if A == 150; m = 149.945950; end
%    if A == 151; m = 150.951080; end
% end 
% 
% if Z == 86     % radon 
%    if A == 222; m = 222.017577269; end  
%    if A == 226; m = 226.030852862; end  
% end 
% 
% if Z == 88     % radium  
%    if A == 226; m = 226.025409353; end  
% end
% 
% if Z == 90     % thorium  
%    if A == 226; m = 226.024891; end   
%    if A == 228; m = 228.028734823; end  
%    if A == 234; m = 234.043602450; end
% end
% 
% if Z == 91     % protactinium 
%    if A == 231; m = 231.035884277; end 
%    if A == 237; m = 237.051147110; end  
% end
% 
% if Z == 92     % uranium  
%    if A == 217; m = 217.024382652; end  
%    if A == 218; m = 218.023522502; end
%    if A == 219; m = 219.025016432; end
%    if A == 222; m = 222.026081   ; end
%    if A == 223; m = 223.027737909; end
%    if A == 224; m = 224.027603882; end
%    if A == 225; m = 225.029389831; end
%    if A == 226; m = 226.029337817; end
%    if A == 227; m = 227.031155482; end
%    if A == 228; m = 228.031373109; end
%    if A == 229; m = 229.033505084; end
%    if A == 230; m = 230.033938943; end
%    if A == 231; m = 231.036293977; end
%    if A == 232; m = 232.037149836; end
%    if A == 233; m = 233.039636574; end
%    if A == 234; m = 234.040953616; end
%    if A == 235; m = 235.043931368; end
%    if A == 236; m = 236.045569468; end
%    if A == 237; m = 237.048731636; end
%    if A == 238; m = 238.050789466; end
%    if A == 239; m = 239.054294518; end
%    if A == 240; m = 240.056593384; end
%    if A == 242; m = 242.062933   ; end
% end 
%    
%   m = m - Z*me;
% end

  
  

