% mpNPreactions.m

% Calculation of the disintegration energy Q and mass defect dM for
%   nuclear reactions  Reactants --> Products
% Masses are given in amu
% Mass: me electron me / mp proton (hydrogen nucleus) / mn neutron

% Nuclear mass: atomic number Z / mass number A
%  calls the function  mass(Z,A) / function returns mass of nucleus
%  Function mass: stored mass values for isoptoes (Z,A)
%    A new values for isotope mass easily added to the data base
%    Mass values can be found at
%     https://wwwndc.jaea.go.jp/NuC/


% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/mp001.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
%   Download the Script(function) mass.m

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 1906708   Matlab 2018b


clc
close all
clear

% MASSES ==============================================================
% Electron mass [amu] 
  me = 5.485799093287202e-04;
% Proton mass [amu]  mass of hydrogen nucleus Z = 1, A = 1
  mp = 1.007276467284985; 
% Neutron mass [amu] 
  mn = 1.008664915821551; 
 

% INPUTS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Initial state: mass of reactants = [mass(Z,A) ... ]
% Final state:   mass of products  = [mass(Z,A ... ] 
% Example: 4He2 + 14N7 --> 17O8 + 1H1
%          Reactants = [mass(2,4), mass(7,14)] ; 
%          Products =  [mass(8,17),mass(1,1) ];
  Reactants = [ mass(90,232)  ] ; 

  Products =  [ mass(88,228) mass(2,4)  ];

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% CALCULATIONS ========================================================
% Sum: mass of reactants
  mReactants = sum(Reactants);
% Sum: mass of products
  mProducts = sum(Products);
  
% Conversion factor   u --> MeV    u --> kg
  u_MeV = 931.4940954;
  u = 1.66053904e-27;
     
% Mass defect  [u]
  dM = mReactants - mProducts;

% Q-value for nuclear reaction  [MeV]
%   Q > 0  exoenergetic reaction: energy released
%   Q < 0  endoenergetic reaction: energy absorbed
  Q = dM * u_MeV;

% OUTPUT ==============================================================  
disp('  ')
disp('Mass: Reactants')
fprintf('     %3.6f  u  \n',Reactants);
fprintf('     %3.6f  u  \n',mReactants);
disp('  ')
disp('Mass: Products')
fprintf('     %3.6f  u  \n',Products);
fprintf('     %3.6f  u  \n',mProducts);
disp('  ')
disp('Mass defect')
fprintf('     dM = %3.6f  u  \n',dM);
disp('Disintegration value Q')
fprintf('     Q = %3.6f  MeV  \n',Q);


  