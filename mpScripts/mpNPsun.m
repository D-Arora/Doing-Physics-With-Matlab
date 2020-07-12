% mpNRsun.m

% Calculation ENERGY OUTPUT OF SUN
%  Total energy from core of Sun from proton-proton cycle
%  Energy output of Sun per year
%  Lifetime of the Sun

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


% DATA ================================================================

% Electron charge [C]
    e = 1.6021766208e-19;
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

   
% CALCULATIONS ========================================================   

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

 
% OUTPUT ==============================================================    

  disp('  ')
  fprintf('Energy output the from core of the Sun  Q_core = %3.2e  J \n',Q_core);
  disp('  ')
  fprintf('Energy output of Sun in one year  E_Sun = %3.2e  J/y  \n',E_sun);
  disp('  ')
  fprintf('Lifetime of Sun = %3.1e  years  \n',LifeTime);
  
    