% mathConstants.m

% Ian Cooper
% School of Physics, Sydney University, ian.cooper@sydney.edu.au
% ../mphome.htm
% 19 sept 2017


%% PHYSICAL & USEFUL CONSTANTS [S.I. units] ==============================

c = 299792458;	                % Speed of light [m/s]
e = 1.6021766208e-19;           % Elementary charge [C]
h = 6.626070040e-34;            % Planck constant [ J.s]
hbar = 1.054571800e-34;         % Reduced Planck's constant
Rydberg = 1.097e7;             % Rydbery constant [1/m]

me = 9.10938356e-31;            % Electron mass [kg]
mp = 1.672621898e-27;           % Proton mass [kg]
mn = 1.674927471e-27;           % Neutron mass [kg]
amu = 1.660539040e-27;          % atomic mass unit [kg]

me_u = me/amu;                  % electron mass [amu]
mp_u = mp/amu;                  % proton mass [amu]
mn_u = mn/amu;                  % neutron mass [amu]
mAlpha = 4.0015061791276;         % alpha particle [amu]

NA = 6.022140857e23;            % Avogadro's constant
kB = 1.38064852e-23;            % Boltzmann constant [J/K]
I0 = 1362 ;                     % Solar constant  [W.m^-2]

eps0 = 8.854187817e-12;         % Permitivitty of free spac [C2/N.m2]
mu0 = 4*pi*e-7;                 % Permeability of free space [T.m/A]
ke = 8.987551e9;                 % Coulomb constant 1/4piepso0 [N.m2.C-2]


Rgas = 8.3144598;               % Universal gas constant [J/mol.K]
sigma = 5.670367e-8;            % Stefan-Boltzmann constant [W/m2.K4]
b = 2.8977729e-3;               % Wien Displacement constant m/K] 
Patm = 101325;                  % Standard atmospheric pressure  [Pa]

g = 9.81;                       % Acceleration due to gravity [m/s2]
G = 6.67408e-11;	            % Gravitational constant [N.m2/kg2]

ME = 5.97e24;                   % Earth's mass   [kg]
RE = 6.38e6;                    % Earth's mean radius   [m]
MM = 7.35e22;                   % Moon's mass   [kg]
RM = 1.74e6;                    % Moon's mean radius   [m]
MS = 1.998e30;                  % Sun's mass   [kg]
RS = 6.96e8;                    % Sun's mean radius   [m]
RES = 1.496e11;                 % Earth-Sun distance (mean) [m]
REM = 3.846e8;                  % Earth-Moon distance (mean) [m]



%% CALCULTIONS ===========================================================

Ratio = G*ME*(27.3*24*3600)^2/(4*pi^2*REM^3)
aG = G*ME/REM^2

aC = 4*pi^2*REM/(27.3*24*3600)^2



%% Wien Displacement law

%T = 2000;
% wLpeak = 450e-9;
% T = b/wLpeak
% wLpeak = b/T

%% Rydberg equation
% nI=7
% nF=2
% wL = 1/(Rydberg*(1/nF^2 - 1/nI^2))

%% proton proton cycle
format long
clc
MBefore = 4 * mp_u
MAfter = mAlpha + 2*me_u

dm = MBefore - MAfter

