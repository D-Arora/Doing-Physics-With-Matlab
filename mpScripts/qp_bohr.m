% qp_constants.m

% QUANTUM PHYSICS NUMERICAL VALUES

% Values assigned to the physical constants called by other scripts
% All values are given in SI. units unless stated otherwise


close all
clear all

% Constants
c = 2.99e8;            % speed of light (m/s) 
h = 6.6256e-34;         % Planck's constant (J.s)
hbar = 1.055e-34;      % Planck's constant hbar (J.s)
e = 1.60210e-19;         % electron charge (C)
me = 9.109e-31;        % electron mass (kg)
mp = 1.6726e-27;       % mass of proton (kg)
mn = 1.6749e-27;
eps0 = 8.8543e-12;      % permittivity of free space (F/m)

% Inputs
Z = 3;                 % atomic number

m = mp;               % hydrogen - 1
%  m = 2*(mp + mn);       % helium - 4
% m = 3*(mp + mn);      % lithium - 6

% Calculatins
mr = m * me / (m + me);
mr = me;
E0 = - mr * Z^2 * e^4 / (8 * eps0^2* h^2);
a0 = h^2 * eps0 / (pi * me * e^2)
% Outputs
n = 1: 10;
Rn = (n.^2 .* a0) .* 1e10 ./Z;
En = (E0 ./ n.^2) ./ e;

Rn(1)
En(1)
disp('       n       Rn       En')
disp([n; Rn; En]')









