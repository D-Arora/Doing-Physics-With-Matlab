% qp_alpha_theory.m

% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Calculation of the theorical value of the half-life of an alpha emitter
% All units are S.I. unless stated otherwise

close all
clear 
clc

% INPUT PARAMETERS =======================================================
A = 212;                % mass number of parent nucleus
Z = 92;                 % atomic number of parent nucleus
EA_MeV = 8.95;           % kinetic energy of alpha particle in MeV

% Constants ==============================================================
e = 1.602176565e-19;    % elementary charge
hbar = 1.054571726e-34; % hbar Planck's constant
mA = 6.64465675e-27;    % alpha particle mass 
eps0 = 8.854187e-12;    % permittivity of free space

% Calculations ===========================================================
EA = EA_MeV *1e6*e;     % kinetic energy ofalpha particle in joules

R0 = 1.07e-15*(4^(1/3) + (A-4)^(1/3));    % Radius of nucleus R0
%R0 = 1.07e-15*((A)^(1/3));    % Radius of nucleus R0

vA = sqrt(2*EA/mA);     % escape velocity of alpha particle

f = vA / (2*R0);        % frequency of collisions between alpha / barrier

K11 = (4*e/hbar)*(mA/(pi*eps0))^0.5;
K1 = K11 * R0^0.5 * (Z-2)^0.5; 

K22 = -(e^2/(hbar * eps0)) * (mA/2)^0.5;
K2 = K22 * (Z-2) * EA^(-0.5);

P = exp(K1 +K2);        % probability

gamma = f * P;          % decay constant

t_half = log(2) / gamma;  % half-life


t_half
f
P
R0

% Test data for Po Z = 84 ================================================
   % Measured values
      hLPo = [183 0.158 1.83e-3 1.6e-4 4.2e-6 3e-7 0.52 1.2e7 3.16e9 9.15e7];
      EPo = [6.12 6.774 7.365 7.68 8.34 8.776 7.59 5.30 4.88 5.11];
   % Computed half-lifes by running the program for each isoptope of Po
      hLPo_c = [2.7e6  5.4574e-08  9.6196e-05 0.0012 0.158 364.7975];
   % Computed data by solving the Schrodinger Equation by FDM   
      hLPo_f = [2.9e6 3.7e-8 4.6e-4 2.4e-5 0.45 602 ];
   % Computed data by solving the Schrodinger Equation by FDM   
      hLPo_f1 = [1.1e7 7.7e-8 8.7e-4 5.3e-5 1.7 2.4e3 ];   
      
% Test data for U Z = 92 ================================================
   % Measured values
      hLU = [1.42e17 7.57e14 2.25e16 7.95e12 5.11e12 2.34e9 1.81e6];
      EU = [4.2 4.5 4.4 4.77 4.82 5.32 5.89];      
   
% Variation of half-life with E(alpha) ===================================
num = 100;
Emin = 4;               % MeV
Emax = 10;              % MeV
E = linspace(Emin,Emax,num);
Ecc = (1e6*e) .* E;     % E in joules
Kcc = -(e^2/(hbar * eps0)) .* (mA/2)^0.5 .* (Z-2) .* Ecc.^(-0.5);

Pcc = exp(K1 +Kcc);        % probability

gammacc = f .* Pcc;          % decay constant

t_halfcc = log(2) ./ gammacc;  % half-life

xp = E.^(-0.5); yp = log10(t_halfcc);

figure(1)
set(gca,'fontsize',12);
plot(xp,yp,'linewidth',2)

hold on
xp = (EPo).^(-0.5); yp = log10(hLPo);
plot(xp,yp,'o','MarkerFaceColor',[0 0 1],'MarkerSize',8,'MarkerEdgeColor',[0 0 1])
% yp = log10(hLPo_c);
% h_plot = plot(xp,yp,'o','MarkerFaceColor',[0 0 1],'MarkerSize',6);
% yp = log10(hLPo_f);
% h_plot = plot(xp,yp,'+','MarkerEdge',[1 0 1],'MarkerSize',9);
% yp = log10(hLPo_f1);
% h_plot = plot(xp,yp,'d','MarkerFaceColor',[0 0 0],'MarkerEdge',[0 0 0],'MarkerSize',8);

xp = (EU).^(-0.5); yp = log10(hLU);
plot(xp,yp,'o','MarkerFaceColor',[1 0 0],'MarkerSize',8,'MarkerEdgeColor',[1 0 0])


xlabel('KE_{alpha} ^{-0.5}','fontsize',14);
ylabel('log_{10}( t_{1/2} )','fontsize',14);
grid on

