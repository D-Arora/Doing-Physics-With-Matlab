% da_constants.m

% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% mscript for the numerical values of commonly used physical constants
% All units are S.I. unless stated otherwise

%
clear all
close all
clc

%%

cL = 2.99792458e8;      % speed of light

h = 6.62606957e-34;     % Planck's constant

hbar = 1.054571726e-34; % hbar Planck's constant

kB = 1.3806488e-23;     % Boltzman constant

me = 9.10938291e-31;    % electron mass

mp = 1.672621777e-27;   % proton mass

mn = 1.674927351e-27;   % neutron mass

mA = 6.64465675e-27;    % alpha particle mass 

e = 1.602176565e-19;    % elementary charge

eps0 = 8.854187e-12;    % permittivity of free space

mu0 = 4*pi*1e-7;        % permeability of free space

kC = 8.987552e9;        % Coulomb constant

G = 6.67384e-11;        % Universal Gravitation constant



%%
L = 4e-10
n = [4 5]
En = n.^2. * h^2/(8*me*L^2*e)
mean(En)
wL = 2*L ./ n
w = e*En./hbar
T = 2*pi./w
Eavh = (9*En(1)+En(2))/10


%% square well
E = e * 150.34
w = (3*pi/2)*hbar / sqrt(2*me* E)
ws = w/4e-12
ws/2

%%
w =  [33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 ... 
    57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83  ];
pT = [47 84 81 43 30 32 51 85 77 41 30 34 55 86 75 40 30 36 58 85 69 38 30 38 ... 
    61 84 65 38 31 40 64 82 62 37 32 42 66 80 59 36 33 45 68 78 57 36 34 47 69 75 55];


figure(99)
plot(w*dx/wL,pT,'linewidth',2)
ylabel('T %','fontsize',14);
xlabel('w / \lambda','fontsize',14);
grid on
set(gca,'fontsize',14);


%%
w = [40:80] * dx;
E = 59.37*e;
k = sqrt(2*me*E)/hbar;
k1 = 600*e; k2 = k1+E;
p1 = sin(k.*w).^2;
ylabel('T %');
xlabel('x  [m]');
z = 1 + 0.25 .* k1^2 .* p1 ./(k2*(k2-k1));
z = 1./z;

figure(98)
plot(w,p1)

figure(97)
plot(w./wL,z)
