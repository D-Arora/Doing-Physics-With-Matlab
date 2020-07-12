% alpha2.m

% Calls:  pot.m (potential well)
% Calls:  funct.m (Schrodinger equation for ODE)

close all
clear all
clc

global Cse Cpot U1 x1

% Input parameters -------------------------------------------------------
E = 4e6;         % Energy of alpha particle (eV)
xMin = 0;           % Range for X-axis (m)
xMax = 2.5e-13;

Z = 84;             % atomic number of nucleus

% Potential Well inputs
x1 = 8.0986e-15;         % well - width of nucleus (m)
U1 = -40e6;           % well - depth (eV)%global k

% Constants --------------------------------------------------------------
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
eps0 = 8.854e-12;      % F/m
m = 6.648e-27;         % mass of alpha particle (kg)

% Setup ------------------------------------------------------------------
Cse = -2*m*e/hbar^2;            % constant for SE: funct.m
Cpot = 2*(Z-2)*e/(4*pi*eps0);   % constant for potential energy: pot.m

k2 = (2*m*E*e/hbar^2);               % square wave number
k = sqrt(k2);                       % wave number (m^-1)
x_c = 2*(Z-2)*e^2/(4*pi*eps0*E*e);  % Point where E = U
v = hbar*k/m                  % velocity of free alpha particle of energy E
% Solve Schrodinger Equation using ode -----------------------------------
tspan = [xMax xMin];           % range of x values: integrate backward
y0 = [ 1; -i*k];               % initial values for wave function and slope

options = odeset('reltol',1e-6);      
[x y] = ode45('funct',tspan,y0,options,E);

psi = real(y(:,1));           % wave function
prob_density = psi.^2;        % probability density (un-normalized)

pd_max = max(prob_density);    % max probabilty density function

figure(1)                     % plot un-normalized probability density
plot(x,prob_density)
grid on
xlabel('distance  x  (m)');
ylabel('prob density  (m^{-1})');
title('un-normalized wave function');

% Calculate the potential energy from x data from ode --------------------
for cn = 1 : length(x)
 U(cn) = pot(x(cn));
end  % for cn

% Normalize wavefunction -------------------------------------------------
area = abs(trapz(x,prob_density));
prob_density = prob_density/area;

% Calculate half life from flux ------------------------------------------
index = min(find(x<x_c))      % point where E > U or K > 0

esc_prob_density = prob_density(index)   

flux = v * prob_density(index)    % flux of escaping alpha beam
decay_rate = flux
t_half_s = log(2)/decay_rate      % half life (s)
t_half_y = t_half_s/(3600*24*365) % half life  (y)

% Calculate half life from max prob density (un-normalized)
particle_density = x1/pd_max         % particles / m
decay_rate1 = v * particle_density   % decay rate  (s-1)
t_half_s1 = log(2)/decay_rate        % half life (s)
t_half_y1 = t_half_s/(3600*24*365)   % half life  (y)

% Graphics ---------------------------------------------------------------
figure(2)
plot(x,prob_density)
grid on
xlabel('distance  x  (m)');
ylabel('prob density  (m^{-1})');
title('normalized wave function');
xlim([1e-13 2e-13])

figure(3)
plot(x,U)
hold on
grid on
plot([x(1) x(end)], [E E],'r');
xlabel('distance  x  (m)');
ylabel('potential energy  U  (eV)');
title('Potential Well');


function ydot = funct(x,y,flag,E);

global Cse
%k = pi;

%ydot = [y(2); -k^2*(y(1)+ pot(x))];

ydot = [y(2); Cse*(E-pot(x))*y(1)];

%ydot = [y(2); Cse*(E*y(1))];

end

function U = pot(x);

global Cpot U1 x1
U = U1;

if x > x1
    U = Cpot/x;
end

end

