% se_alphaB1.m

% Solutions to [1D]time independent Schrodinger Equation
% Finite Difference Method: ALPHA DECAY
% Nuclear Potentials Models
% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/

close all
clear 
clc

% INPUTS ==============================================================
% Atomic number
Z = 92;
% Mass number
A = 238;
% Potential well depth: U0 > 0 [MeV]   
U0 = 40;
% Max radial coordinate [fm] 
xMax = 100;
% Number of x gris points (odd number)                 
N = 99009;                   


% CONSTANTS ===========================================================
   h = 6.62607004e-34;         % Planck constant [J.s]
   e = 1.602e-19;              % Fundamanetal charge [C]
   eps0 = 8.854187e-12;        % permittivity of free space [S.I. units]
   mA = 6.64465675e-27;        % mass of alpha particle [kg]
   
   
% SETUP ===============================================================
   Ese = 1e6*e;                         % Energy scaling factor  
   Lse = 1e-15;                         % Length scaling factor
   hbar = h/(2*pi);                     % hbar
   Cse = -hbar^2/(2*mA) / (Lse^2*Ese);  % Schrodinger Eq constant    
   k = 1/(4*pi*eps0);                   % Coulomb constant
   K = k*2*(Z-2)*e^2/(Ese*Lse);         % Cacluation constant
   R0 = 1.26 * (4^(1/3) + (A-4)^(1/3)); % Potential Well width
      
% Radial coodinates [fm]
   x = linspace(0, xMax,N);
   dx = x(2)-x(1);
% Potential [meV] 
   U = -U0 .* ones(1,N); 

   % Angular momentum quantum number L
figure(1)
set(gcf,'color',[1 1 1]);
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.1 0.1 0.3 0.4]) 
set(gca,'fontsize',14);

a = 2;
U =  -U0./(1+exp((x-R0)./a));   
for  L = 0:4  
  for n = 2:N
    if x(n) >  R0
      U(n) = K/x(n) + Cse*L*(L+1)/x(n)^2 + U0 / (1+ exp((x(n)-R0))/a);
    end
  end

  plot(x,U,'lineWidth',2)
  hold on
end

xlabel('position  x  (fm)','fontsize',14)
ylabel('U (MeV)','fontsize',14);
set(gca,'fontsize',14);
box on
grid on



