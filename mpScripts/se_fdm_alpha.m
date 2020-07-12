% se_fdm_alpha.m
% Solutions to [1D]time independent Schrodinger Equation
% Finite Difference Method
% Half-life calculation by solving S.E.
% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm

close all
clear all
clc

% INPUTS =================================================================
disp('Enter data for alpha particle and the parent nucleus:')
disp('  ')
E = input('Enter the alpha particle energy E in MeV     ');
disp('  ')
A = input('Enter the mass number A for the parent nucleus     ');
disp('  ')
Z = input('Enter the atomic number Z for the parent nucleus     ');
disp('   ')

U0 = -40;                   % depth of potential well [MeV]
L = 120;                    % max radial coordinate ([fm] 
N = 59009;                   % number of x values - must be ODD

% CONSTANTS ===============================================================
hbar = 1.055e-34;           % J.s
e = 1.602e-19;              % C
mA = 6.64465675e-27;        % kg
eps0 = 8.854187e-12;        % permittivity of free space


Ese = 1.6e-13;                       % Energy scaling factor  
Lse = 1e-15;                         % Length scaling factor
Cse = -hbar^2/(2*mA) / (Lse^2*Ese);  % Schrodinger Eq constant       
 
% SETUP CALCULATIONS ==============================
% x coodinates
x_min = 0;
x_max = L;
x = linspace(x_min,x_max,N);
dx = x(2)-x(1);
 
U = U0 .* ones(1,N);

psi = zeros(1,N);         % initial conditions

x0 = 1.07 * (4^(1/3) + (A-4)^(1/3));
x1 = 2*(Z-2)*e^2/(4*pi*eps0*E*Ese*Lse);

% FDM =============================================
for n = 2:N
   if x(n) >  x0, U(n) = 9e9 * 2 * (Z-2) * e^2 / (Ese * Lse * x(n)); end;
end

psi(N) = 1; psi(N-1) = 2; 

for n = N-1:-1:2    
   SEconst = Lse^2 .* Ese .* (2*mA/hbar^2).*(E - U(n)).* dx^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
end

% normalize wavefunction
Area = simpson1d(psi.*psi,x_min,x_max);
psi = psi ./sqrt(Area);

index = min(find(x>x1));      % index for x position when E = U outside well
index1 = min(find(x>x0));
vA = sqrt(2*E*Ese/mA);         % velocity of escaped alpha particle
f =  vA / (2*x0*Lse);          % frequency at which particles strike barrier on the RHS of well
Ain  = max(abs(psi(1:index1)));     % amplitude of wavefunction inside nucleues
Aout = max(abs(psi(index:N)));       % amplitude of wavefunction outside nucleus
Aout = max(abs(psi(2 * index:N)));
P = (Aout / Ain)^2;            % probability of alpha particle escaping 
gamma = f * P;                 % decay constant
half_life = log(2) / gamma;    % half-life

% OUTPUTS ==============================================================
E
A
Z
x0
x1
vA
f
Ain
Aout
P
gamma
half_life


% GRAPHICS ==============================================================
close 
figure(1)
set(gcf,'color',[1 1 1]);
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.1 0.1 0.6 0.8]) 
set(gca,'fontsize',14);

subplot(3,1,1)
plot(x,U,'r','lineWidth',3)
hold on
plot([x(1) x(end)],[E E],'b','linewidth',2);
%xlabel('position  x  (fm)','fontsize',14)
ylabel('E, U (MeV)','fontsize',14);
set(gca,'fontsize',14);
plot([x(index) x(index)],[-1.2*max(U) 1.2*max(U)],'color',[0.7 0.7 0.7],'lineWidth',0.5);
plot([0 0],[-40 40],'r','lineWidth',3);

subplot(3,1,2)
set(gcf,'color',[1 1 1]);
plot(x,psi,'lineWidth',3)
hold on
set(gca,'Ylim',[-1.2*max(psi) 1.2*max(psi)]);
plot([x(index) x(index)],[-1.2*max(psi) 1.2*max(psi)],'color',[0.7 0.7 0.7],'lineWidth',0.5);
set(gca,'fontsize',14);
%xlabel('position  x  (fm)','fontsize',14)
ylabel('\psi ','fontsize',14);

subplot(3,1,3)
Mf = 1e20;
set(gcf,'color',[1 1 1]);
plot(x,Mf.*psi,'lineWidth',3)
hold on
plot([x_min x_max],[0 0],'k');
plot([x(index) x(index)],[-1.2*Mf*max(psi(index:N)) 1.2*Mf*max(psi(index:N))],'color',[0.7 0.7 0.7],'lineWidth',0.5);
set(gca,'Ylim',[-1.2*Mf*max(psi(index:N)) 1.2*Mf*max(psi(index:N))]);
set(gca,'fontsize',14);
xlabel('radial coordinate   \it{r}  (fm)','fontsize',14)
ylabel('\psi ','fontsize',14);
title('magnification x10^{20}','fontsize',14);



