% se_fdm_barrier.m
% Solutions to [1D]time independent Schrodinger Equation
% Finite Difference Method - Potential Barrier - Tunneling
% Ian Cooper
% School of Physics, University if Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

close all
clear 
clc

% INPUTS ==========================================
disp('Height of potential barrier =  100 eV')
disp('   ')
disp('Enter a value of total energy E (E < 100 eV');
disp('   ');
Ev = input('E  =  ');

Uv = 100;                    % height of potential barrier [eV]
L = 0.2;                     % width of potential barrier [m] 
N = 1001;                    % number of x values  (must be odd)

% CONSTANTS =======================================
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
 
Ese = 1.6e-19;                       % Energy scaling factor  
Lse = 1e-9;                          % Length scaling factor
Cse = -hbar^2/(2*me) / (Lse^2*Ese);  % Schrodinger Eq constant       
 
% SETUP CALCULATIONS ==============================
% x coodinates
x_min = -5*L;
x_max = -x_min;
x = linspace(x_min,x_max,N);
dx = x(2)-x(1);

U0 = Uv;                % potential well depth [J]        
E = Ev;                 % total energy [J]
U = zeros(1,N);

psi = zeros(1,N);       % initial conditions
psi(N-1) = 1;


% FDM =============================================
for n = N-1: -1: 2
   if abs(x(n)) < L/2, U(n) = U0; end
   SEconst = Lse^2 .* Ese .* (2*me/hbar^2).*(E - U(n)).* dx^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
end

% Normalize wavefunction
A = simpson1d(psi.*psi,x_min,x_max);
psi = psi ./sqrt(A);

% GRAPHICS ========================================
close 
figure(1)
set(gcf,'color',[1 1 1]);
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.3 0.2 0.6 0.4]) 
set(gca,'fontsize',14);
subplot(1,2,1)
plot(x,U,'r','lineWidth',3)
hold on
plot([x(1) x(end)],[Ev Ev],'b','linewidth',2);
xlabel('position  x  (nm)','fontsize',14)
ylabel('potential energy U (eV)','fontsize',14);
set(gca,'fontsize',14);

subplot(1,2,2)
set(gcf,'color',[1 1 1]);
plot(x,psi,'lineWidth',3)
hold on
plot([x_min x_max],[0 0],'k');
col = [1 0 0];
plot([-L/2 -L/2],[-max(psi) max(psi)],'color',col,'lineWidth',0.5);
plot([L/2 L/2],[-max(psi) max(psi)],'color',col,'lineWidth',0.5);
set(gca,'Ylim',[-1.2*max(psi) 1.2*max(psi)]);
axis off

