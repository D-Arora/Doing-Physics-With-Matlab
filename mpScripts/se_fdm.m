% se_fdm.m
% Solutions to [1D]time independent Schrodinger Equation
% Finite Difference Method
% Enter E - shooting method
% Ian Cooper
% School of Physics, University if Sydney
% email: ian.cooper@sydney.edu.au


close all
clear 
clc


% INPUTS ==========================================
%Ev = input('Enter the energy E in eV    ');
Ev =  -295;           %            -373.7954 -296.63 -173.5 -21;
Uv = -400;                   % depth of potential well [eV]
L = 1e-10;                   % depth of potential well [m] 
N = 1001;                    % number of x values

% CONSTANTS =======================================
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg

disp('Enter a value of E so that psi(end) = 0   default value = -295');
disp('To end the App enter a positive number eg 9)');
disp('   ');

while (Ev < 0)

Ev = input('E  =  ');

 if Ev > 0, break, end

% SETUP CALCULATIONS ==============================
% x coodinates
x_min = -1*L;
x_max = -x_min;
x = linspace(x_min,x_max,N);
dx = x(2)-x(1);

U0 = e * Uv;              % potential well depth [J]        
E = e* Ev;                % total energy [J]
U = zeros(1,N);

psi = zeros(1,N);         % initial conditions
psi(2) = 1;

% FDM =============================================
for n = 2:N-1
   if abs(x(n)) < L/2, U(n) = U0; end
   SEconst = (2*me/hbar^2).*(E - U(n)).* dx^2;
   psi(n+1) = (2 - SEconst) * psi(n) - psi(n-1);
end

% normalize wavefunction
A = simpson1d(psi.*psi,x_min,x_max);
psi = psi ./sqrt(A);

% No. of crossings
cross = 0;
for c = 2 : N
    if psi(c-1)*psi(c) < 0, cross = cross + 1; end
end

fprintf('No. of crossing for psi = %d \n',round(cross)) 
fprintf('End value of psi  =  %0.3g \n',psi(end))
disp('  ');
disp('  ');

% GRAPHICS ========================================
close 
figure(1)
set(gcf,'color',[1 1 1]);
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.3 0.2 0.6 0.4]) 
set(gca,'fontsize',12);
subplot(1,2,1)
plot(x*1e9,U/e,'r','lineWidth',3)
hold on
plot([x(1)*1e9 x(end)*1e9],[Ev Ev],'b');
xlabel('position  x  (nm)')
ylabel('potential energy U (eV)');
set(gca,'fontsize',12)

subplot(1,2,2)
set(gcf,'color',[1 1 1]);
plot(x,psi,'lineWidth',3)
hold on
plot([x_min x_max],[0 0],'k');
plot([x_min/2 x_min/2],[-1.2*max(psi) max(psi)],'r','lineWidth',2);
plot([x_max/2 x_max/2],[-1.2*max(psi) max(psi)],'r','lineWidth',2);
plot([x_min x_min/2],[max(psi) max(psi)],'r','lineWidth',2);
plot([x_max/2 x_max],[max(psi) max(psi)],'r','lineWidth',2);
plot([x_min/2 x_max/2],[-1.2*max(psi) -1.2*max(psi)],'r','lineWidth',2);
axis off




end