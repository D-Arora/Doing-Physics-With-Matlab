% qp_h_pot.m

clear all 
close all
clc

% Inputs ------------------------------------------------------------------
L = 3;              % Orbital quantum number

num = 5001;         % Number of data points

r_min = 1e-12;
r_max = 25e-10;

% Constants ---------------------------------------------------------------
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
eps0 = 8.854e-12;      % F/m

% Setup -------------------------------------------------------------------
% radial position
r = linspace(r_min,r_max,num);

% Calculate potnetial energies
U_c = -(e/(4*pi*eps0))./r;             % Coulomb interaction          

U_L = (hbar^2*L*(L+1)/(2*me*e))./r.^2; % Angular momentum

U_eff = U_c + U_L;                     % Effective potential energy

% Graphics ----------------------------------------------------------------
% black & white plot
figure(1)
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.2 0.15 0.2 0.2]) 
set(gca,'fontsize',8);
numS = [200:1000:num];
plot(r(numS),U_c(numS),'ko')
hold on
plot(r(numS),U_L(numS),'ks')
plot(r,U_eff,'k','LineWidth',2)
h_L = legend('\itU_c', '\itU_L', '\itU_{eff}','Orientation','horizontal');
set(h_L,'Box','off');
grid on
plot(r,U_c,'k')
plot(r,U_L,'k')
axis([r_min r_max -20 20])

tm1 = '{\itL} = '
tm2 = num2str(L,1);
tm = [tm1 tm2];
title(tm,'FontSize',12);
xlabel('radial position  {\itr}  (m)','FontSize',12);
ylabel('potential energy  {\itU}  (eV)','FontSize',12')

% colour plot
figure(2)
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.5 0.15 0.2 0.2]) 
set(gca,'fontsize',8);

plot(r,U_c,'b')
hold on
plot(r,U_L,'r')
plot(r,U_eff,'k','LineWidth',2)
axis([r_min r_max -20 20])
h_L = legend('\itU_c', '\itU_L', '\itU_{eff}','Orientation','horizontal');
set(h_L,'Box','off');
tm1 = '{\itL} = ';
tm2 = num2str(L,1);
tm = [tm1 tm2];
title(tm,'FontSize',12);
xlabel('radial position  {\itr}  (m)','FontSize',12);
ylabel('potential energy  {\itU}  (eV)','FontSize',12')
grid on
