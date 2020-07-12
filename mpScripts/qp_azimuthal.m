% qp_azimuthal.m
% Ian Cooper, School of Physics, The University of Sydney
% Plot for the solution to the Azimuthal differential equation

clear all; close all; clc;
ml = 1;      % change the value of ml to the required value
phi = linspace(0,1,500).* (2*pi);
PHI = exp(j .* ml .* phi);

figure(1)
set(gcf,'Units','Normalized') 
set(gcf,'Position',[0.2 0.15 0.2 0.2]) 
set(gca,'fontsize',8);
x = phi./(2*pi); y1 = real(PHI); y2 = imag(PHI);
plot(x,y1,'linewidth',2);
hold on
plot(x,y2,'r','linewidth',2)
grid on
xlabel('azimuthal angle    \phi / 2\pi')
ylabel('azimuthal   wavefunction     \Phi')
title_m = ['m_l  =  ', num2str(ml)];
title(title_m);