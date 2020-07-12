% se_infWell.m
% Ian Cooper
% School of Physics, University of Sydney
% Calculates the energy eigenvalues for an
%     electron in an infinite square well by a analytical method
%     and compares results with the Matrix Method 

close all
clear 
clc

% Inputs
n = 1:10;        % principle quantum number
L = 0.1e-9;
% Constants
h =  6.63e-34;     % Planck's constant [J.s]
me = 9.1e-31;      % mass of electron [kg]
e = 1.60e-19;      % electron charge [C]

% Calculations

E = (h^2/(8*me*L^2)) .* n.^2;     % total energy [J]

E = E ./ e;                       % change units J to eV


% Data from se_wells.m and se_solve.m for 
%    infinite square well of width 0.1 nm
Emm(1:6) = [37.6863, 150.7446, 339.1732, 602.9691, 942.1283, 1.3566e+003];
Emm(7:10) = [1.8465e+003, 2.4117e+003, 3.0523e+003, 3.7681e+003];

disp(' n   E (eV)    Emm (eV)')
for c = 1 : 10
output = [num2str(c,'%0.0f') '    ' num2str(E(c),'%0.2f') '   ' num2str(Emm(c),'%0.2f')];
disp(output);

end

figure(1)
set(gcf,'Color',[1 1 1]);
set(gca,'Fontsize',12);


plot(n,E,'o')
hold on
plot(n,Emm,'+r')
xlabel('Quantum No. n','Fontsize',12);
ylabel('Total Energy  E  (eV)','Fontsize',12);
legend('analytical','Matrix Method');



