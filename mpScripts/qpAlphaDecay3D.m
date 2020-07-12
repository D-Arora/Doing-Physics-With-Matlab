%qpAlphaDecay3D.m

% ALPHA DECAY HALF-LIFE COMPUTATION 
%   using the Finite Difference Method to solve the [1D] time independent
%   Schrodinger Equation for a nuclear potential. 

% VARIABLES:  potential energy U, kinetic energy T
%             transition probability P, radial coordinates x.

% Badr-Eddine Elfatehy Chfira: badreddine.el@um.es
% University of Murcia (UMU, Spain) 

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney


% DOING PHYSICS ONLINE: 
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% http://www.physics.usyd.edu.au/teach_res/mp/doc/qpAlphaDecay.pdf



close all
clear all
clc
%==========================================================================
L = 500;                        %radial x coordinate
N = 999999;                    % number of x values
% CONSTANTS =======================================
e = 1.602e-19;              % C
mA = 6.64465675e-27;            %kg
h = 1.054571818e-34;        %J*s ,Planck Constant
eps0 = 8.854187e-12;        % permittivity of free space
Ese = 1.60877e-13;                       % Energy scaling factor  
Lse = 1e-15;                             % Length scaling factor  
Cse = h.^2./(Lse.^2.*Ese)/(2.*mA);       % h^2/2m (Schrod Eq.Constant)
% Inputs ===========================================
E = input(' E(MeV) =     ');
disp('  ')
A = input(' A =     ');
disp('  ')
Z = input(' Z =     ');
disp('   ')
U0 = input(' U0  =     '); 
disp('   ')
disp('   ')

%Data ==============================================
R0 = 1.2* (4^(1/3) + (A-4)^(1/3)); %classical turning point r1, when U=E

% SETUP CALCULATIONS ==============================
% x coodinates
x_min = 0;
x_max = L;
x = linspace(x_min,x_max,N);
dx = x(2)-x(1);
%POTENTIAL:
for a = 0.5:0.1:0.9                       % Saxon-Woods Potential parameter
U = -U0./(1+exp((x-R0)./a)).*ones(1,N);                   % Potential depth
for l = 0:3            % for different values of the orbital quantum number
for n =2:N
   if x(n) > R0
    U(n) = 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n))+l.*(l+1).*Cse./((x(n)).^2)-U0./(1+exp((x(n)-R0)./a)) ; end; %Coulomb+Orbital+Saxon-Woods Potential
end
%===================== WAVEFUNCTION =======================================
% Initial conditions
psi = zeros(1,N); 
% Boundary conditions
psi(N-1) = 2;
psi(N) = 1;
%================ Finite Difference Method ================================
for n = N-1:-1:2
   SEconst = (E - U(n))./Cse.* dx.^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
end   

psi = psi ./max(psi);

% ======== Half_Life Calculation ==========================================
% second turning point for each value of the orbital angular momentum, l = 0, 1 ,2, 3
R1 = max(double(solve(@(z) l.*(l+1).*Cse./(z.^2) +2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* z)==E))); % U(R1)=E 
vA = sqrt(2.*(E+U0).*Ese./mA);         % velocity of escaped alpha particle
f =  vA./ (2.*R0.*Lse);          % frequency at which particles strike barrier on the RHS of well
Ain  = max(abs(psi(1:10000)));     % amplitude of wavefunction inside nucleues
Aout = max(abs(psi(666667:N)));       % amplitude of wavefunction outside nucleus
P = (Aout./ Ain).^2;               % probability of alpha particle escaping 
Pa = 0.0341 ;                               % Pre-formation probability : For Po (Table I) : http://hal.in2p3.fr/in2p3-00637971/document                   
gamma = Pa.*f.* P;                 % decay constant
half_life = log(2)./ gamma;    % half-life
log_half_life = log(half_life);    
% OUTPUTS ==============================================================
T = table(l,U0,E,A,Z,R0,vA,f,P,gamma,half_life,log_half_life)
end
end
%======================== GRAPHICS ========================================
figure(1) % Potential variation with the change of the orbital quantum number 
for a = 0.5 %Fixed value
for l= 0:4
for n =1:N
   if x(n) > R0
    U(n) = 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n)) +l.*(l+1).*Cse./((x(n)).^2)-U0./(1+exp((x(n)-R0)./a)) ;
   else
    U(n) = -U0./(1+exp((x(n)-R0)./a));   
   end
end
plot(x,U,'DisplayName',[' l = ' num2str(l)]);
hold on
end
end
ppp = plot([x(1) x(end)],[E E],'b','linewidth',1.5,'DisplayName','E_{alpha}');
lgd=legend('show');
set(lgd,'Fontsize',8,'Location','southeast');
pp=plot([x(1) x(end)], [0 0],'k');
set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
grid on;
xlabel(' r  (fm)','fontsize',9)
ylabel(' U (MeV)','fontsize',9);
title(['Effective Potential (  A = ' num2str(A) '  , Z = ' num2str(Z) ', U0 = ' num2str(U0) ')'],'fontsize',12);
set(gca,'fontsize',14);
hold off

%==========================================================================
figure(2) % Potential variation with the change of the Saxon Woods Parameter
for a = 0.5:0.1:0.9                       % Saxon-Woods Potential parameter
for l= 1                                     % fixed orbital quantum number
for n =1:N
   if x(n) > R0
    U(n) = 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n)) +l.*(l+1).*Cse./((x(n)).^2)-U0./(1+exp((x(n)-R0)./a)) ;
   else
    U(n) = -U0./(1+exp((x(n)-R0)./a));   
   end
end
end
plot(x,U,'DisplayName',[' a = ' num2str(a)]);
hold on
end
ppp=plot([x(1) x(end)],[E E],'b','linewidth',1.5,'DisplayName','E_{alpha}');
lgd=legend('show');
set(lgd,'Fontsize',8,'Location','southeast');
pp=plot([x(1) x(end)], [0 0],'k');
set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
grid on;
xlabel(' r (fm)','fontsize',9)
ylabel(' U (MeV)','fontsize',9);
title(['Effective Potential (  A = ' num2str(A) '  , Z = ' num2str(Z) ', U0 = ' num2str(U0) ')'],'fontsize',12);
set(gca,'fontsize',14);
hold off

%==========================================================================
figure(3) % Wavefunction variation with the change of the Saxon Woods Parameter
psi = zeros(1,N); 
psi(N-1) = 2;
psi(N) = 1;
for a = 0.5:0.1:0.9                       % Saxon-Woods Potential parameter
for l= 1                                     % fixed orbital quantum number
for n =1:N
   if x(n) > R0
    U(n) = 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n)) +l.*(l+1).*Cse./((x(n)).^2)-U0./(1+exp((x(n)-R0)./a)) ;
   else
    U(n) = -U0./(1+exp((x(n)-R0)./a));   
   end
end
for n = N-1:-1:2
   SEconst = (E - U(n))./Cse.* dx.^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
end   
Area = simpson1d(psi.*psi,x_min,x_max);
psi = psi ./sqrt(Area);
end
plot(x,psi,'DisplayName',[' a = ' num2str(a)]);
hold on
end
lgd=legend('show');
set(lgd,'Fontsize',8);
pp=plot([x(1) x(end)], [0 0],'k');
set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
grid on;
xlabel(' r  (fm)','fontsize',9)
ylabel(' \psi   (fm)','fontsize',9);
title([' Wavefunction (  A = ' num2str(A) '  , Z = ' num2str(Z) ', U0 = ' num2str(U0) ')'],'fontsize',12);
set(gca,'fontsize',14);
hold off

%==========================================================================
figure(4) % Wavefunction variation with the change of orbital quantum number
psi = zeros(1,N); 
psi(N-1) = 2;
psi(N) = 1;
for a = 0.5                     
for l= 0:4                                    
for n =1:N
   if x(n) > R0
    U(n) = 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n)) +l.*(l+1).*Cse./((x(n)).^2)-U0./(1+exp((x(n)-R0)./a)) ;
   else
    U(n) = -U0./(1+exp((x(n)-R0)./a));   
   end
end
for n = N-1:-1:2
   SEconst = (E - U(n))./Cse.* dx.^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
end   
Area = simpson1d(psi.*psi,x_min,x_max);
psi = psi ./sqrt(Area);
plot(x,psi,'DisplayName',[' l = ' num2str(l)]);
hold on
end
end
lgd=legend('show');
set(lgd,'Fontsize',8);
pp=plot([x(1) x(end)], [0 0],'k');
set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
grid on;
xlabel(' r  (fm)','fontsize',9)
ylabel(' \psi   (fm)','fontsize',9);
title([' Wavefunction (  A = ' num2str(A) '  , Z = ' num2str(Z) ', U0 = ' num2str(U0) ')'],'fontsize',12);
set(gca,'fontsize',14);
hold off

%==========================================================================
figure(5)                      %  Closer look to the Wavefunction behaviour 
Mf = 1e20;
set(gcf,'color',[1 1 1]);
index = min(find(x>R1));      % index for x position when E = U outside well
index1 = min(find(x>R0));
plot(x,Mf.*psi,'DisplayName','Wavefunction');
hold on;
pp=plot([x(1) x(end)], [0 0],'k');
set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
lgd=legend('show');
set(lgd,'Fontsize',8);
set(gca,'Ylim',[-1.0*Mf*max(psi(index:N)) 1.0*Mf*max(psi(index:N))]);
set(gca,'fontsize',14);
xlabel('radial coordinate   \it{r}  (fm)','fontsize',14)
ylabel('\psi ','fontsize',14);
title('Wavefunction  x10^{20}','fontsize',14);
grid on;