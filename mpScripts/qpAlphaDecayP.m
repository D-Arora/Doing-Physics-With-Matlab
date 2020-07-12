% qpAlphaDecayP.m


%%%%%%%%%%%%%%%%% Alpha Log Halflife vs       Energy^-1/2 %%%%%%%%%%%%%%%%%%%%
% Badr-Eddine Elfatehy Chfira: badreddine.el@um.es
% University of Murcia (UMU, Spain) 

% DOING PHYSICS ONLINE: 
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% http://www.physics.usyd.edu.au/teach_res/mp/doc/qpAlphaDecay.pdf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 300;                        %radial x coordinate
N = 999999;                    % number of x values
% CONSTANTS =======================================
e = 1.602e-19;              % C
mA = 6.64465675e-27;            %kg
h = 1.054571818e-34;        %J*s
eps0 = 8.854187e-12;        % permittivity of free space

Ese = 1.60877e-13;                       % Energy scaling factor  
Lse = 1e-15;                             % Length scaling factor  
Cse = h.^2./(Lse.^2.*Ese)/(2.*mA);       % h^2/2m (Schrod Eq.Constant)
% Inputs ===========================================
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Un = zeros(1,N); Uo = Un; Uc = Un; U = Un;

c = 0; %initial value for the count
for E = 4:10
    c=c+1;% counter to store E and loghalflife data for the plot
a = 0.5;
l = 1;
for n =1:N
    Un(n) =  -U0./(1+exp((x(n)-R0)./a));                                   % U_Saxon-Woods
    Uo(n) = 0.0;
    Uc(n) = 0.0;
   if x(n) > R0
       Uc(n) =  2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n));          % U_Coulomb
       Uo(n) = +l.*(l+1).*Cse./((x(n)).^2);                                % U_Centrifugal 
   end
   U(n) = Un(n)+Uc(n)+Uo(n);                                               % U_effective
end
% second turning point for each value of the orbital angular momentum, l = 0, 1 ,2, 3
R1 = max(double(solve(@(z) l.*(l+1).*Cse./(z.^2) +2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* z)==E))); % U(R1)=E 
% ========================FDM =============================================
% Initial conditions
psi = zeros(1,N); 
% Boundary conditions
psi(N-1) = 2;
psi(N) = 1;
for n = N-1:-1:2
   SEconst = (E - U(n))./Cse.* dx.^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
end   
%  wavefunction normalization
Area = simpson1d(psi.*psi,x_min,x_max);
psi = psi ./sqrt(Area);
% ======== Log Half_Life Calculation ======================================
index = min(find(x>R1));      % index for x position when E = U outside well
index1 = min(find(x>R0));     % index for x position when E = U inside well 
vA = sqrt(2.*(E+U0./(1+exp((-R0)./a))).*Ese./mA);         % velocity of escaped alpha particle
f =  vA./ (2.*R0.*Lse);          % frequency at which particles strike barrier on the RHS of well
Ain  = psi(index1);     % amplitude of wavefunction inside nucleues
Aout = psi(index);       % amplitude of wavefunction outside nucleus
P = (Aout./ Ain).^2;            % probability of alpha particle escaping 
gamma = f.* P;                 % decay constant
half_life = log(2)./ gamma;    % half-life
log_half_life = log(half_life);
% DATA ====================================================================
ss(c) = 1./sqrt(E);
ff(c) = log_half_life;
end

% OUTPUT ==================================================================
figure(1)
  set(gcf,'color',[1 1 1]);
  set(gcf,'Units','Normalized') 
  set(gcf,'Position',[0.1 0.1 0.4 0.4]) 
plot(ss,ff,'*','DisplayName',[' Z = ' num2str(Z)]);
hold on;
lgd=legend('show');
set(lgd,'Fontsize',8);
xlabel(' 1/sqrt(Kinetic Energy)  (MeV^-1/2)','fontsize',10)
ylabel('  Log ( t_{1/2} ) ( s )','fontsize',10);
title(['\alpha-Decay : Log(t_{1/2}) vs  Kinetic Energy ^-1/2 (  A = ' num2str(A) '  , Z = ' num2str(Z) ', U0 = ' num2str(U0) ')'],'fontsize',12);
grid on;
% Lets make a curve fit of the graph
p = polyfit(ss,ff,1);  % 2 for a curve 1 for a linear regression, some fits look like a curve other like a straight line.
yy = polyval(p,ss);
plot(ss,yy,'DisplayName','Data fit');
lgd=legend('show');
hold off
slope=p(1);
intercept=p(2);
table(slope,intercept)

set(gca,'fontsize',12)