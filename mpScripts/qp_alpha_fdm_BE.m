% qp_alpha_fdm_BE
% Solutions to [1D]time independent Schrodinger Equation
% Finite Difference Method
% Enter E - shooting method
% Ian Cooper
% School of Physics, University if Sydney
% email: ian.cooper@sydney.edu.au

close all
clear 
clc

L = 120;                        %radial x coordinate
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
R0 = 1.26* (4^(1/3) + (A-4)^(1/3)); %classical turning point r1, when U=E
% SETUP CALCULATIONS ==============================
% x coodinates
x_min = 0;
x_max = L;
x = linspace(x_min,x_max,N);
dx = x(2)-x(1);
        
%================= Effective potential chosen : U_eff= U_Saxon-Woods + U_orbital + U_Coulomb=========================================== 
for a = 0.25:0.25:1                                                          % For different values of the diffusion parameter of Saxon-Woods Potential
% Initial conditions
U = ( -U0./(1+exp(-R0./a))).* ones(1,N);  % U0>=0  UsaxonWoods   Depth
for Ll = 0:2                                                                % For different orbital angular momentum
for n =2:N
    U(n) =  -U0./(1+exp((x(n)-R0)./a));                                   % U_Saxon-Woods
   if x(n) > R0
       U(n) = -U0./(1+exp((x(n)-R0)./a))+ 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n))+Ll.*(Ll+1).*Cse./((x(n)).^2);
   end
end

% second turning point for each value of the orbital angular momentum, l = 0, 1 ,2, 3
R1 = max(double(solve(@(z) Ll.*(Ll+1).*Cse./(z.^2) +2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* z)==E))); % U(R1)=E 
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
% ======== Half_Life Calculation ==========================================
index = find(x>R1, 1 );      % index for x position when E = U outside well
index1 = find(x<R0, 1, 'last' );     % index for x position when E = U inside well 
vA = sqrt(2.*E.*Ese./mA);         % velocity of escaped alpha particle
f =  vA./ (2.*R0.*Lse);          % frequency at which particles strike barrier on the RHS of well
Ain  = psi(index1);     % amplitude of wavefunction inside nucleues
Aout = psi(index);       % amplitude of wavefunction outside nucleus
P = (Aout./ Ain).^2;            % probability of alpha particle escaping 
gamma = f.* P;                 % decay constant
half_life = log(2)./ gamma;    % half-life
% OUTPUTS ==============================================================
T=table(a,Ll,U0,E,A,Z,R0,R1,vA,f,P,gamma,half_life)
% GRAPHICS ==============================================================
figure(1)
subplot(3,1,1)
plot(x,U)
hold on
plot([x(1) x(end)],[E E],'b','linewidth',2);
xlabel('position  x  (fm)','fontsize',14)
ylabel('E, U (MeV)','fontsize',14);
title('Effective Potential','fontsize',14);
set(gca,'fontsize',14);
plot([x(index) x(index)],[-1.2*max(U) 1.2*max(U)],'color',[0.7 0.7 0.7],'lineWidth',0.5);
plot([0 0],[-60 60],'r','lineWidth',3);
end
subplot(3,1,2)
set(gcf,'color',[1 1 1]);
plot(x,psi,'lineWidth',3,'DisplayName',[' a = ' num2str(a)]);
legend('show');
hold on
set(gca,'Ylim',[-1.2*max(psi) 1.2*max(psi)]);
set(gca,'fontsize',14);
xlabel('position  x  (fm)','fontsize',14)
ylabel('\psi ','fontsize',14);
title('Normalized Wavefunction','fontsize',14);
subplot(3,1,3)
Mf = 1e20;
set(gcf,'color',[1 1 1]);
plot(x,Mf.*psi,'lineWidth',3,'DisplayName',[' a = ' num2str(a)]);
legend('show');
hold on

set(gca,'Ylim',[-1.2*Mf*max(psi(index:N)) 1.2*Mf*max(psi(index:N))]);
set(gca,'fontsize',14);
xlabel('radial coordinate   \it{r}  (fm)','fontsize',14)
ylabel('\psi ','fontsize',14);
title('magnification x10^{20}','fontsize',14);
end





