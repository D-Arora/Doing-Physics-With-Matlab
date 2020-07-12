%coulomb potential
%+WKB approximation to solve the non analytical SE.(centrifugal+nuclear
%potentials
% based on the se_fdm_alpha.m script of  Ian Cooper,School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% modification made by Badr-Eddine Elfatehy Chfira: badreddine.el@um.es

close all
clear all
clc

% INPUTS =================================================================
E = input(' E(MeV) =     ');
disp('  ')
A = input(' A =     ');
disp('  ')
Z = input(' Z =     ');
disp('   ')
l = input(' l (if l=0 then Ucentrifugal=0) =     '); % angular momentum quantum number
disp('   ')
U0 = input(' U0 (if Uo=0 then Saxon-Woods pot=0) =     '); 
disp('   ')

L = 120;                    % max radial coordinate ([fm] 
N = 999999;                   % number of x values - must be ODD

% CONSTANTS ===============================================================
hbar = 1.054571818e-34;           % J.s
e = 1.602e-19;              % C
mA = 6.64465675e-27;        % kg

eps0 = 8.854187e-12;        % permittivity of free space
k = 8.9876e9;                %coulomb's constant
Ese = 1.6e-13;                       % Energy scaling factor  
Lse = 1e-15;                         % Length scaling factor  

% SETUP CALCULATIONS ==============================
% x coodinates
x_min = 0;
x_max = L;
x = linspace(x_min,x_max,N);
dx = x(2)-x(1);
%Values for Saxon-Woods Potential
R0 = 1.2.* (4^(1/3) + (A-4)^(1/3)); %classical turning point r1, when U=E
a = 0.5; %0.662fm surface thickness of the nucleus
%=====================WKB METHOD==========================================

% case l =/= 0
%Turning points R1 and R2

if l ~= 0
re2=solve(@(z) (l+1./2).^2.*hbar.^2./(2.*mA.*Ese.*(z.*Lse).^2)+2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* z)==E);%Turning points >R0
R1 = double((max(re2)));            % this is the classical turning point r2
b = R1-R0;

%Effective potential, U = U_centrifugal+U_Saxon-Woods+U_Coulombian, respectively
for n = 1:N
   Un(n) =  -U0./(1+exp((x(n)-R0)./a));                                    % Saxon Woods Potential
   Uo(n) = 0;                  % Centrifugal Potential
   Uc(n) = 0;
   if x(n) > R0
       Uc(n) = + 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n));         %Coulomb potential (r>R0)
       Uo(n) = (l+1./2).^2.*hbar.^2./(2.*mA.*Ese.*(x(n).*Lse).^2);
   end;
   U(n) = Un(n)+Uc(n)+Uo(n);                                               %Effective Potential 
end
%======Graphics============================================================
plot(x,Un);
hold on;
plot(x,Uc);
hold on;
plot(x,Uo);
hold on;
plot(x,U); 
line([0,L],[E,E]);  %Plot line of alpha particles' energy
legend('Saxon-Woods','Coulomb ','Centrifugal','U_eff','E');
hold on;

xlim([0.1 L]);
ylim([-60 60]);

%Half life calculation for l=/=0 using WKB approximation :

Schrct= Lse^2* Ese*(2*mA/hbar^2);                                          % 2m/hbar^2
P2   =@(g) abs(sqrt(Schrct*(abs(E - ((l+1./2).^2.*hbar.^2./(2.*mA.*Ese.*(g.*Lse).^2)-U0./(1+exp((g-R0)./a))+2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* g))))));   % Momentum P(r)
P3 = @(g) 1./abs(sqrt(Schrct*(abs(E - ((l+1./2).^2.*hbar.^2./(2.*mA.*Ese.*(g.*Lse).^2)-U0./(1+exp((g-R0)./a))+2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* g))))));
psi  = integral( P2,R0, R1);                                               % Integral for psi
norm = sqrt(Schrct.*integral(P3,0,R0));
TT = exp(-2.*psi);                                                            %Transmission Coefficient,the square of it give the probability 
vA = sqrt(2.*(E+U0./(1+exp((-R0)./a))).*Ese./mA);                                             % velocity of escaped alpha particle
f =  vA ./ (2.*(R0).*Lse);                                              % frequency at which particles strike barrier on the RHS of well 
P = TT;                                                                 %Probability
lambda = f.*P;                                                             % Decay Constant
T = log(2)./lambda;                                                        % Half Life
T1 =table(l,R0,R1,b,vA,f,P,lambda,T)                                 % table of results for every 
end

% case l = 0:
if l == 0   
re1=solve(@(p) 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* p)==E);%Turning point >R0
R1 = double(re1(1));             % this is the classical turning point r2
b = R1-R0;
%Effective potential, U_centrifugal=0
for n = 1:N
   Un(n) =  -U0./(1+exp((x(n)-R0)./a));                                    % Saxon Woods Potential
   Uo(n) = 0.0;                                                            % Centrifugal Potential
   Uc(n) = 0.0;
   if x(n) > R0
       Uc(n) = + 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* x(n));         %Coulomb potential (r>R0)
   end;
   U(n) = Un(n)+Uc(n)+Uo(n);                                               %Effective Potential 
end
%==========Graphics========================================================
plot(x,Un);
hold on;
plot(x,Uc);
hold on;
plot(x,Uo);
hold on;
plot(x,U); 
line([0,L],[E,E]);  %Plot line of alpha particles' energy
legend('Saxon-Woods','Coulomb ','Centrifugal','U_eff','E');
hold on;

xlim([0.1 L]);
ylim([-60 60]);

%Half life calculation for l==0 using WKB approximation :


re2=solve(@(z) +2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* z)==E);%Turning points >R0
R2 = double((max(re2)));            % this is the classical turning point r2
b = R2-R0;  
Schrct= Lse^2* Ese*(2*mA/hbar^2);                                          % 2m/hbar^2
P2   =@(g) abs(sqrt(Schrct*(abs(E - (-U0./(1+exp((g-R0)./a))+2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* g))))));   % Momentum P(r)
P3 = @(g) 1./abs(sqrt(Schrct*(abs(E - (-U0./(1+exp((g-R0)./a))+2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* g))))));
psi  = integral( P2,R0, R2);                                               % Integral for psi
TT = exp(-2.*psi);                                                            %Transmission Coefficient,the square of it give the probability 
vA = sqrt(2.*(E+U0./(1+exp((-R0)./a))).*Ese./mA);                                             % velocity of escaped alpha particle
f =  vA ./ (2.*(R0).*Lse);                                              % frequency at which particles strike barrier on the RHS of well 
P = TT;                                                                 %Probability
lambda = f.*P;                                                             % Decay Constant
T = log(2)./lambda;                                                        % Half Life
T2 =table(l,R0,R2,b,vA,f,P,lambda,T)                                 % table of results for every     

end   
