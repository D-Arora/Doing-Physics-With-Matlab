function [ EB Psi r ] = qp_fh( pqn, L, m_L )
% Solution of Schrodinger Equation for hydrogen atom

% INPUTS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
r_max = 20e-10;
% if pqn > 3
%     r_max = 40e-10;
% end

Z = 1;
num = 901;     % must be odd number

% CONSTANTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
eps0 = 8.854e-12;      % F/m

% POTENTIAL WELL +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% potential energy in electron volts (eV)
r_min = 1e-15;    % default 1e-15
r = linspace(r_min,r_max, num);
dr = r(2)-r(1);
dr2 = dr^2;

% Coulomb term
K = -Z*e/(4*pi*eps0);
U_c = K./r;

% Angular momnetum term
U_L = (hbar^2*L*(L+1)/(2*me*e))./r.^2;

% Effective potential energy
U = U_c + U_L;

U(U < -2000) = -2000;     % default value 1000
U(U > 1000) = 1000;

for cn =1:(num-2)
    U_matrix(cn,cn) = U(cn+1);
end;


% SOLVE SCHRODINGER EQUATION +++++++++++++++++++++++++++++++++++++++++++++
% Make Second Derivative Matrix 
off     = ones(num-3,1);                 
SD_matrix = (-2*eye(num-2) + diag(off,1) + diag(off,-1))/dr2;

% Make KE Matrix
K_matrix = -hbar^2/(2*me*e) * SD_matrix;            

% Make Hamiltonian Matrix
H_matrix = K_matrix + U_matrix;

% Find Eignevalues E_n and Eigenfunctions psi_N 
[e_funct e_values] = eig(H_matrix);

% All Eigenvalues 1, 2 , ... n  where E_N < 0
flag = 0;
n = 1;
while flag == 0
    E(n) = e_values(n,n);
    if E(n) > 0, flag = 1; end; % if
    n = n + 1;
end  % while
%E(n-1) = [];
%n  = n-2;

% Corresponding Eigenfunctions 1, 2, ... ,n: Normalizing the wavefunction
for cn = 1 : n
   psi(:,cn) = [0; e_funct(:,cn); 0];
   area = simpson1d((psi(:,cn).* psi(:,cn))',r_min,r_max);
   psi(:,cn) = psi(:,cn)/sqrt(area);
end % for

qn = pqn-L;

Psi = psi(:,qn)';
EB = -E(pqn-L);

end

