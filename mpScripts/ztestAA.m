


clear 
close all
clc


% % Elementary charge [C]  
%    e = 1.6021766208e-19;
% %   
% % % Planck constant [ J.s] / Reduced Planck's constant
% %   h = 6.626070040e-34;            
% %   hbar = 1.054571800e-34; 
% % 
% % % Electron mass [kg / amu / MeV] 
% %   m = 9.10938356e-31; 
% %   
% %   L = 0.1e-9;
% %   N = 51;
% %   
% % %   x = linspace(0,L,N);
% % %   dx = x(2) - x(1);
% %   
% %   K = -hbar^2/(2*m);
% %   
%   
% % % Three-point finite-difference representation of Laplacian
% % Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) ...
% % + diag(ones((N-1),1),-1))/(dx^2);
% % % Next modify Lap so that it is consistent with f(0) = f(L) = 0
% % Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; % So that f(0) = 0
% % Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0;% So that f(L) = 0  
% % 
% % H = K.*Lap;
% % [V,E] = eig(H);
% % 
% % E = diag(E)./e
% m = 1; hbar = 1;
% % Parameters for solving problem in the interval 0 < x < L
% L = pi; % Interval Length
% N = 100; % No. of coordinate points
% x = linspace(0,L,N)'; % Coordinate vector
% dx = x(2) - x(1); % Coordinate step
% % Two-point finite-difference representation of Derivative
% D=(diag(ones((N-1),1),1)-diag(ones((N-1),1),-1))/(2*dx);
% % Next modify D so that it is consistent with f(0) = f(L) = 0
% D(1,1) = 0; D(1,2) = 0; D(2,1) = 0; % So that f(0) = 0
% D(N,N-1) = 0; D(N-1,N) = 0; D(N,N) = 0; % So that f(L) = 0
% % Three-point finite-difference representation of Laplacian
% Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) ...
% + diag(ones((N-1),1),-1))/(dx^2);
% % Next modify Lap so that it is consistent with f(0) = f(L) = 0
% Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; % So that f(0) = 0
% Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0;% So that f(L) = 0
% 
% % Total Hamiltonian where hbar=1 and m=1
% %hbar = 1; m = 1; 
% 
% H = -(1/2)*(hbar^2/m)*Lap;
% % Solve for eigenvector matrix V and eigenvalue matrix E of H
% [V,E] = eig(H);
% % Plot lowest 3 eigenfunctions
% plot(x,V(:,3),'r',x,V(:,4),'b',x,V(:,5),'k'); shg;
% E % display eigenvalue matrix
% E = diag(E) % display a vector containing the eigenvalues  
 

% % Parameters for solving problem in the interval -L < x < L
% % PARAMETERS:
% L = 5; % Interval Length
% N = 1000; % No of points
% x = linspace(-L,L,N)'; % Coordinate vector
% dx = x(2) - x(1); % Coordinate step
% % POTENTIAL, choose one or make your own
% U = 1/2*100*x.^(2); % quadratic harmonic oscillator potential
% %U = 1/2*x.^(4); % quartic potential
% % Finite square well of width 2w and depth given
% %w = L/50;
% %U = -500*(heaviside(x+w)-heaviside(x-w));
% % Two finite square wells of width 2w and distance 2a apart
% %w = L/50; a=3*w;
% %U = -200*(heaviside(x+w-a) - heaviside(x-w-a) ...
% % + heaviside(x+w+a) - heaviside(x-w+a));
% % Three-point finite-difference representation of Laplacian
% % using sparse matrices, where you save memory by only
% % storing non-zero matrix elements
% e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
% % Total Hamiltonian
% hbar = 1; m = 1; % constants for Hamiltonian
% H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
% % Find lowest nmodes eigenvectors and eigenvalues of sparse matrix
% nmodes = 4; options.disp = 0;
% [V,E] = eigs(H,nmodes,'sa',options); % find eigs
% [E,ind] = sort(diag(E));% convert E to vector and sort low to high
% V = V(:,ind); % rearrange corresponding eigenvectors
% % Generate plot of lowest energy eigenvectors V(x) and U(x)
% Usc = U*max(abs(V(:)))/max(abs(U)); % rescale U for plotting
% plot(x,V,x,Usc,'--k'); % plot V(x) and rescaled U(x)
% % Add legend showing Energy of plotted V(x)
% lgnd_str = [repmat('E = ',nmodes,1),num2str(E)];
% legend(lgnd_str) % place lengend string on plot
% shg



% N = 57;
% L = 1;
% x = linspace(0,L,N+2);
% dx = x(2)-x(1);
% 
% 
% off = ones(N-1,1);
% A = -2*eye(N) + diag(off,1) + diag(off,-1);
% 
% A = A./dx^2;
% 
% [eignFN, eignV] = eig(A);
% 
% K = diag(eignV);
% 
% m = 1;
% y = zeros(N+2,1);
% y(2:N+1) = eignFN(:,m);
% y = y ./max(y);
% 
% 
% figure(1)
% plot(x,y)


% INPUTS  =============================================================
% Boundary conditions:
%     fixed fixed BC = 1 / fixed free BC = 2 / free free BC = 3 
BC = 3; 
% Mode number for graphical outputs  BC = 3 --> m > 1
m = 5;
% Number of grid points
N = 599;
% Length of rod [m]
L = 1;

% CALCULATIONS  ======================================================


% Spatial domian  [m]
%  x = (0:N+1).*(L/(N+1));
  x = linspace(0,L,N+2);
  dx = x(2) - x(1);
% Eigenvalue Matrix A: eigenfunctions (eignFN) / eigenvalues (eignV) 
  off = ones(N-1,1);
 
  A = - diag(off,-1) + 2*eye(N) - diag(off,1) ;
  
  
  if BC == 2; A(N,N) = 1; end               % fixed free
  if BC == 3; A(1,1) = 1; A(N,N) = 1; end   % free free

  [eignFN, eignV] = eig(A);
  
% Spatial Wavefunction  US
  US = zeros(N+2,1);
  US(2:N+1) = eignFN(:,m);
  US = US ./max(US);

  if BC == 2; US(N+2) = US(N+1); end
  if BC == 3; US(1)   = US(2); US(N+2) = US(N+1); end
  
   K = diag(eignV);
   % Time dependent wavefunction UT
  % propagation constant  [1/m]
  %  k = sqrt(eignV(m,m)/dx^2);
    k = sqrt(K(m))/dx;
    wL = 2*pi/k
    
    
figure(1)
plot(x,US)
grid on

