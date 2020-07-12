% qm2DAA






clear
close all
clc



% SETUP  ==============================================================
% Mass  [kg]
  m = 9.10938356e-31; 
% Reduced Planck's constant  [J.s]
  h = 6.626070040e-34;            
  hbar = 1.054571800e-34;
% Elementary charge [C]  
  e = 1.6021766208e-19;
  
% 2D region  [m]
  Nx = 501; Ny = 501;
  Lx = 2e-9; Ly = 2e-9;
  x = linspace(-Lx,Lx,Nx); y = linspace(-Ly,Ly,Ny);
  dx = x(2)-x(1); dy = y(2)-y(1);
  
  [xx, yy] = meshgrid(x,y);
  
% Simulation time  [s] and stability condition;
  Nt = 500;
  K = 0.1;
  dt = K*m*dx^2/hbar;
  t = 0:dt:Nt*dt;

% Potential energy function
  V = zeros(Nx,Ny);
  U = (dt/hbar).*V;
  
% Wavefunction
%  psiR = zeros(Nx,Ny,Nt);
%  psiI = zeros(Nx,Ny,Nt);
  psi = zeros(Nx,Ny,Nt);
  k = 10e10;
  
  
  
% Initial wave function Gaussian
  yC = 0;%-1.5e-9;
  xC = 0;
  s = 0.5e-9; s2 = -1/(2*s^2);
  for nx = 1:Nx
    for ny = 1:Ny
      psi(nx,ny,1) = exp(1i*k*y(ny)) .* exp(s2*(x(nx)-xC).^2).* exp(s2*(y(ny)-yC).^2);
     % psiR(nx,ny,1) = cos(k*x(nx)) .* exp(s2*((x(nx)-xC).^2 + (y(ny)-yC)^2));
     % psiI(nx,ny,1) = sin(k*x(nx)) .* exp(s2*((x(nx)-xC).^2 + (y(ny)-yC)^2));
    end
  end
% Solve [2D] Schrodinger equation
     psiR = real(psi); psiI = imag(psi);
     
for nt = 2 : Nt
    for nx = 2:Nx-1
        for ny = 2:Ny-1
         
         psiI(nx,ny,nt) = psiI(nx,ny,nt-1)...
             + K*(psiR(nx+1,ny,nt-1) - 2*psiR(nx,ny,nt-1) + psiR(nx-1,ny,nt-1))...
             + K*(psiR(nx,ny+1,nt-1) - 2*psiR(nx,ny,nt-1) + psiR(nx,ny-1,nt-1)...
             - U(nx,ny)*psiR(nx,ny,nt-1));   
            
         psiR(nx,ny,nt) = psiR(nx,ny,nt-1)...
             - K*(psiI(nx+1,ny,nt-1) - 2*psiI(nx,ny,nt-1) + psiI(nx-1,ny,nt-1))...
             - K*(psiI(nx,ny+1,nt-1) - 2*psiI(nx,ny,nt-1) + psiI(nx,ny-1,nt-1)...
             + U(nx,ny)*psiI(nx,ny,nt-1));
         
        end
    end
end

 

%
figure(1)
for c = 1:10: Nt
zP = sqrt(psiR(:,:,c).^2 + psiI(:,:,c).^2);
zP = zP./max(max(zP));
pcolor(x,y,zP);
%colormap(copper)
colorbar
shading interp
axis square
pause(0.1)
zlim([0 1])
end
  
  