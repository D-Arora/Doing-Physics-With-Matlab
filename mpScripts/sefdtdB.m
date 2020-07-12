% se_fdtd_A.m


% Solving the Time Dependent Schrodinger Equation using the FDTD Method

% Ian Cooper
% School of Physics, University of Sydney
% documentation: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts


%181014  Matlab 2018b

close all
clear 
%clc
tic

% ========================================================================
% PARAMETERS 
% ========================================================================
Nx = 1001;               %  ODD number - number of grid points [1001]
Nt = 500;               %  Number of time steps [10000]
L = 4e-9;                %  Width of X domain [4e-9] 
me = 9.10938291e-31;     % electron mass
hbar = 1.054571726e-34;  % hbar Planck's constant
e = 1.602176565e-19;     % elementary charge


% =======================================================================
% SETUP
% ======================================================================
  x = linspace(0,L,Nx); dx = x(2) - x(1);
  k1 = -hbar^2/(2*me);
  C1 = 1/10;
  dt = C1 * 2 * me * dx^2 / hbar;
  C2 = e*dt/hbar;
  C3 = -hbar^2 / (2 * me * dx^2 * e);
   
% Wavefunction y yR yI / prob density pd  / potential energy U 
  yR = zeros(Nt,Nx); yI = yR; y = yR + 1i.*yR; pd = yR;

% Potential energy function
  U = zeros(1,Nx);

% ========================================================================
% INITIAL WAVE PACKET
% ========================================================================
   nx0 = round(Nx/4);     % pulse centre   [round(Nx/4)]
   s = L/25;              % pulse width    [L/25]
   wL = 1.6e-10;
  
% Real and Imaginary parts of wavefunction  
    %yyyR(1,:) = exp(-0.5.*((x-x(nx0))./s).^2).*cos(2*pi.*(x-x(nx0))./wL);
%   yI(1,:) = exp(-0.5.*((x-x(nx0))./s).^2).*sin(2*pi.*(x-x(nx0))./wL);
 
   yR(1,:) = exp(-0.5*((x-x(nx0))./s).^2).*cos(2*pi.*(x-x(nx0))./wL);
   yI(1,:) = exp(-0.5*((x-x(nx0))./s).^2).*sin(2*pi.*(x-x(nx0))./wL);
   for nx = 1:Nx   
   yyR(nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2).*cos(2*pi*(x(nx)-x(nx0))/wL);
  % yyR(nx) =  exp(-0.5*((x(nx)-x(nx0))/s)^2)*cos(2*pi*(x(nx)-x(nx0))/wL);
  % yyI(1,nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2)*sin(2*pi*(x(nx)-x(nx0))/wL);
   end
  yR(1,200:209)
  yyR(200:209)
  
   
   
   
   % Normalize wavefunction  
   [yR(1,:), yI(1,:),pd(1,:)] = normalize(yR(1,:),yI(1,:),L);
   psiNorm = simpson1d(pd(1,:),0,L);
   
% Kinetic energy
   K1 = real(KE(yR(1,:),yI(1,:),dx,k1,L)/e);
 
   
% =====================================================================
% Solve Schrodinger Equation: FDTD Method
% =====================================================================
for nt = 1 : Nt
   for nx = 2 : Nx - 1
      yR(nt+1,nx) = yR(nt,nx) - C1*(yI(nt,nx+1) - 2*yI(nt,nx)+yI(nt,nx-1)) ...
          + C2*U(nx)*yI(nt,nx);
   end
   
   for nx = 2 : Nx-1
      yI(nt+1,nx) = yI(nt,nx) + C1*(yR(nt,nx+1) - 2*yR(nt,nx)+yR(nt,nx-1)) ...
          - C2*U(nx)*yR(nt,nx);
   end
end  




% =====================================================================
% GRAPHICS
% =====================================================================

%%
  FS = 12;
  
figure(1)
  pos = [0.05 0.05 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
 % set(gcf,'color','w'); 
  xP = x;
 
  
  for cc = Nt
  subplot(3,1,1)
  yP = yR(cc,:);
  plot(xP,yP,'b','linewidth',2)
  grid on
  box on
  set(gca,'fontsize',FS)
  ylim([-1e5 1e5])
  
  subplot(3,1,2)
  yP = yI(cc,:);
  plot(xP,yP,'r','linewidth',2)
  grid on
  box on
  set(gca,'fontsize',FS)
  ylim([-1e5 1e5])
  
  subplot(3,1,3)
  yP = yR(cc,:).^2 + yI(cc,:).^2;
  plot(xP,yP,'k','linewidth',2)
  
  grid on
  box on
  set(gca,'fontsize',FS)
 ylim([0 1e10])
  
  
  end
  
  
  %%
  
  % ===================================================================
  % FUNCTIONS
  % ===================================================================
  
 function [f,g,pd] = normalize(F,G,L)
    M = F.^2 + G.^2;
    A = simpson1d(M,0,L);
    f = F ./ sqrt(A); g = G ./ sqrt(A);
    pd = f.^2 + g.^2;
 end

 function [f] = KE(F,G,dx,k1,L)
   Y = F + 1i.*G;
   secDerY = 4*del2(Y,dx);
   fn = conj(Y).*secDerY;
   f = k1*simpson1d(fn,0,L);
 end
 