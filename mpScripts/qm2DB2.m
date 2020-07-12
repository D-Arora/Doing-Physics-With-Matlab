% qm2DB.m

% QUANTUM MECHANICS
% Finite Difference Time Development Method
% [2D] Schrodinger Equation
% Propagation of a [2D] Gausssian Pulse
% Arbitrary units are used

% Pree propgation / Potential Hill / Potential Cliff
%    Single Slit / Double Slit
% The variable flagU (1,2,3,4,5) is used to change the potential energy

% Animation can be saved as a gif file (flag
% Script - modified version(Kevin Berwick)
%          'Computational Physics' by Giordano Nakanishi

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 200429 / Matlab version R2020a

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/qm2DA.htm
% Download Location
%    http://www.physics.usyd.edu.au/teach_res/mp/mscripts/


clear
close all
clc 

% Setup for saving images =============================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   flagG = 1;%f_gif = 1;
   ag_name = 'ag_qm.gif';
% Delay in seconds before displaying the next image  
   delay = 0;  
% Frame counter start
   nt = 1; 

% Setup [2D] grid  ====================================================
  N = 200;
  xG = linspace(0,1,N);
  [x, y] = meshgrid(xG);
  
  dx = xG(2) - xG(1);
  
% Time increment / Number of time Steps
  dt = 1e-5;
  nT = 81;

% Constant in S.E.
  f =  dt/(2*dx^2); 
  
  
% Potential Energy Function  ==========================================
  flagU = 5;
  U = zeros(N,N);
  G = round(N/2);
  switch flagU
     case 1                           % Free propagation
   
     case 2
       U(:,G:N) =  1e3;      % Potential Hill 
     case 3
       U(:,G:N) = -1e3;      % Potential Well  
     case 4                  % Single Slit  
       U(1:G-5,G-3:G+3)   = 15e3;
       U(G+5:N,G-3:G+3)   = 15e3;
     case 5                  % Double Slit  
       U(1:G-10,G-3:G+3)   = 15e3;
       U(G+10:N,G-3:G+3)   = 15e3;  
      U(G-5:G+5,G-3:G+3)   = 15e3;  
  end
  
 
% [2D] GAUSSIAN PULSE (WAVE PACKET) ==================================
% Initial centre of pulse
  x0 = 0.20;  y0 = 0.5;
% Initial amplitude of pulse 
  A = 10;
% Pulse width: sigma squared
  s = 5e-3;
% Wavenumber
  k0 = 50; 

% Envelope
  psiE = A*exp(-(x-x0).^2/s).*exp(-(y-y0).^2/s);
% Plane wave propagation in +X direction
  psiP = exp(1i*k0*x);
% Wavefunction
  psi1 = psiE.*psiP;
% Probability Density  
  prob1 = conj(psi1).*psi1;
% Extract Real and Imaginary parts
  R1 = real(psi1);  I1 = imag(psi1);
  
   
% UPDATE WAVEFUNCTION & GRAPHICS  ====================================
% psi1 (current value) psi2 (next value)

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.1 0.30 0.70]);
  set(gcf,'color','w');
  
for c = 1:nT
% Update real prt of wvefunction    
  R2 = psiR(N, R1, I1, dt, U, f);
  R1 = R2;

% Update imaginary part of wavefunction
  I2 = psiI(N, I1, R2, dt, U, f);

% Probability Density Function  
  prob2 = R2.^2 + I1.*I2;

  I1 = I2;

subplot(2,1,1) 
% Graph only updated after a specified number of time steps
% Probability Density values scaled to shoow very small values 
if rem(c, 8) == 0
   surf(x,y, abs(prob2.^0.05))    
   title('Probability Density (scaled)','fontweight','normal')
   xlabel('x')
   ylabel('y')
   zlabel('ps*psi');
   axis([0 1 0 1 0 2])
   view(44,55)
   axis off
   light
   lighting phong
   camlight('left')
   shading interp
  % colorbar
   
subplot(2,1,2)
   pcolor(x,y, abs((prob2).^0.05));
   %title('Probability density','fontweight','normal');
   xlabel('x')
   ylabel('y')
   axis off
   shading interp
 % colorbar
   axis square
   pause(0.00001)
end

   if flagG > 0 
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
         if nt == 9
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
         nt = nt+1;
    end   
end


figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.40 0.1 0.25 0.45]);
  set(gcf,'color','w');
 
subplot(2,1,1)
    pcolor(x,y,prob1)
    shading interp
    set(gca,'xtick',-100:100)
    set(gca,'ytick',-100:100)
    axis square
    title('Initial Probability Density','fontweight','normal')
    xlabel('x')
    ylabel('y')
    set(gca,'fontsize',10)
    
 subplot(2,1,2)
    surf(x,y,prob1)
    view(34,18)
    shading interp
    set(gca,'xtick',-100:100)
    set(gca,'ytick',-100:100)
    xlabel('x')
    ylabel('y')
    set(gca,'fontsize',10)
 
    
 figure(3)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.7 0.1 0.2 0.30]);
  set(gcf,'color','w');

  if flagU == 4 || flagU == 5
    pcolor(x,y,U)
  else
    surf(x,y,U)
    view(34,18)
  end
  shading interp
  hold on
  
  set(gca,'xtick',-100:100)
  set(gca,'ytick',-10000:10000)
  axis square
  title('Potential Energy Function','fontweight','normal')
  xlabel('x')
  ylabel('y')
  set(gca,'fontsize',10)
   
    
% FUNCTIONS  ==========================================================
% psi1 (current value at time t)--> psi2 (next value at time t+dt/2 & t+dt
% functioin pisI at times t+dt/2, t+3dt/2, ....
% function  psiR at times t+dt, t+2dt, 

function I2 = psiI(N, I1, R1, dt, U, f)
  I2 = zeros(N,N);
 % f =  dt/(2*dx^2);
  x = 2:N-1;
  y = 2:N-1;

  I2(x,y) = I1(x,y) + f*(R1(x+1,y)-2*R1(x,y)+R1(x-1,y)+R1(x,y+1)-2*R1(x,y)+R1(x,y-1))...
           - dt*U(x,y).*R1(x,y);
end

function R2 = psiR(N, R1, I1, dt, U,  f)
  R2 = zeros(N,N);
 % f = dt/(2*dx^2);
  x = 2:N-1;
  y = 2:N-1;

  R2(x,y) = R1(x,y) - f*(I1(x+1,y)-2*I1(x,y)+I1(x-1,y)+I1(x,y+1)-2*I1(x,y)+I1(x,y-1))...
           + dt*U(x,y).*I1(x,y);
end