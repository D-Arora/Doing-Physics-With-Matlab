% QMG2B.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG02C.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230512   Matlab R2021b

clear; close all;clc

% SETUP  ===========================================================
  a = 1e-9;
  NX = 201;      % [201]   must be an odd number - number of grid points
  NT = 99;   % [5e5]  number of time steps

  me = 9.10938291e-31;    % electron mass
  hbar = 1.054571726e-34; % hbar Planck's constant
  h = 6.62607015e-34; 
  e = 1.602176565e-19;    % elementary charge
  
  x1 = 0; x2 = a;
  x = linspace(x1,x2,NX); dx = x(2)-x(1);
  
  A = sqrt(30/a)*(2/pi)^3;
  k1=  pi/a;
  w = pi^2*hbar/(2*me*a^2);
  P = 2*pi/w;
  t = linspace(0,2*P,NT);
  k2 = -1i*w;
  PSI = zeros(NX,NT);

% Initial wavefunction  t = 0;
  for n = 1:2:10
   S = A*sin(n*k1*x)/n^3;
   PSI(:,1) = PSI(:,1) + S';
  end

% Time evolution of wavefunction
  for c = 2 : NT
    for n = 1:2:10
       S = A*sin(n*k1*x)/n^3 .* exp(n^2*k2*t(c));
       PSI(:,c) = PSI(:,c) + S';
    end
  end

% Probability density / check normalization at t = 0 
  probDensity = conj(PSI).*PSI;
  check = simpson1d(probDensity(:,1)',x1,x2);

 
% OUTPUT  ===========================================================
  fprintf('probability(t = 0) = %2.4f eV \n', check)
  fprintf('ang. frequency(n = 1) = %2.4e rad/s \n', w)
  fprintf('period(n = 1) = %2.4f fs \n', P/1e-15)
  fprintf('frequency(n = 1) = %2.4e Hz \n', 1/P)

% ANIMATION SETUP =======================================================
% (0 no)  (1 yes) for flag1
   flag1 = 1;    
% file name for animated gif   
    ag_name = 'agQMG2B.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.01; 
% Frame to start
    frame1 = 0;  

% GRAPHICS =========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.05 0.2 0.4]);
  set(gcf,'color','w');
  FS = 14;

  xP = x; 
  for c = 1:NT
subplot(2,1,1)      
      yP =  real(PSI(:,c));
      plot(xP./1e-9,yP,'b','LineWidth',2)
      ylim([-8e4 8e4])
      grid on
      xlabel('x  [ nm ]')
      ylabel('\Psi(x,t)')
      title('Initial wavefunction','FontWeight','normal')
      set(gca,'FontSize',FS)
subplot(2,1,2)      
      yP =  real(probDensity(:,c));
      plot(xP./1e-9,yP,'b','LineWidth',2)
      ylim([-2.5e9 2.5e9])
      grid on
      xlabel('x  [ nm ]')
      ylabel('|\Psi(x,t|^2')
      title('Prob. Density','FontWeight','normal')
      set(gca,'FontSize',FS)
      pause(0.1)

     if flag1 > 0
         frame1 = frame1 + 1;
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
      % On the first loop, create the file. In subsequent loops, append.
         if frame1 == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
    end 
  end

figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.30 0.05 0.2 0.2]);
  set(gcf,'color','w');
  FS = 14;
  xP = t; yP = real(PSI(101,:));
  plot(xP./1e-15,yP,'b','LineWidth',2)
  grid on
  xlabel('t  [ fs ]')
  ylabel('\psi(x = 0.5 nm)')
  set(gca,'FontSize',FS)



