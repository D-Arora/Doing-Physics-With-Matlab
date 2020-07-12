% qmElectronBeamB.m

% QUANTUM MECHANICS 
% Electron Beams: Step Potential  TOTAL energy < step height (E < U0)

% Calls external function: simpson1d.m

% Calculations use S.I. units
% conversion of units:  eV <---> J   nm <---> m

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/qmElectronBeamB.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b
% 200105

clear
close all
clc

tic


global E K U0s
tic

% INPUTS ==============================================================
% Energy of particles  [100 eV];
  Es = 100;
  
% Step Potential: height of step E < U0   [200 eV]
  U0s = 200; 
  
% X axis  [-0.5  0.5 nm]
  xMin = -0.5;
  xMax = 0.5;
  
% Grid points 
  N = 15501;

% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_qm.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 
  
  
% CONSTANTS ===========================================================
  e = 1.6021766208e-19;           % Elementary charge [C]
  h = 6.626070040e-34;            % Planck constant [ J.s]
  hbar = 1.054571800e-34;         % Reduced Planck's constant
  me = 9.10938356e-31;            % Electron mass [kg]

  
% SETUP ===============================================================  
% Conversion of units:  eV --> J and nm --> m
  X = 1e-9;
  E =  e*Es;
  U0 = e*U0s;
  x = X.*linspace(xMin,xMax,N);   
  dx = x(2) - x(1);
 
% Wave number  [1/m]
  k = sqrt(2*me*E)/hbar;
  K = 2*me/hbar^2;
% Wavelength  [m]
  wL = 2*pi/k;
% Momentum  [kg.m/s]  
  p = hbar*k;  
% Angular frequency  omega  [1/s]  
  w = E/hbar;  
  
% Period T [s]
  T = 2*pi/w;
% Time grid
  nT = 201;
  t = linspace(0,1*T,nT);
  psiT = exp(-1i*w*t);
  
% POTENTIAL ENERGY  [J]
  U = zeros(N,1);
  for c = 1:N  
    U(c) = pot(x(c));
  end
     
  
% ode45
  xSpan = flip(x);
  u0 = [1,1i*k];
  options = odeset('RelTol',1e-6);
  
  [xR,u] = ode45(@FNode, xSpan, u0, options);
  
% Eignenfunction
  psi = u(:,1);
  psiS = conj(psi);
  psi_max = max(real(psi));
  
% Probability Density
  PD = conj(psi).*psi;
  PD = PD./max(PD);

% Numerical estimate of wavelength from peak #2 and peak #(N-1)
  [pks, locs] =  findpeaks(real(psi),x);
  LEN = length(locs);
  wL_N = (locs(3) - locs(2) ) / (1);
 
 
  toc


% GRAPHICS ==========================================================
 
 figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.52 0.05 0.3 0.3]);
    set(gcf,'color','w');
      
    for c = 1 : nT
      xP = xR./X;
      yP = real(psi) .* real(psiT(c));
      yP = yP./psi_max;
      plot(xP,yP,'b','lineWidth',2)
      ylim([-1.1 1.1])
      set(gca,'xtick',-0.5:0.1:0.5)
      set(gca,'fontsize',12)
      ylabel('\psi_R','fontsize',16)
      grid on 
      xlabel('x  [ nm ]','linewidth',2)
      
      if f_gif > 0 
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
         if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
         nt = nt+1;
      end  
      
       pause(0.01)
    end
    
 figure(2)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.02 0.05 0.4 0.5]);
    set(gcf,'color','w');
    
    
 subplot(3,1,1)
    xP = x./X; yP = U;
    plot(xP,yP,'k','lineWidth',2)
    hold on
    yP = (E/e).*ones(N,1);
    plot(xP,yP,'r','lineWidth',2)
    grid on
    box on
    ylabel('E & U  [ eV ]')
    
   txt = sprintf('E = %3.0f eV   U_0 = %3.0f eV  \\lambda = %3.3f nm   \\lambda_N = %3.3f nm' ...
       ,Es,U0s, wL/X,wL_N/X);
   title(txt,'fontweight','normal')
    
    set(gca,'fontsize',12)
    grid on
    
subplot(3,1,2)
    xP = xR./X; yP = real(psi);
    plot(xP,yP,'b','lineWidth',2)
    hold on
    yP = imag(psi); 
    plot(xP,yP,'r','lineWidth',2)
    box on
    legend('R','Im')
    set(gca,'fontsize',12)
    ylabel('\psi','fontsize',16)
    grid on
    
subplot(3,1,3)
    xP = xR./X; yP = PD;
    area(xP,yP)
    hold on
    plot([0 0],[0 1],'y')
    grid on
    box on
    xlabel('x  [ nm ]','linewidth',2)
    set(gca,'fontsize',12)
    ylabel('\Psi^* \Psi  [ m^{-1} ]','fontsize', 16)
    
    
 
% FUNCTIONS ==========================================================
 
 function du = FNode(x,u)
 global E K
 
  U = pot(x);
  U = U*1.6e-19;
  y = u(1);
  ydot = u(2);
  du = zeros(2,1);
  du(1) = ydot;
  du(2) = -K*(E-U)*y;
 end


 function U = pot(x)
 % Potential energy function  [eV] 
   global U0s
    
   U = 0;
    if x > 0
        U = U0s;
    end
 end
 
 
 
