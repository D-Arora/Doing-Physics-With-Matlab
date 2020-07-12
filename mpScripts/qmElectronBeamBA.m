% qmElectronBeamBT.m

% QUANTUM MECHANICS 
% Electron Beams: Step Potential  TOTAL energy > step height (E > U0)
% Analytical Solution of Schrodinger Equation

% Calls external function: simpson1d.m

% Calculations use S.I. units
% conversion of units:  eV <---> J   nm <---> m

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/qmElectronBeamB.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b
% 200108

clear
close all
clc

tic


% INPUTS ==============================================================
%  Energy of particles  [250 eV];
   Es = 4*200/3;
 
% Step Potential: height of step E > U0   [200 eV]
  U0s = 200; 

% Amplitude of incident wave  [10^8]
  A = 10^8;
   
% X axis  [-0.5  0.5 nm]
  xMin = -0.5;
  xMax = 0.5;
   
% Grid points and indices for region x < 0 and x > 0
  N = 5501;            % must be an odd number
  N1 = 1:floor(N/2);
  N2 = N1(end)+1:N;

  
% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_qm1.gif';
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
% nm <--> m
  X = 1e-9;
% Length of cylinder  [m]
  L = X*(xMax - xMin);
% X grid   [m]
  x = X.*linspace(xMin,xMax,N);   
  dx = x(2) - x(1);
  
% eV <---> J  
  E =  e*Es;
  U0 = e*U0s;
 
% Wave number  [1/m]
  k1 = sqrt(2*me*E)/hbar;
  k2 = sqrt(2*me*(E - U0))/hbar;
  K = 2*me/hbar^2;
  
% Wavelength lambda [m]
  wL1 = 2*pi/k1;
  wL2 = 2*pi/k2;
  
% Momentum  [kg.m/s]  
  p1 = hbar*k1; 
  p2 = hbar*k2;
  
% Velocity  [m/s]  
  v1 = hbar*k1/me; 
  v2 = hbar*k2/me;  
  
% Angular frequency  omega  [1/s]  
  w = E/hbar;  
  
% Period T [s]
  Tp = 2*pi/w;
% Time grid
  nT = 201;
  t = linspace(0,1*Tp,nT);
  psiT = exp(-1i*w*t);
  
% Amplitudes: B reflected wave / C transmitted wave
  B = (k1-k2)/(k1+k2)*A;
  C = 2*k1/(k1+k2)*A;

% Potential energy function: step potential E > U0
  U = zeros(N,1);
  U(x>0) = U0s;

% Eigenfunction  psi(x)
  psi = zeros(N,1);
  for c = 1 : N
    psi(c) = A*exp(1i*k1*x(c)) + B*exp(-1i*k1*x(c));
    if x(c) > 0
      psi(c) = C.*exp(1i*k2*x(c));
    end
  end

 % Probability density  [1/m]
   psiS = conj(psi); 
   PD = conj(psi).*psi;
   psi_max = max(real(psi));
  
% Amplitude:  numerical estimates  A -->   B <--   C -->
   Amax = max(max((abs(psi))));
   Amin = max(min((abs(psi))));
   An = 0.5*(Amax + Amin);
   Bn = 0.5*(Amax - Amin);
   Cn = 0.5*(max((real(psi(N2)) - min(real(psi(N2)) ))));
  
% Reflection and transmission coefficients
  R = B^2/A^2;
  T = 4*k1*k2/(k1+k2)^2; 
    
% Number of particles in cylinder: x < 0 num1 and x > 0 num2 
  num1 = simpson1d(PD(1:N2(1))',x(1),x(N2(1)));
  num2 = simpson1d(PD(N2(1):N2(end))',x(N2(1)),x(N2(end)));
    
% Probabilty current   [1/s]
  J = psiS.*gradient(psi,dx) - psi.*gradient(psiS,dx);
  J = -(1i*hbar/(2*me)) .* J;
  J_AC = mean(J);
  J_AB = v1*A^2 - v1*B^2;
  J_C = v2*C^2;
  
% Numerical estimate of wavelength from peak #2 and peak #(N-1)
  [~, locs] =  findpeaks(real(psi(N1)),x(N1));
   wL1_N = (locs(2) - locs(1) ) / (1);
 
  [~, locs] =  findpeaks(real(psi(N2)),x(N2));
   wL2_N = (locs(2) - locs(1) ) / (1);
   

% GRAPHICS ==========================================================

figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.05 0.05 0.25 0.3]);
    set(gcf,'color','w');

    for c = 1 : nT
    subplot(2,1,1)  
      xP = x./X; yP = real(psi .* psiT(c));
      yP = yP./psi_max;
      plot(xP,yP,'b','lineWidth',2)
      ylim([-1.2 1.2])
      grid on 
      set(gca,'xtick',-0.5:0.1:0.5)
      set(gca,'fontsize',12)
      ylabel('\psi_R','fontsize',16)
      
    subplot(2,1,2)  
      yP = imag(psi .* psiT(c));
      yP = yP./psi_max;
      plot(xP,yP,'r','lineWidth',2)
      ylim([-1.2 1.2])
      grid on 
      set(gca,'xtick',-0.5:0.1:0.5)
      xlabel('x  [ nm ]','linewidth',2)
      set(gca,'fontsize',12)
      ylabel('\psi_I','fontsize',16)  
      
      
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
    set(gcf,'position',[0.32 0.05 0.25 0.6]);
    set(gcf,'color','w');
    xP = x./X;
    
subplot(4,1,1)   
    yP = U;
    plot(xP,yP,'b','lineWidth',2)
    hold on
    yP = Es.*ones(N,1);
    plot(xP,yP,'r','lineWidth',2)
    grid on
    ylabel('E   U') 
    set(gca,'fontsize',12)
     
    
subplot(4,1,2)   
    yP = real(psi);
    plot(xP,yP,'b','lineWidth',2)
    hold on
    yP = imag(psi);
    plot(xP,yP,'r','lineWidth',2)
    set(gca,'xtick',-0.5:0.1:0.5)
    grid on
    set(gca,'fontsize',12)
    ylabel('\psi_R  \psi_I','fontsize',16)  
    
subplot(4,1,3)
   yP = PD;
   plot(xP,yP,'b','lineWidth',2)
   set(gca,'xtick',-0.5:0.1:0.5)
   grid on
   set(gca,'fontsize',12)
   ylabel('\psi^* \psi','fontsize',16)  
   
subplot(4,1,4)
   yP = J;
   plot(xP,yP,'k','lineWidth',2)
   set(gca,'xtick',-0.5:0.1:0.5)
   grid on
   ylim([0 1.2*J_AC])
   ylabel('\psi_I','fontsize',16)  
   xlabel('x  [ nm ]')
   ylabel('J  [ s_{-1} ]')
   set(gca,'fontsize',12)
   
   
% DISPLAY RESULTS -----------------------------------------------------

figure(7)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.58 0.05 0.25 0.75]);
   set(gcf,'color','w');

   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -18;
   
   txt = 'ELECTRON BEAM: STEP POTENTIAL  E > U_0';
   Htext = text(0,hh,txt,'fontsize',14);
   set(Htext,'color','b')
   
   hh = hh+dh;
   txt = sprintf('Electron energy  E = %3.0f  eV', Es);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Step Height  U_0 = %3.0f  eV', U0s);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('\\omega = %3.2e s^{-1}         period = %3.2e s', w, Tp);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('A = %3.2e               A_N = %3.2e', A, An);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('B = %3.2e               B_N = %3.2e', B, Bn);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('C = %3.2e               C_N = %3.2e', C, Cn);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh; 
   txt = sprintf('\\lambda_1 = %3.4f nm            \\lambda_{N1} = %3.4f nm', wL1/X,wL1_N/X);
   text(0,hh,txt,'fontsize',12)
  
   hh = hh+dh; 
   txt = sprintf('\\lambda_2 = %3.4f nm            \\lambda_{N2} = %3.4f nm', wL2/X,wL2_N/X);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('v_1 = %3.2e m.s^{-1}       v_2 = %3.2e m.s^{-1}', v1, v2);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Electrons: x < 0  N_1 = %3.2e   x > 0  N_2 = %3.2e', num1, num2);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Prob current  J = %3.2e  s^{-1}', J_AC);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Reflection coefficient  R = %3.2f ', R);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Transmission coefficient  T = %3.2f ', T);
   text(0,hh,txt,'fontsize',12)
  
   axis off
   
toc


 
 
