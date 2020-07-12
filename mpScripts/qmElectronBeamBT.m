% qmElectronBeamBT.m

% QUANTUM MECHANICS 
% Electron Beams: Step Potential  TOTAL energy > step height (E > U0)

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
% Energy of particles  [250 eV];
  Es = 4*200/3;
  
% Step Potential: height of step E > U0   [200 eV]
  U0s = 200; 
  
% X axis  [-0.5  0.5 nm]
  xMin = -0.5;
  xMax = 0.5;
  
% Grid points and indices for region x < 0 and x > 0
  N = 5501;
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
% Conversion of units:  eV --> J and nm --> m
  X = 1e-9;
  E =  e*Es;
  U0 = e*U0s;
  x = X.*linspace(xMin,xMax,N);   
  dx = x(2) - x(1);
 
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
  
  
% POTENTIAL ENERGY  [J]
  U = zeros(N,1);
  for c = 1:N  
    U(c) = pot(x(c));
  end
     
  
% ode45
  xSpan = flip(x);
  u0 = [1,1i*k1];
  
  options = odeset('RelTol',1e-8);
  [xR,u] = ode45(@FNode, xSpan, u0, options);
 
% Eignenfunction and complex conjugate
  psi = flip(u(:,1));
  psiS = conj(psi);
   
  psi_max = max(real(psi));
  
  psiR = 1.333e8*real(psi)./max(real(psi));
  psiI = 1.333e8*imag(psi)./max(imag(psi));
  
  psi = psiR + 1i*psiI;
  psiS = conj(psi);
  
  

  
% Probability Density
  PD = conj(psi).*psi;
  
% Numerical estimate of wavelength from peak #2 and peak #(N-1)
  [~, locs] =  findpeaks(real(psi(N1)),x(N1));
   wL1_N = (locs(2) - locs(1) ) / (1);
  [~, locs] =  findpeaks(real(psi(N2)),x(N2));
  wL2_N = (locs(2) - locs(1) ) / (1);
  
% Amplitude:  numerical estimates  A -->   B <--   C -->
   Amax = max(max((abs(psi))));
   Amin = max(min((abs(psi))));
   An = 0.5*(Amax + Amin);
   Bn = 0.5*(Amax - Amin);
   Cn = 0.5*(max((real(psi(N2)) - min(real(psi(N2)) ))));

% Relection and Transmission coeffficents   
  R = ((k1-k2)/(k1+k2))^2;
  Rn = Bn*Bn/(An*An);
  T  = 1 - R;
  Tn = 1 - Rn;

% Number of particles in cylinder: x < 0 num1 and x > 0 num2 
  num1 = simpson1d(PD(1:N2(1))',x(1),x(N2(1)));
  num2 = simpson1d(PD(N2(1):N2(end))',x(N2(1)),x(N2(end)));
  
  
% Probability current  [1/s]
  J = psiS.*gradient(psi,dx) - psi.*gradient(psiS,dx);
  J = -(1i*hbar/(2*me)) .* J;
   
  toc


% GRAPHICS ==========================================================
 
figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.05 0.05 0.25 0.3]);
    set(gcf,'color','w');

    for c = 1 : nT
    subplot(2,1,1)  
      xP = x./X; yP = real(psi .* psiT(c));
      plot(xP,yP,'b','lineWidth',2)
      grid on 
      set(gca,'xtick',-0.5:0.1:0.5)
      set(gca,'fontsize',12)
      ylabel('\psi_R','fontsize',16)
      ylim(max(yP).*[-1.2 1.2])
      
    subplot(2,1,2)  
      yP = imag(psi .* psiT(c));
      plot(xP,yP,'r','lineWidth',2)
      ylim(max(yP).*[-1.2 1.2])
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
   ylim([0 1.2*J(1)])
   ylabel('\psi_I','fontsize',16)  
   xlabel('x  [ nm ]')
   ylabel('J  [ s_{-1} ]')
   set(gca,'fontsize',12)
    
%%
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
   txt = sprintf('A_N = %3.2e',An);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('B_N = %3.2e',Bn);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('C_N = %3.2e',Cn);
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
   txt = sprintf('Prob current  J = %3.2e  s^{-1}', J(1));
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Reflection coefficient  R = %3.2f   R_N = %3.2f ', R, Rn);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Transmission coefficient  T = %3.2f   T_N = %3.2f', T, Tn);
   text(0,hh,txt,'fontsize',12)
  
   axis off
   
 %%  
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
 
 
 
