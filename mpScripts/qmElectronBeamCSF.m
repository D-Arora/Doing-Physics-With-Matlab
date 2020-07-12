% qmElectronBeamCSF.m

% QUANTUM MECHANICS 
% Electron Beams: Sigmoid Shaped Potential
% Shape of Sigmoid function determined by b value

% Calls external function: simpson1d.m

% Calculations use S.I. units
% conversion of units:  eV <---> J   nm <---> m

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/qmElectronBeamC.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b
% 200111

clear
close all
clc

tic


global E K U0s b
tic

% INPUTS ==============================================================
% E > U0s
% Energy of particles  [250 eV];
  Es =  210; %4*200/3;
  
% Step Potential: height of step E > U0   [200 eV]
% Steepness of sigmoid functionb = 5e10 / step function  b = 100e10]
  U0s = 200; 
  b   = 10e10;
  
% X axis  [-0.5  0.5 nm]
  xMin = -0.6;
  xMax = 1.0;
  
% Grid points and indices for region x < 0 and x > 0
  N = 5501;

  
% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
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
  flagE = 1;
  if Es < U0s; flagE = 0; end

% Conversion of units:  eV --> J and nm --> m
  X = 1e-9;
  E =  e*Es;
  U0 = e*U0s;
  x = X.*linspace(xMin,xMax,N);   
  dx = x(2) - x(1);
  
  % Index when x = 0
  N0 = find(x>0,1);
  if mod(N0,2) == 0, N0 = N0-1; end
      
  N1 = 1:N0;
  N2 = N0:N;
  
  k1 = sqrt(2*me*E)/hbar;
  K = 2*me/hbar^2;
   
% Angular frequency  omega  [1/s]  
  w = E/hbar;  
% Frequency
  f = w/(2*pi);
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
    
  psiMax = max(real(psi));
  
  psiR = 1.333e8*real(psi)./max(real(psi));
  psiI = 1.333e8*imag(psi)./max(imag(psi));
  
  psi = psiR + 1i*psiI;
  psiS = conj(psi);
  
% Probability Density
  PD = conj(psi).*psi;

  
if flagE == 1 
   
% THEORETICAL VALUES:  region 1 U = 0 / region 2 U = U0
% Kinetic Energy eV and J
  KEs(1) = Es;
  KEs(2) = Es - U0s;
  KE = e .* KEs;
  
% Momentum  [kg.m/s]
  p = sqrt(2*me*KE);

% Velcoity  [m/s]
  v = p./me;
    
% Wave number  [1/m]
  k = p./hbar;

% Wavelength   [m]
  wL = 2*pi./k;
  
% Numerical estimate of wavelength from peak #2 and peak #(N-1)
  [~, locs] =  findpeaks(real(psi(N1)),x(N1));
   wL1_N = (locs(2) - locs(1) ) / (1);
  [~, locs] =  findpeaks(real(psi(N2)),x(N2));
  wL2_N = (locs(end) - locs(end-1) ) / (1);
 
% Amplitude:  numerical estimates  A -->   B <--   C -->
   Amax = max(max((abs(psi))));
   Amin = max(min((abs(psi))));
   An = 0.5*(Amax + Amin);
   Bn = 0.5*(Amax - Amin);
   Cn = 0.5*(max((real(psi(N2)) - min(real(psi(N2)) ))));

% Reflection and Transmission coefficients   
  Rn = Bn*Bn/(An*An);
  Tn = 1 - Rn;

% Number of particles in cylinder: x < 0 num1 and x > 0 num2 
  num(1) = simpson1d(PD(1:N2(1))',x(1),x(N2(1)));
  num(2) = simpson1d(PD(N2(1):N2(end))',x(N2(1)),x(N2(end)));
  
% Number density
  n(1) = num(1)/abs(x(1));
  n(2) = num(2)/abs(x(N));
  
% Probability current  [1/s]
  J = psiS.*gradient(psi,dx) - psi.*gradient(psiS,dx);
  J = -(1i*hbar/(2*me)) .* J;
    
  
end

toc


% GRAPHICS ==========================================================
 
figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.05 0.05 0.25 0.5]);
    set(gcf,'color','w');

    Nind = 1:10:N;
    
        
    for c = 1 : nT
    subplot(3,1,1)  
      xP = x(Nind)./X; yP = real(psi(Nind) .* psiT(c));
      plot(xP,yP,'b','lineWidth',2)
      xlim([xMin xMax])
      ylim([-2e8 2e8])
      set(gca,'xtick',xMin:0.2:xMax)
      set(gca,'fontsize',12)
      ylabel('\psi_R','fontsize',16)
      grid on      
      
    subplot(3,1,2)  
      yP = imag(psi(Nind) .* psiT(c));
      plot(xP,yP,'r','lineWidth',2)
      xlim([xMin xMax])
      ylim([-2e8 2e8])
      set(gca,'xtick',xMin:0.2:xMax)
      xlabel('x  [ nm ]','linewidth',2)
      set(gca,'fontsize',12)
      ylabel('\psi_I','fontsize',16)  
      grid on
     
     subplot(3,1,3)  
      yP = PD(Nind);
      %plot(xP,yP,'r','lineWidth',2)
      area(xP,yP)
      xlim([xMin xMax])
      set(gca,'xtick',xMin:0.2:xMax)
      xlabel('x  [ nm ]','linewidth',2)
      set(gca,'fontsize',12)
      ylabel('\psi^* \psi','fontsize',16)  
      grid on
      
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
    
%%
if flagE == 1
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
    xlim([xMin xMax])
    set(gca,'xtick',xMin:0.2:xMax)
    grid on
    ylabel('E   U') 
    set(gca,'fontsize',12)
     
    
subplot(4,1,2)   
    yP = real(psi);
    plot(xP,yP,'b','lineWidth',2)
    hold on
    yP = imag(psi);
    plot(xP,yP,'r','lineWidth',2)
    xlim([xMin xMax])
    set(gca,'xtick',xMin:0.2:xMax)
    grid on
    set(gca,'fontsize',12)
    ylabel('\psi_R  \psi_I','fontsize',16)  
    
subplot(4,1,3)
   yP = PD;
  % plot(xP,yP,'b','lineWidth',2)
    area(xP,yP)
    xlim([xMin xMax])
    set(gca,'xtick',xMin:0.2:xMax)
    grid on
   set(gca,'fontsize',12)
    ylabel('\psi^* \psi','fontsize',16)  
   
subplot(4,1,4)
   yP = J;
   plot(xP,yP,'k','lineWidth',2)
   xlim([xMin xMax])
   set(gca,'xtick',xMin:0.2:xMax)
   grid on
   ylim([0 1.2*J(1)])
   ylabel('\psi_I','fontsize',16)  
   xlabel('x  [ nm ]')
   ylabel('J  [ s^{-1} ]')
   set(gca,'fontsize',12)
    
%%

figure(7)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.58 0.05 0.35 0.75]);
   set(gcf,'color','w');

   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -18;
   
   txt = 'ELECTRON BEAM: SIGMOID  E > U_0';
   Htext = text(0,hh,txt,'fontsize',14);
   set(Htext,'color','b')
   
   hh = hh+dh;
   txt = sprintf('Electron energy  E = %3.0f  eV', Es);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Step Height  U_0 = %3.0f  eV', U0s);
   text(0,hh,txt,'fontsize',12)
   
   txt = sprintf('Sigmoid Function  b = %3.2e  ', b);
   text(52,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('KE_1 = %3.0f eV      KE = %3.0f eV', KEs(1), KEs(2));
   text(0,hh,txt,'fontsize',12)
      
   hh = hh+dh;
   txt = sprintf('\\omega = %3.2e s^{-1}        period T_P = %3.2e s        frequency f = %3.2e Hz' ... 
       , w, Tp, f);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('A_N = %3.2e',An);
   text(0,hh,txt,'fontsize',12)
   
   
   txt = sprintf('B_N/A_N = %3.2f',Bn/An);
   text(40,hh,txt,'fontsize',12)
   
   txt = sprintf('C_N/A_N = %3.2f',Cn/An);
   text(80,hh,txt,'fontsize',12)
   
    hh = hh+dh; 
   txt = sprintf('\\lambda_1 = %3.4f nm       \\lambda_2 = %3.4f nm', ...
       wL(1)/X, wL(2)/X);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh; 
   txt = sprintf('\\lambda_{N1} = %3.4f nm     \\lambda_{N2} = %3.4f nm', ...
       wL1_N/X, wL2_N/X);
   text(0,hh,txt,'fontsize',12)
    
   hh = hh+dh;
   txt = sprintf('v_1 = %3.2e m.s^{-1}       v_2 = %3.2e m.s^{-1}', v(1), v(2));
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Number density: x < 0  n_1 = %3.2e m^{-1}   x > 0  n_2 = %3.2e  m^{-1}'  ...
       , n(1), n(2));
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Prob flux  J = %3.2e  s^{-1}', J(1));
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Reflection coefficient  R_N = %3.2f ',Rn);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Transmission coefficient  T_N = %3.2f',Tn);
   text(0,hh,txt,'fontsize',12)
  
   axis off
   
end   
   
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
    global U0s b
   
    U = U0s./(1 + exp(-b.*x));
 end
 
 
 
