% qmElectronBeamB.m





clear
close all
clc


global E K
tic

% INPUTS ==============================================================

% Energy of particles  [eV];
  Es = 100;
% X axis  [nm]
  xMin = 0;
  xMax = 0.5;
% Grid points 
  N = 15001;

  % CONSTANTS ===========================================================

  e = 1.6021766208e-19;           % Elementary charge [C]
  h = 6.626070040e-34;            % Planck constant [ J.s]
  hbar = 1.054571800e-34;         % Reduced Planck's constant
  me = 9.10938356e-31;            % Electron mass [kg]

  s1 = 1e9;                       % m --> nm
  s2 = e;                          % J --> eV
 
   
% SETUP ===============================================================
  
  K = 2*me/hbar^2;
  E = s2*Es;                              % [J]
  xs = linspace(xMin,xMax,N)./s1;   % [m]

  
  % FREE PARTICLE =======================================================

% Wavelength  [n nm] and Propagation constant  [m^-1]
  wL = 2*pi/sqrt(K*E);              % [m]
  wLs = s1*wL;                      % [nm]
  k = 2*pi/wL;                      % [m^-1]
% Momentum  [kg.m.s^-1]
  p = hbar*k;                       %
% Angular frequency  [s^-1]
  w = E/hbar;
  
  
% POTENTIAL ENERGY  [eV]
  UMin = -100;
  UMax = 100;
  U = zeros(N,1);
  for c = 1:N  
    U(c) = pot(xs(c));
  end
    
   
  
% ode45
  xSpan = linspace(xMin,xMax,N)./s1;   % [m]
  dx = xSpan(2)-xSpan(1);
  xSpan = xMax/s1: -dx : xMin/s1;
%  xSpan = [xMin, xMax];
 u0 = [1,1i*k];
  options = odeset('RelTol',1e-6);
  
%  options1 = odeset('Refine',1);
%  options2 = odeset(options1,'NonNegative',1);
  
  [x,u] = ode45(@FNode, xSpan, u0, options);
  
  psi = -u(:,1)./max(u(:,1));
  toc

  

% GRAPHICS ==========================================================
  
  figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.02 0.05 0.4 0.5]);
    set(gcf,'color','w');
    hold on
    
    %subplot(2,1,1)
    subplot('position',[0.1 0.68 0.8 0.3])
    xP = s1.*xs; yP = U;
    xP = flip(xP); yP = flip(yP); 
    plot(xP,yP,'k','lineWidth',2)
    hold on
    yP = Es.*ones(N,1);
    plot(xP,yP,'r','lineWidth',2)
    grid on
    box on
    xlabel('x  [ nm ]')
    ylabel('E & U  [ eV ]')
    set(gca,'fontsize',12)
    
   % subplot(2,1,2)
   subplot('position',[0.1,0.1 0.8 0.48])
    xP = s1.*x; yP = -real(psi);
 %   xP = flip(xP); %yP = flip(yP); 
      plot(xP,yP,'b','lineWidth',2)
    hold on
    yP = imag(psi);
      plot(xP,yP,'r','lineWidth',2)
%     yP = abs(psi);
%       plot(xP,yP,'k','lineWidth',2)
%     yP = conj(psi).*psi;
%       plot(xP,yP,'g','lineWidth',2)
    grid on
    box on
    xlabel('x  [ nm ]')
    set(gca,'fontsize',12)
    ylabel('\psi','fontsize',20)
    
   
 
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
   U = 0;
    if x > 0.20*1e-9
        U = 200;
    end
   
%    if x > 0.3*1e-9
%        U = 0;
%    end
   
 end