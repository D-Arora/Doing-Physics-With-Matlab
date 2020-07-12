% qm001.m



clear
close all
clc


tic

% INPUTS ==============================================================

% Energy of particles  [eV];
  Es = 250;
% X axis  [nm]
  xMin = 0;
  xMax = 0.5;
% Grid points 
  N = 5501;

  % CONSTANTS ===========================================================

  e = 1.6021766208e-19;           % Elementary charge [C]
  h = 6.626070040e-34;            % Planck constant [ J.s]
  hbar = 1.054571800e-34;         % Reduced Planck's constant
  me = 9.10938356e-31;            % Electron mass [kg]

  s1 = 1e-9;                       % m --> nm
  s2 = e;                          % J --> eV
 
   
% SETUP ===============================================================
  
  K = 2*me/hbar^2;
  E = s2*Es;                    % [J]
  xs = linspace(xMin,xMax,N);   % [nm]
  x  = s1.*xs;                  % [m]
  dx = x(2) - x(1);               % [m]
  
  
% FREE PARTICLE =======================================================

% Wavelength  [n nm] and Propagation constant  [m^-1]
  wL = 2*pi/sqrt(K*E);              % [m]
  wLs = s1*wL;                      % [nm]
  k = 2*pi/wL;                      % [m^-1]
% Momentum  [kg.m.s^-1]
  p = hbar*k;                       %
% Angular frequency  [s^-1]
  w = E/hbar;
  
  
% POTENTIAL ENERGY  [eV  J]
  UMin = -100;
  UMax = 100;
  Us = zeros(N,1);
  for c = 1:N  
    Us(c) = 0;
    if xs(c) > 0.2
        Us(c) = 200;
    end
   
%     if xs(c) > 0.3
%         Us(c) = 0;
%     end 
%     
  end
  
  U = s2.*Us; 
  
% FINITE DIFFERENCE METHOD: Wave Function psi =========================
  Kse = 2 + (dx^2*K).*(U - E);
  
  u = zeros(N,1);
  u(N) = 1; u(N-1) = 0;
    
  for c = N-1 : -1: 2
   u(c-1) = Kse(c) * u(c) - u(c+1);
  end
  
  psi = u./max(u);
  
  


  

% GRAPHICS ==========================================================
  
  figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.02 0.05 0.4 0.5]);
    set(gcf,'color','w');
    hold on
    
    %subplot(2,1,1)
    subplot('position',[0.1 0.68 0.8 0.3])
    xP = xs; yP = Us;
    plot(xP,yP,'k','lineWidth',2)
     hold on
    yP = Es.*ones(N,1);
    plot(xP,yP,'r','lineWidth',2)
    grid on
    box on
    xlabel('x  [ nm ]')
    ylabel('U  [ eV ]')
    set(gca,'fontsize',12)
    
   % subplot(2,1,2)
   subplot('position',[0.1,0.1 0.8 0.48])
    xP = xs; yP = psi;
       plot(xP,yP,'b','lineWidth',2)
%     hold on
%      yP = psi.*psi;
%        plot(xP,yP,'r','lineWidth',2)
%     yP = abs(psi);
%       plot(xP,yP,'k','lineWidth',2)
% %     yP = conj(psi).*psi;
%       plot(xP,yP,'g','lineWidth',2)
    grid on
    box on
    xlabel('x  [ nm ]')
    set(gca,'fontsize',12)
    ylabel('\psi','fontsize',18)
    
   figure(2)
   xP = xs; yP = psi.*psi;
    plot(xP,yP,'k','lineWidth',2)
