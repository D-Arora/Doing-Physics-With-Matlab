% emElectronBeamA.m

% QUANTUM MECHANICS 
% Electron Beam in a Field-Free Space

% Calls external function: simpson1d.m

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/qmElectronBeamA.htm
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
% Energy of electron beam  [400 eV]
  Es = 400;

% Amplitudes:incident beam (left to right) A / reflected beam (right to left) B
% Defaults: A = 10^8   B_A = 0  or  B_A = 0.5
  A = 10^8;
% Ratio B / A  
  B_A = 0.5;
  
  B = B_A*A;
% Wave right to left A = 0 and B = A (uncomment / comment)
%   B = A; A = 0; 
  
% Define X axis  [L = xMax   0.5e-9 m]
  xMin = 0;
  xMax = 0.5e-9;
% Number of grid pints: must be an odd number  [999]  
  nX = 999;
  
% Number of time steps for simulation [301]
  nT = 301;
% Number of periods   [6]  
  num = 6;
  
  
% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_qm.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 
  
  
% SETUP ==============================================================
% Fundamental Charge  [C]
  e = 1.60218e-19;
% Electron mass [kg]
  mE = 9.10939e-31; 
% Planck constant [J.s]
  h = 6.626070040e-34; 
% Reduced Planck constant [J.s]
  hbar = 1.054571800e-34; 
  
% X grid  [m]
  x = linspace(xMin,xMax,nX);
  dx = x(2) - x(1);
  
  
% ELECTRON BEAM PARAMETERS ============================================  
% Energy of electron beam  [J]
  E = e*Es;
% Angular frequency omega  [1/s]
  w = E/hbar;
% Wave number  [1/m]
  k = sqrt(2*mE*E)/hbar;
% Wavelength  lambda  [m]
  wL = 2*pi/k;
% Momentum  [kg.m/s]
  p = hbar*k;
% Velocity  [m/s]
  v = p/mE;
% Period  T  [s]
  T = 1/w;
% frequency  [Hz]
  f = 1/T;  
  
% Simulation time   [s]
  t = linspace(0,num*T,nT);

  
% WAVE FUNCTION =======================================================
  PSI_A = zeros(nT,nX);       % --->
  PSI_B = zeros(nT,nX);       % <---
  PSI   = zeros(nT,nX);       % ---> + <---
  PSIs  = zeros(nT,nX);       % complex conjugate
  J     = zeros(nT,nX);       % probabilty current  [s^-1]

for c = 1:nT
    PSI_A(c,:) = A.*exp(1i.*( k.*x - w*t(c)));
    PSI_B(c,:) = B.*exp(1i.*(-k.*x - w*t(c)));
    PSI(c,:)   = PSI_A(c,:) + PSI_B(c,:);
    PSIs(c,:) = conj(PSI(c,:));
    J(c,:) = PSIs(c,:).*gradient(PSI(c,:),dx) - PSI(c,:).*gradient(PSIs(c,:),dx);
end
 
% Probability Density  [m^-1]   
  ProbDensity = conj(PSI).*PSI;
  
% Number in electrons in length L  
  N = simpson1d(ProbDensity(1,:),0,xMax); 
% Number density  [1/m
  n = N/xMax;
  
% Probability current  [s^-1]    
  J = -(1i*hbar/(2*mE)) .* J;
% max value  
  Jmax = max(max(J));
  
% Amplitude numerical estimates:  A ---> and B  <---
  Amax = max(max((abs(PSI))));
  Amin = min(max((abs(PSI))));
  An = 0.5*(Amax + Amin);
  Bn = 0.5*(Amax - Amin);
  

% GRAPHICS ============================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.30 0.65]);
  set(gcf,'color','w');

  xP = 1e9.*x;    %  [nm]
  
for c = 1:5:nT
    
  subplot(4,1,1)  
    yP = real(PSI_A(c,:));
    plot(xP,yP,'b','linewidth',2)
    hold on
    yP = real(PSI_B(c,:));
    plot(xP,yP,'r','linewidth',2)
    yP = real(PSI_A(c,:)+ PSI_B(c,:)) ;
    plot(xP,yP,'k','linewidth',1)
    
    hold off
    set(gca,'fontsize',12)
    ylabel('real(\psi  \psi_A  \psi_B ) ','fontsize',16)
    grid on
    ylim(max(max(abs(PSI))).*[-1 1])
  
    
  subplot(4,1,2)
    yP = real(PSI(c,:));
    plot(xP,yP,'b','linewidth',2)
    set(gca,'fontsize',12)
    ylabel('real(\psi)','fontsize',16)
    grid on
    ylim(max(max(abs(PSI))).*[-1 1])
    
  subplot(4,1,3)
    yP = imag(PSI(c,:));
    plot(xP,yP,'r','linewidth',2)
    set(gca,'fontsize',12)
    grid on
    ylabel('imag(\psi)','fontsize',16)
    grid on
    ylim(max(max(abs(PSI))).*[-1 1])
  
  if A ~= B  
  subplot(4,1,4)
    
    yP = J(c,:);
    plot(xP,yP,'k','linewidth',2)
    hold on
    grid on
    xlabel('x  [ nm ] ')
    ylabel('J   [ s^{-1} ]')
    set(gca,'fontsize',12)
    grid on
    
    if max(max(J)) > 0
      ylim([0 1.1*max(max(J))])
    else
     ylim([1.1*min(min(J)) 0])
    end
   end    
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
    
   pause(0.002)
end


figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.33 0.05 0.25 0.25]);
   set(gcf,'color','w');

   xP = 1e9.*x;
   
   yP = ProbDensity(1,:);
 %  plot(xP,yP,'b','linewidth',2)
   area(xP,yP)
   xlabel('x  [ m ]')
   ylabel('\Psi^* \Psi  [ m^{-1} ]','fontsize', 16)
   set(gca,'fontsize',12)
   grid on
 
   ytickformat('%3.2f')
    
 %%   

% DISPLAY RESULTS -----------------------------------------------------

figure(7)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.60 0.05 0.25 0.75]);
   set(gcf,'color','w');

   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -18;
   
   txt = 'ELECTRON BEAM IN A FIELD-FREE SPACE';
   Htext = text(0,hh,txt,'fontsize',14);
   set(Htext,'color','b')
   hh = hh+dh;
   txt = sprintf('Electron energy  E = %3.0f  eV', Es);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Electron momentum  p = %3.2e  kg.m.s^{-1}', p);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Electron velocity  v = %3.2e m.s^{-1}', v);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('deBroglie wavelength  \\lambda = %3.2e m', wL);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Wave number  k = %3.2e m^{-1}', k);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Angular frequency  \\omega = %3.2e s^{-1}', w);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('frequency  f = %3.2e Hz', f);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Period  T = %3.2e s', T);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Amplitudes:');
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('   B/A = %3.2f     A = %3.2e     B = %3.2e', B_A,A,B);
   text(2,hh,txt,'fontsize',12)
   
   hh = hh+dh; 
   txt = sprintf('   A_n = %3.2e     B_n = %3.2e', An,Bn);
   text(2,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('Number of electrons: N = %3.2e   n = %3.2e m^{-1} ', N, n);
   text(0,hh,txt,'fontsize',12)
   
   if A == B; Jmax = 0; end
   hh = hh+dh;
   txt = sprintf('Probability current J_{max} = %3.2e s^{-1}', Jmax);
   text(0,hh,txt,'fontsize',12)
   
   
   
   
   axis off
   
   
   

  
  
  