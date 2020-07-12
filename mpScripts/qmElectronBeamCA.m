% qmElectronBeamCA.m

% QUANTUM MECHANICS 
% Electron Beams: Step Potential  TOTAL energy > step height (E > U0)
% Analytical Solution of Schrodinger Equation

% Calls external function: simpson1d.m

% Calculations use S.I. units
% conversion of units:  eV <---> J   nm <---> m

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/qmElectronBeamC.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b
% 200112

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
  xMin = -0.6;
  xMax = 1;
   
% Grid points and indices for region x < 0 and x > 0
  N = 5501;            % must be an odd number
   
% ANIMATION SETUP   ===================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif =0;
   ag_name = 'ag_qm.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 
  
  
% CONSTANTS ===========================================================
  e = 1.6021766208e-19;           % Elementary charge [C]
  h = 6.626070040e-34;            % Planck constant [ J.s]
  hbar = 1.054571800e-34;         % Reduced Planck's constant
  mE = 9.10938356e-31;            % Electron mass [kg]
  K = 2*mE/hbar^2;
  
% SETUP ===============================================================  
% nm <--> m
  X = 1e-9;
% Length of X domain  [m]
  L = X*(xMax - xMin);
% X grid   [m]
  x = X.*linspace(xMin,xMax,N);   
  dx = x(2) - x(1);
  
% Grid points and indices for region x < 0 and x > 0
  N0 = find(x>0,1);
  if mod(N0,2) == 0, N0 = N0-1; end
      
  N1 = 1:N0;
  N2 = N0:N;
  
% eV <---> J  
  E =  e*Es;
  U0 = e*U0s;

% Potential energy function: step potential E > U0
  U = zeros(N,1);
  U(x>0) = U0s;

  
% MODEL PARAMETERS  ==================================================  
% Region 1 and Region 2: first and second elements of a vector  

% Angular frequency  omega  [1/s]  
  w = E/hbar; 
% Frequency  [Hz] 
  f = w/(2*pi);
% Period T [s]
  Tp = 1/f;
% Time grid
  nT = 101;
  t = linspace(0,2*Tp,nT);
  psiT = exp(-1i*w*t);  
%   figure(9)
%   plot(t,psiT)
  

% Kinetic Energy  [eV and J]
  KEs(1) = Es;
  KEs(2) = Es - U0s;
  KE = e.*KEs;
  
% Momentum  [kg.m/s];
  p = sqrt(2*mE*KE);

% Velocity  [m/s]
  v = p./mE;

% Wavelength  lambda  [m and nm]
  wL = h./p;
  wLs = wL./X;

% Wave number  [1/m]
  k = p./hbar;

% Amplitudes: B reflected wave / C transmitted wave
  B = ( k(1)-k(2)) / (k(1)+k(2 ))*A;
  C = 2*k(1) /( k(1)+k(2) )*A;
  
% Reflection and transmission coefficients
  R = B^2/A^2;
  T = 1 - R;
  

% Eigenfunction  psi(x)
  psi = zeros(N,1);
  for c = 1 : N
    psi(c) = A*exp(1i*k(1)*x(c)) + B*exp(-1i*k(1)*x(c));
    if x(c) > 0
      psi(c) = C.*exp(1i*k(2)*x(c));
    end
  end

% Probability density  [1/m]
   psiS = conj(psi); 
   PD = conj(psi).*psi;
   psi_max = max(real(psi));
   
% Number of particles in X domain: x < 0 num1 and x > 0 num2 
  num(1) = simpson1d(PD(1:N2(1))',x(1),x(N2(1)));
  num(2) = simpson1d(PD(N2(1):N2(end))',x(N2(1)),x(N2(end)));
  
% Number density
  n(1) = num(1)/abs(x(1));
  n(2) = num(2)/abs(x(N));
    
% Probability current   [1/s]
  J = psiS.*gradient(psi,dx) - psi.*gradient(psiS,dx);
  J = -(1i*hbar/(2*mE)) .* J;
  Jnet = mean(J);
  
  Jinc   = v(1)*A^2;
  Jrefl  = -v(1)*B^2;
  Jtrans = v(2)*C^2;
 
  J1 = Jinc + Jrefl;
  J2 = Jtrans;
  

% GRAPHICS ==========================================================

figure(1)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.02 0.05 0.35 0.7]);
    set(gcf,'color','w');

    xP = x./X;
    
    for c = 1 : nT
        
    subplot(4,1,1)   
      yP = U;
      plot(xP,yP,'b','lineWidth',2)
      hold on
      yP = Es.*ones(N,1);
      plot(xP,yP,'r','lineWidth',2)
      grid on
      xlim([xMin xMax])
      set(gca,'xtick',xMin:0.2:xMax)
      ylabel('E   U  [eV]') 
      set(gca,'fontsize',12)     
      xtickformat('%2.1f')  
        
    subplot(4,1,2)  
      yP = real(psi .* psiT(c));
      yP = yP./psi_max;
      plot(xP,yP,'b','lineWidth',2)
      ylim([-1.2 1.2])
      xlim([xMin xMax])
      set(gca,'xtick',xMin:0.2:xMax)
      grid on 
      set(gca,'fontsize',12)
      ylabel('\psi_R','fontsize',16)
      xtickformat('%2.1f') 
      
    subplot(4,1,3)  
      yP = imag(psi .* psiT(c));
      yP = yP./psi_max;
      plot(xP,yP,'r','lineWidth',2)
      ylim([-1.2 1.2])
      xlim([xMin xMax])
      set(gca,'xtick',xMin:0.2:xMax)
      grid on 
      set(gca,'fontsize',12)
      ylabel('\psi_I','fontsize',16)  
      xtickformat('%2.1f') 
      
     subplot(4,1,4)  
      yP = PD;
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

    
%% Reflection  and Transmission coeffficients as functions of E/U0
figure(9)
    set(gcf,'units','normalized');
    set(gcf,'position',[0.38 0.05 0.25 0.3]);
    set(gcf,'color','w');
    
  % E/U0 < 0
    xP = linspace(0,1,201); yP1 = ones(1,201); yP2 = zeros(1,201);
    plot(xP,yP1,'b','linewidth',2)
    hold on
    plot(xP,yP2,'r','linewidth',2)
    
  % E/U0 > 0
    xP = linspace(1,4,201);
    yP1 = ( (1 - sqrt(1-1./xP)) ./ (1 + sqrt(1 - 1./xP)) ).^2;
    yP2 = 1- yP1;
    plot(xP,yP1,'b','linewidth',2)
    hold on
    plot(xP,yP2,'r','linewidth',2) 
    
    xP = [Es/U0s Es/U0s]; yP = [0 1];
    plot(xP,yP,'k','linewidth',0.5) 
    xP = [0 Es/U0s]; yP = [R R];
    plot(xP,yP,'b','linewidth',0.5) 
    xP = [0 Es/U0s]; yP = [T T];
    plot(xP,yP,'r','linewidth',0.5) 
    grid on
    box on
    xlabel('E / U_0')
    ylabel('R  and  T')
    legend('R','T')
    txt = sprintf('E/U_0 = %3.3f    R = %3.3f    T = %3.3f     ',Es/U0s,R,T);
    title(txt,'fontweight','normal');
    set(gca,'fontsize',12)
    
    
%%

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
   txt = sprintf('E = %3.2f eV    U_0 = %3.2f eV    E/U_0 = %3.2f', Es,U0s,Es/U0s);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
%    txt = sprintf('Step Height  U_0 = %3.0f  eV', U0s);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('\\omega = %3.2e s^{-1}         period = %3.2e s', w, Tp);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('A = %3.2e               A_N = %3.2e', A, An);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('B = %3.2e               B_N = %3.2e', B, Bn);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('C = %3.2e               C_N = %3.2e', C, Cn);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh; 
%    txt = sprintf('\\lambda_1 = %3.4f nm            \\lambda_{N1} = %3.4f nm', wL1/X,wL1_N/X);
%    text(0,hh,txt,'fontsize',12)
%   
%    hh = hh+dh; 
%    txt = sprintf('\\lambda_2 = %3.4f nm            \\lambda_{N2} = %3.4f nm', wL2/X,wL2_N/X);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('v_1 = %3.2e m.s^{-1}       v_2 = %3.2e m.s^{-1}', v1, v2);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('Electrons: x < 0  N_1 = %3.2e   x > 0  N_2 = %3.2e', num1, num2);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('Prob flux  J = %3.2e  s^{-1}', J_AC);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('Reflection coefficient  R = %3.2f ', R);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('Transmission coefficient  T = %3.2f ', T);
%    text(0,hh,txt,'fontsize',12)
%   
%    axis off
   
toc


 
 
