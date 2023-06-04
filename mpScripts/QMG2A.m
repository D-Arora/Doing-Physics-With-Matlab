% QMG2A.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG02.pdf
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230422   Matlab R2021b

% GRIFFITH Problem 2.1

 clear; close all;clc

% =======================================================================   
% Setup for saving images animated gif file and
% =======================================================================
   flagA =0;
   ag_name = 'ag_QMG21C.gif';
   delay = 0.1;  %  Delay in seconds before displaying the next image  
   nt = 1;%  Frame counter start
   frame1 = 0;

% SETUP =================================================================
% Input quantum numbers
  n(1) = 3; n(2) = 4;
% Constants
  h = 6.62607015e-34;
  hbar = 1.05457182e-34;
  m = 9.1093837e-31;
  e = 1.60217663e-19;
  L  = 5e-10;
  cS = sqrt(2/L);
  c1 = 3; c2 = 1;

  k0 = pi/L;
  NX = 299;
  NT = 169;
  tmax = 1e-15;
  x = linspace(-L/2,L/2,NX);
  t = linspace(0,tmax,NT);
% Wavefunctions: stationary states
  PSI1 = zeros(NX,NT);
  PSI2 = zeros(NX,NT);
  E    = zeros(2,1);      % total energy
  w    = zeros(2,1);      % angular frequency
  wL   = zeros(2,1);      % wavelength
  psi1 = cS.*sin(n(1)*k0*(L/2-x));
  psi2 = cS.*sin(n(2)*k0*(L/2-x));
  for ct = 1:2
    [E(ct),w(ct),wL(ct)] = energy(n(ct),hbar,m,L);
  end
% Check normalization
  check(1) = simpson1d(psi1.^2,-L/2,L/2);
  check(2) = simpson1d(psi2.^2,-L/2,L/2);
% Probabilty density
  pD1 = psi1.^2; pD2 = psi2.^2;
% Time dependence
  for ct = 1:NT
      PSI1(:,ct) = psi1.*exp(-1i*w(1)*t(ct));
      PSI2(:,ct) = psi2.*exp(-1i*w(2)*t(ct));
  end
% WAVEFUNCTION: linear combination
  PSI = c1.*PSI1 + c2.*PSI2;
% Normalize function
  A = simpson1d((conj(PSI(:,1)).*PSI(:,1))',-L/2,L/2);
  AN = 1/sqrt(A);
  PSI = PSI./sqrt(A);
  CHECK = simpson1d((conj(PSI(:,1)).*PSI(:,1))',-L/2,L/2);
  pD = conj(PSI).*PSI;
  c1 = c1*AN; c2 = c1*AN;

% EXPECTATION VALUE position <x>
  xavg = zeros(1,NT);
  for ct = 1: NT 
    fn = x'.*pD(:,ct);
    xavg(ct) = simpson1d(fn',-L/2,L/2);
  end

% Compound state - oscillation of prob. desnity
  W = (E(2) - E(1))/hbar;   % angular frequency
  T = (2*pi/W)/1e-15;               % period
  [pks,locs] = findpeaks(xavg);
  Tp = (t(locs(2)) - t(locs(1)))/1e-15;

% OUTPUT
   disp('Stationary states')
   fprintf('n = %2.0f   E = %2.4f eV  \n', n(1), E(1)/e)
   fprintf('n = %2.0f   E = %2.4f eV  \n', n(2), E(2)/e)
   disp('Compound state: period of oscillation')
   fprintf('Simulation: T = %2.3f fs  \n', Tp)
   fprintf('Theory:     T = %2.3f fs  \n', T)
   disp('Check normalization')
   fprintf('#1 --> %2.2f   #2 --> %2.2f eV  \n', check(1), check(2))
   fprintf('c1*psi1 + c2*psi2  --> %2.2f  \n', CHECK)
   fprintf('c1 = %2.2f   c2 = %2.2f  \n', c1, c2)
   fprintf('c1^2 + c2^2 = %2.2f  \n', c1^2 + c2^2)
   disp('Probability for total energy E measurement')
   fprintf('State #1, prob(E) = %2.2f   State #2, prob(E) = %2.2f  \n', c1^2, c2^2)

% GRAPHICS ==========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.40 0.60]);
  set(gcf,'color','w');
  FS = 14;

for ct = 1:NT  

subplot(3,2,1)
  xP = x.*1e9; yP = real(PSI1(:,ct));
  plot(xP,yP,'b','LineWidth',2)
  grid on
  ylim([-10 10].*1e4)
  xlim([-L/2 L/2].*1e9)
%  xlabel('x  [nm]')
  ylabel('\Psi_1')
  txt = sprintf('\\Psi_1(x,t)   n_1 = %2.0f \n',n(1));
  title(txt,'FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(3,2,2)
  xP = x.*1e9; yP = pD1;
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlim([-L/2 L/2].*1e9)
%  xlabel('x  [nm]')
  ylabel('|\Psi_1|^2')
  title('|\Psi_1(x)|^2','FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(3,2,3)
  xP = x.*1e9; yP = real(PSI2(:,ct));
  plot(xP,yP,'b','LineWidth',2)
  xlim([-L/2 L/2].*1e9)
  ylim([-10 10].*1e4)
  grid on
%  xlabel('x  [nm]')
  ylabel('\Psi_2')
  txt = sprintf('\\Psi_2(x,t)   n_2 = %2.0f \n',n(2));
  title(txt,'FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(3,2,4)
  xP = x.*1e9; yP = pD2;
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlim([-L/2 L/2].*1e9)
 % xlabel('x  [nm]')
  ylabel('|\Psi_2|^2')
  title('|\Psi_2(x)|^2','FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(3,2,5)
  xP = x.*1e9; yP = real(PSI(:,ct));
  plot(xP,yP,'b','LineWidth',2)
  ylim([-10 10].*1e4)
  grid on
  xlim([-L/2 L/2].*1e9)
  xlabel('x  [nm]')
  ylabel('\Psi')
  title('\Psi(x,t)','FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(3,2,6)
  xP = x.*1e9; yP = pD(:,ct);
  plot(xP,yP,'b','LineWidth',2)
  xlim([-L/2 L/2].*1e9)
  ylim([0 10].*1e9)
  grid on
  xlabel('x  [nm]')
  ylabel('|\Psi|^2','FontWeight','normal')
  title('|\Psi(x,t)|^2','FontWeight','normal')
  set(gca,'FontSize',FS)  
  pause(0.001)

   if flagA == 1
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
  set(gcf,'position',[0.35 0.05 0.2 0.2]);
  set(gcf,'color','w');
  FS = 14;
  plot(t./1e-15,xavg./1e-9,'b','LineWidth',2)
  grid on; box on
  axis tight
  xlabel('t  [ fs ]')
  ylabel('<x>  [ nm ]','FontWeight','normal')
  title('<x>','FontWeight','normal')
  set(gca,'FontSize',FS)


% FUNCTIONS ========================================================

function  [E,w,wL] = energy(n,hbar,m,L)
     E = n^2*pi^2*hbar^2/(2*m*L^2);
     w = E/hbar;
     wL = n*(L/2);
end
