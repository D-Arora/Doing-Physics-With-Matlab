% QMG2BB.m


% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG02C.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230426   Matlab R2021b

clear; close all;clc

% SETUP  ===========================================================
  a = 1e-9;
  NX = 159;      % [201]   must be an odd number - number of grid points
  NT = 68000;   % [5e5]  number of time steps

  me = 9.10938291e-31;    % electron mass
  hbar = 1.054571726e-34; % hbar Planck's constant
  h = 6.62607015e-34; 
  e = 1.602176565e-19;    % elementary charge
  
  x1 = 0; x2 = a;
  x = linspace(x1,x2,NX); 
  dx = x(2) - x(1);

  y = zeros(NT,NX);
  yR = zeros(NT,NX);
  yI = zeros(NT,NX);
  y(1,:) = x.*(a-x);

  fn = y(1,:).^2;
  Area = simpson1d(fn,x1,x2);
  y(1,:) = y(1,:)./sqrt(Area);
  fn = y(1,:).^2;
  check = simpson1d(fn,x1,x2);
  
  C = 1/2;
  dt = C*2*me*dx^2/hbar;
  t = 0:dt:NT*dt;
   
  yR(1,:) = y(1,:);

% SOLVE SCHRODINER EQUATION ===========================================
for nt = 2 : NT
for nx = 2:NX-1    
      yR(nt,nx) = yR(nt-1,nx) - C*( yI(nt-1,nx+1) - 2*yI(nt-1,nx) + yI(nt-1,nx-1) ) ;
end
for nx = 2:NX-1          
      yI(nt,nx) = yI(nt-1,nx) + C*( yR(nt,nx+1) - 2*yR(nt,nx) + yR(nt,nx-1) ) ;
end      
end
 
% Probabilitiy density
   probDensity = conj(yR+1i.*yI).*(yR+1i.*yI);
 
% FOURIER TRANSFORM CALCULATIONS ========================================
   nF = 299;
   H = zeros(1,nF); 
   FT= yR(1:999:NT,100)';
   f = linspace(0,5e14,nF);
   tF = t(1:999:NT);
% Fourier Transform  H(f)
   for c = 1:nF
     g = FT.* exp(1i*2*pi*f(c)*tF);
     H(c) = simpson1d(g,0,max(t));
   end
     psd = conj(H).*H;

% EXPECTATION VALUE  kinetic energy = total energy E ==================
  psi = yR(1,:);
  psi1 = gradient(psi,dx);
  psi2 = gradient(psi1,dx);
  probD = psi.*psi2;
  fn = probD;
  int1 = simpson1d(fn,x1,x2);
  E = -(hbar^2/(2*me))*int1/e;
  ind = find(psd == max(psd));
  fpeak = f(ind(1));
  Ef = h*fpeak/e;
  ET = (h/(2*a))^2/(2*me*e);

% EXPECTATION VALUES & UNCERTAINTIES  
  psi = yR(1,:) + 1i.*yI(1,:);
  fn = conj(psi).*x.*psi;
  xAVG = simpson1d(fn,x1,x2);
  fn = conj(psi).*x.^2.*psi;
  x2AVG = simpson1d(fn,x1,x2);
  xSTD = sqrt(x2AVG - xAVG^2);

  psi1 =  gradient(psi,dx);
  psi2 = gradient(psi1,dx);
  fn = conj(psi).*psi;
  pAVG = -1i*hbar*simpson1d(fn,x1,x2);
  fn = conj(psi).*psi2;
  p2AVG = -hbar^2*simpson1d(fn,x1,x2);
  pSTD = sqrt(p2AVG - conj(pAVG)*pAVG);
  dxdp = xSTD*pSTD;

% OUTPUT  ===========================================================
  disp('Fourier transform')
  fprintf('   fpeak = %2.2e Hz  \n', fpeak)
  fprintf('   period = %2.2f fs   \n', 1e15 / fpeak)
  disp('  ')
  disp('Total energy E')
  fprintf('   expectational value, <E> = %2.2f eV \n', E)
  fprintf('   Fourier transform, Ef = h fpeak = %2.2f eV \n', Ef)
  fprintf('   Theory, ET = p^2/2m = %2.2f eV \n', ET)
  disp('  ')
  disp('Uncertainty Principle')
  fprintf('   dx = xSTD = %2.2e m \n', xSTD)
  fprintf('   dp = pSTD = %2.2e m \n', pSTD)
  fprintf('   dx dp = %2.2e N.s \n', dxdp)
  fprintf('   hbar/2 = %2.2e N.s.m \n', hbar/2)
  disp('   dx dp > hbar/2 ')


% ANIMATION SETUP =======================================================
% (0 no)  (1 yes) for flag1
   flag1 = 1;    
% file name for animated gif   
    ag_name = 'agQMG2BB.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.00; 
% Frame to start
    frame1 = 0;  


% GRAPHICS =========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.05 0.2 0.4]);
  set(gcf,'color','w');
  FS = 14;

  for c = 1: 1000:NT
subplot(2,1,1)  
   xP = x; yP = yR(c,:);
  plot(xP./1e-9,yP,'b','LineWidth',2)
  grid on
  xlabel('x  [ nm ]')
  ylabel('real(\psi)')
 % title('Probability density','FontWeight','normal')
  set(gca,'FontSize',FS)
  ylim([-8e4 8e4])
  xticks(0:0.25:1)

subplot(2,1,2)
  xP = x; yP = probDensity(c,:);
  plot(xP./1e-9,yP,'b','LineWidth',2)
  grid on
  xlabel('x  [ nm ]')
  ylabel('|\psi|^2')
%  title('Probability density','FontWeight','normal')
  set(gca,'FontSize',FS)
  ylim([-4e9 4e9])
  xticks(0:0.25:1)

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
   pause(0.1)
   end
 
figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.2 0.2]);
  set(gcf,'color','w');
  FS = 14;
  xP = t(1:end-1)./1e-15; yP = yR(:,100)';
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlabel('t  [ fs ]')
  ylabel('\psi(x= 0.5)')
  title('Wavefunction','FontWeight','normal')
  ylim([-5e4 5e4])  
 % axis tight
  set(gca,'FontSize',FS)

figure(3)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.65 0.05 0.2 0.2]);
  set(gcf,'color','w');
  FS = 14;
  xP = f; yP = 10.*log10((conj(H).*H)./max(conj(H).*H));
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlim('tight')
  xticks((0:1:5).*1e14)
  xlabel('f  [ Hz ]')
  ylabel('psd  [ dB ]')
  title('Fourier Transform','FontWeight','normal')
  set(gca,'FontSize',FS)
