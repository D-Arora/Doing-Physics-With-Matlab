% QMG2B.m


% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG?.pdf
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230426   Matlab R2021b

clear; close all;clc

% SETUP  ===========================================================
  a = 1e-9;
  NX = 201;      % [201]   must be an odd number - number of grid points
  NT = 5e5;   % [5e5]  number of time steps

  me = 9.10938291e-31;    % electron mass
  hbar = 1.054571726e-34; % hbar Planck's constant
  h = 6.62607015e-34; 
  e = 1.602176565e-19;    % elementary charge
  
  x1 = 0; x2 = a;
  x = linspace(x1,x2,NX); dx = x(2)-x(1);
  C = 1/10;
  dt = C * 2 * me * dx^2 / hbar;
  t = 0:dt:(NT-1)*dt;

  y = zeros(NT,NX);
  yR = zeros(NT,NX);
  yI = zeros(NT,NX);
  y(1,:) = x.*(a-x);

  fn = y(1,:).^2;
  Area = simpson1d(fn,x1,x2);
  y(1,:) = y(1,:)./sqrt(Area);
  fn = y(1,:).^2;
  check = simpson1d(fn,x1,x2);
  
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
   f = linspace(0,10e14,nF);
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
  f0 = f(ind(1));
  Ef = h*f0/e;
  ET = (h/(2*a))^2/(2*me*e);

% OUTPUT  ===========================================================
  disp('Expectation value: total energy eV')
  fprintf('total energy, E = %2.2f eV \n', E)
  fprintf('total energy, Ef = hf0 = %2.2f eV \n', Ef)
  fprintf('total energy, ET = (h/wL)/(2m) = %2.2f eV \n', ET)


% ANIMATION SETUP =======================================================
% (0 no)  (1 yes) for flag1
   flag1 = 1;    
% file name for animated gif   
    ag_name = 'agQMG2B.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.01; 
% Frame to start
    frame1 = 0;  
% GRAPHICS =========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.05 0.25 0.65]);
  set(gcf,'color','w');
  FS = 14;

  for c = 1:round(NT/500):NT
subplot(2,1,1)
  xP = x; yP = yR(c,:);
  plot(xP./1e-9,yP,'b','LineWidth',2)
  grid on
  xlabel('x  [ nm ]')
  ylabel('\psi')
  title('Wavefunction','FontWeight','normal')
  set(gca,'FontSize',FS)
  ylim([-5e4 5e4])
  yticks(-5e4:1e4:5e4)
subplot(2,1,2)
  xP = x; yP = probDensity(c,:);
  plot(xP./1e-9,yP,'b','LineWidth',2)
  grid on
  xlabel('x  [ nm ]')
  ylabel('|\psi|^2')
  title('Probability density','FontWeight','normal')
  set(gca,'FontSize',FS)
  ylim([0 2.5e9])
 % yticks(-5e4:1e4:5e4)
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
  pause(0.001)
  end

figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.35 0.05 0.25 0.35]);
  set(gcf,'color','w');
  FS = 14;
  xP = t./1e-15; yP = yR(:,100);
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlabel('t  [ fs ]')
  ylabel('\psi(x= 0.5)')
  title('Wavefunction','FontWeight','normal')
  axis tight
  set(gca,'FontSize',FS)

figure(3)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.65 0.05 0.25 0.35]);
  set(gcf,'color','w');
  FS = 14;
  xP = f; yP = 10.*log10((conj(H).*H)./max(conj(H).*H));
  plot(xP,yP,'b','LineWidth',2)
  grid on
  xlabel('f  [ Hz ]')
  ylabel('psd  [ dB ]')
  title('Fourier Transform','FontWeight','normal')
  set(gca,'FontSize',FS)
