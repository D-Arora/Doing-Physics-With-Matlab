% math_ft_kx1.m

% Fourier Transform by direct integration
%  Variables x and k

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/maths_ft01.pdf
% Sripts: Download Location
%    http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

clear 
close all
clc

tic
% INPUTS / FUNCTION h(x) ==============================================
% PositionTime: sample points (must be an odd integer)  [2001]
% Wavefunction: sample points (must be an odd integer)   [2001]
   nX = 2001;
   nK = 2001;

% Gaussian Function h(x)   
  xMin = -3;
  xMax = 3;
  kMin = -3;
  kMax = 3;
  s = 1;
    
  x = linspace(xMin,xMax,nX);
  k = linspace(kMin,kMax,nK);
  
  h = exp(-0.5*((x)/s).^2);
  
  H = zeros(1,nK);
  hI = zeros(1,nX);
  HT = zeros(1,nK);

% FOURIER TRANSFORM
for c = 1:nK
   g = h.* exp(-1i*k(c)*x);
   H(c) = simpson1d(g,xMin,xMax);
end
   H = (1/sqrt(2*pi)).*H;
   P = conj(H).*H;
   
% FOURIER SERIES
  NK = 9;
  kn = linspace(kMin,kMax,NK);
  cn = zeros(NK,1);
  for c = 1:NK
    u = (1/sqrt(2*xMax)) .*exp(1i*kn(c)*x);
    cn(c) = simpson1d(u.*h,xMin,xMax);
  end
   Pc = conj(cn).*cn;
   Pc = Pc./max(Pc);
   cn_avg = sqrt(sum(Pc));
   
   f = zeros(nX,1);
  for c = 1:NK
   f = f + cn(c).*exp(1i*kn(c)*x)';
  end
   f = (1/sqrt(2*xMax)).*f;
 

% Generate the plot, title and labels.
figure(1);
   pos = [0.05 0.05 0.30 0.40];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

subplot(3,1,1)
   xP = x; yP = h;
   plot(xP,yP,'lineWidth',2);
  % title('Function h(x)')
   hold on
   plot(xP(1:100:end),real(f(1:100:end)./max(f)),'ro')
 %  plot(xP,abs(f))
   xlabel('position x [m]')
   ylabel('h(x)')
   grid on
   set(gca,'fontsize',12)

subplot(3,1,2)
   xP = k; yP = real(H);
   plot(xP,yP,'lineWidth',2);
   hold on
   xP = k; yP = imag(H);
   plot(xP,yP,'lineWidth',1)
   %title('Function H(y)')
   xlabel('wavenumber k [m^{-1}]')
   ylabel('H(k)')
   grid on
   set(gca,'fontsize',12)
   
 subplot(3,1,3)  
   xP = k; yP = P;
   plot(xP,yP,'lineWidth',2);
  % title('Function h(x)')
   xlabel('wavenumber k [m^{-1}]')
   ylabel('H^* H')
   grid on
   set(gca,'fontsize',12)

figure(2);
   pos = [0.37 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = kn; yP = abs(cn);
   bar(xP,yP,0.2)
   
   hold on
   xP = x; yP = abs(H);
   plot(xP,yP,'b','lineWidth',2);
     
   xlabel('k_n  [m^{-1}]');
   ylabel('c_n  &  H');
   grid on
   set(gca,'fontsize',12)
    
   
figure(3);
   pos = [0.65 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   hold on
   z= zeros(nX,1);
   for c = 1:9
   xP = x; yP = (1/sqrt(2*xMax)).*cn(c).*exp(1i*kn(c)*x)';
   z = z+yP;
   plot(xP,real(yP),'lineWidth',2)
   end
   plot(xP,real(z))
   title('Double Sided FFT (no shift)');
   xlabel('Normalized frequency');
   ylabel('power P');
   grid on   

