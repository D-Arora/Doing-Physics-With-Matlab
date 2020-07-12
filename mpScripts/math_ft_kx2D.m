% math_ft_kx2D.m

% Fourier Transform by direct integration
%  Variables x and k
%  Fourier Transform of a wave packet


% Ian Cooper
% email: matlabvisualphysics@gmail.com
%
% 200427 / Matlab version R2020aa

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/qm2DA.htm
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
  s = 1e-3; k0 = 100; 
  xMin = -0.2; xMax = 0.2;
  kMin = -50; kMax = 250;
  
  x = linspace(xMin,xMax,nX);
  k = linspace(kMin,kMax,nK);
  
  h = 10 .* exp(-x.^2./s) .* exp((1i*k0).*x);
    
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
   
   
% Generate the plot, title and labels.
figure(1);
   pos = [0.05 0.05 0.30 0.40];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

subplot(3,1,1)
   xP = x; yP = real(h);
   plot(xP,yP,'lineWidth',2);
   title('Wave Packet / Fourier Transform')
  
   xlabel('position x')
   ylabel('h(x)')
   grid on
   set(gca,'fontsize',12)

subplot(3,1,2)
   xP = k; yP = real(H);
   plot(xP,yP,'lineWidth',2);
   hold on
   xP = k; yP = imag(H);
   plot(xP,yP,'lineWidth',1)
   xlabel('wavenumber k')
   ylabel('H(k)')
   grid on
   set(gca,'fontsize',12)
   
 subplot(3,1,3)  
   xP = k; yP = P;
   plot(xP,yP,'lineWidth',2);
   xlabel('wavenumber k')
   ylabel('H^* H(k)')
   grid on
   set(gca,'fontsize',12)


   


