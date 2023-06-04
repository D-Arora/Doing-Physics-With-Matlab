% maths_ft_04.m

% SPECTROGRAMS 

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 230407 / Matlab version R2021ba

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/?.pdf


% https://towardsdatascience.com/brief-introduction-of-hamming-and-hanning-function-as-the-preprocessing-of-discrete-fourier-8b87fe538bb7

% https://www.egr.msu.edu/classes/me451/me451_labs/Fall_2013/Understanding_FFT_Windows.pdf


clear
close all
clc

% SETUP =======================================================
% Time domain
  nT = 999; tMin = 0; tMax = 0.1;
  t = linspace(tMin,tMax,nT);
  f0 = 100;
  h = sin(2*pi*f0*t);
  w = hamming(nT)';

% Freqeuncy domain
  nF = 2999; fMax = 200; fMin = -fMax;
  f = linspace(fMin,fMax, nF);

  if mod(nT,2) == 0; nT = nT-1; end
  h = h(1:nT);
  [n,m] = size(h);
  if n > 1; h = h'; end 
 % h = h - mean(h);
 
% Fourier transform
  H = fourier(h,t,f);

% Inverse Fourier transform  
  hI = IFT(H,t,f);

% Hamming window
  hw = h.*w;
  Hw = fourier(hw,t,f);
  
  hIw = IFT(Hw,t,f);

% GRAPHICS  ========================================================
  fs = 14;

figure(1)   
   pos = [0.02 0.05 0.27 0.67];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

subplot(2,1,1)
   xP = t; yP = h;
   plot(xP,yP,'b','lineWidth',2);
   title('signal  h(t)')
   xlabel('t')
   ylabel('h(t)')
   grid on
   set(gca,'fontsize',fs)

subplot(2,1,2)
   xP = t; yP = real(hw);
   plot(xP,yP,'b','lineWidth',2);
   title('Hamming window h_{w}(t)')
   xlabel('t')
   ylabel('hI(t)')
   grid on
   set(gca,'fontsize',fs)

% subplot(3,1,3)
%    xP = t; yP = real(hIw);
%    plot(xP,yP,'rs','lineWidth',2);
%    hold on
%    yP = hw;
%    plot(xP,yP,'bo','lineWidth',2);
%    title('Hamming: IFT  h_{inv}(t)')
%    xlabel('t')
%    ylabel('hI(t)')
%    grid on
%    set(gca,'fontsize',fs)   
  
figure(2)   
   pos = [0.32 0.05 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = f; yP = abs(H);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   yP = abs(Hw);
   plot(xP,yP,'r','lineWidth',2);
   axis tight
   xlim([0 fMax])
   title('F.T.  H(f)')
   xlabel('t')
   ylabel('h(t)')
   grid on
   set(gca,'fontsize',fs)

% FUNCTIONS ===========================================================

function H = fourier(h,t,f)  % Fourier transform H(f)
   nT = length(t); tMin = min(t); tMax = max(t);
   nF = length(f);
   H = zeros(1,nF);
   for c = 1:nF
       g = h.* exp(1i*2*pi*f(c)*t);
       H(c) = simpson1d(g,tMin,tMax);
    end
end

function hI = IFT(H,t,f)  % inverse Fourier transform hI(t)
   nT = length(t); fMin = min(f); fMax = max(f);
   hI = zeros(1,nT);
   for c = 1:nT
       g = H.* exp(-1i*2*pi*f*t(c));
       hI(c) = simpson1d(g,fMin,fMax);
    end
end
