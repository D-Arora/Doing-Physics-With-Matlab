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

clear
close all
clc



% IMPORT SOUND
load splat

signal = y;

% SPECTROGRAM
  sound(y, Fs)

  % Time domain   
  Ls = length(signal);
  dt = 1/Fs;
  ts = linspace(0,Ls*dt,Ls);
       
% Function h(t) 
  txt_h = 'h(t)';

%   index1 = round(Ls/2);
%   dI = 2000;
%   index2 = index1+dI;
  
  index1 = 1; index2 = Ls;

  h = signal(index1:index2);
  nT = length(h);
  if mod(nT,2) == 0; nT = nT-1; end
  h = h(1:nT);
  [n,m] = size(h);
  if n > 1; h = h'; end 
  h = h - mean(h);
  tMax = (nT-1)*dt;
  tMin = 0;
  t = linspace(tMin,tMax,nT);
      
 % Frequency domain     
   nF = 2999; fMax = 4500; fMin = -fMax;
   f = linspace(fMin,fMax,nF);

 % FOURIER TRANSFORM CALCULATIONS 
   H = zeros(1,nF); hI = zeros(1,nT);
   
% Fourier Transform  H(f)
   for c = 1:nF
     g = h.* exp(1i*2*pi*f(c)*t);
     H(c) = simpson1d(g,tMin,tMax);
   end


% INVERSE Fourier Transform  hI(t)  
 for c = 1:nT
   g = H.* exp(-1i*2*pi*t(c)*f);
   hI(c) = simpson1d(g,fMin,fMax);
 end  
  
% One-sided power spectral density PSD  Ph(f)
   psd = 2.*conj(H).*H;
   
% Total power  PT (time domain) and PF (frequenct domain)   
   PT = simpson1d(h.^2,tMin,tMax);
   PF = simpson1d(psd,fMin,fMax)./2;

% Output   
   fprintf('PT = %4.4f  \n \n',PT);
   fprintf('PF = %4.4f  \n \n',PF);
   
% GRAPHICS
figure(1)
   pos = [0.02 0.05 0.30 0.62];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   fs = 14;

subplot(3,1,1)   
   xP = ts; yP = signal;
   plot(xP,yP,'b','lineWidth',2);
%   title(txt_s,'FontWeight','normal')
   xlabel('t')
   ylabel('signal')
   grid on
   axis tight
  % xlim(XLIMS)
  % xticks(XTICKS)
   set(gca,'fontsize',fs)

subplot(3,1,2)   
   xP = t; yP = h;
   plot(xP,yP,'b','lineWidth',2);
   title(txt_h,'FontWeight','normal')
   xlabel('t')
   ylabel('h(t)')
   grid on
   axis tight
  % xlim(XLIMS)
  % xticks(XTICKS)
   set(gca,'fontsize',fs)
 
subplot(3,1,3)  
   xP = t; yP = real(hI);
   plot(xP,yP,'r','lineWidth',2);
%   title(txt_I,'FontWeight','normal')
   xlabel('t')
   ylabel('h_I(t)')
   grid on
   axis tight
  % xlim(XLIMS)
  % xticks(XTICKS)
   set(gca,'fontsize',fs)


figure(2)  % 222222222222222222222222222222222222222222222222
   pos = [0.35 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   subplot(2,1,1)   
   xP = f; yP = abs(H)./max(abs(H));
   plot(xP,yP,'b','lineWidth',2);
%   title(txt2,'FontWeight','normal')
   xlabel('f')
   ylabel('|H|')
   grid on
   xlim([0, fMax])
   set(gca,'fontsize',fs)
  % xlim(XLIMSF); xticks(XTICKSF)

subplot(2,1,2)  
   xP = f; yP = psd./max(psd);
   plot(xP,yP,'r','lineWidth',2);
   title('psd','FontWeight','normal')
   xlabel('f')
   ylabel('psd')
   grid on
   set(gca,'fontsize',fs)
 %  xlim(XLIMSF); xticks(XTICKSF)
   xlim([0, fMax])
 %  xlim([990 1010])
 %  xticks([990 1000 1010])

figure(3)
  pos = [0.35 0.55 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   spectrogram(y, ones(128,1), 0, 128, Fs, 'yaxis')


%% 
clear
close all
clc

r = [zeros(1,50), ones(1, 100), zeros(1, 50)];
rFT = abs(fft(r));
rFT = rFT./max(rFT); 
plot(rFT(1:100)) 


%%  
   clear
   close all
   clc

   global f
% Fourier transform

  h = [zeros(1,50), ones(1, 101), zeros(1, 50)];
  nT = length(h);
  w = hamming(nT);
  h = w';

  H = fourier(h);
 
  plot(f,abs(H))
  xlim([0 0.2])


function H = fourier(h)
  global f
% FOURIER TRANSFORM CALCULATIONS 
   nT = length(h);
   

   tMin = 0; tMax = 201;
   t = linspace(tMin,tMax,nT);
% Frequency domain     
   nF = 2990; fMax = 0.2; fMin = -fMax;
   f = linspace(fMin,fMax,nF);
   
   H = zeros(1,nF); hI = zeros(1,nT);
   
% Fourier Transform  H(f)
   for c = 1:nF
     g = h.* exp(1i*2*pi*f(c)*t);
     H(c) = simpson1d(g,tMin,tMax);
   end
end