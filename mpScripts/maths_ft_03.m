% maths_ft_01.m

% Fourier Transform by direct integration: decomposing a time series into
%  its frequency components
% MANIPULATING SPECTRUM DATA


% Ian Cooper
% email: matlabvisualphysics@gmail.com
% School of Physics, University of Sydney
% 230402 / Matlab version R2021ba

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/maths_ft01.pdf

% Number of time grid points nT (must be odd)
% Number of frequency grid points  nF (must be odd)
% Signal function h(t)    (must be a column vector)
% Fourier transform H(t)  (must be a column vector)
% Script uses the function simpson1d.m

clear; close all; clc

tic
% ===================================================================
% Filtering: flagF = 1 --> ON   /    flagF = 0 --> OFF
  flagF = 1;

% Filter cut-off frequency (high pass filter)
  fC = 1000;

% Audio signals (select: comment / uncomment
  flagS = 6;

  switch flagS

  case 1   
        [signal, Fs] = audioread('audioClarinet1.wav');
        txt_s = 'Clarinet 220 Hz';
  case 2
       [signal, Fs] = audioread('audioGuitar1.wav');
       txt_s = 'Guitar 220 Hz';
  case 3
       [signal, Fs] = audioread('audioVoice1.wav'); 
       txt_s = 'Voice 220 Hz';
 case 4
       [signal, Fs] = audioread('audio440.wav');
       txt_s = 'audio signal 440 Hz Hz'; 
 case 5   
      [signal, Fs] = audioread('Train.wav');
      txt_s = 'train whistle'; 
case 6
      [signal, Fs] = audioread('audio_secret_message.wav');
      txt_s = 'Secret message'; 
   
end  

% Play original signal
   sound(signal, Fs);
% Time domain   
  Ls = length(signal);
  dt = 1/Fs;
  ts = linspace(0,Ls*dt,Ls);
       
% Function h(t) 
  txt_h = 'h(t)';

  index1 = round(Ls/2);
  dI = 2000;
  index2 = index1+dI;
  
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
   nF = 2999; fMax = 1500; fMin = -fMax;
   f = linspace(fMin,fMax,nF);

  
% ==================================================================
% FOURIER TRANSFORM CALCULATIONS 
   H = zeros(1,nF); hI = zeros(1,nT);
   
% Fourier Transform  H(f)
   for c = 1:nF
     g = h.* exp(1i*2*pi*f(c)*t);
     H(c) = simpson1d(g,tMin,tMax);
   end

% Filtering 
    txt_I = 'h_{inv}(t)   Filter OFF';
    txt2 = 'inv FT Filter OFF';
    if flagF == 1
        H(f > 500) = 0;
        txt_I = 'h_{inv}(t)   Filter ON';
        txt2 = 'inv FT Filter OFF';
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
  

% PEAKS ================================================================
  [A, fS, w, p] = findpeaks(abs(H)./max(abs(H)), f, 'MinPeakProminence',0.05, 'MinPeakDistance', 80);
  audioData = table(fS', A', w', p', 'VariableNames', {'Frequency', 'PeakHeight', 'PeakWidth', 'Prominence'})


% CREATE and PLAY SOUND
%%
% Frequency inputs
%  f1 = 660;   A1 = 1;
%  f2 = 1099;  A2 = 0.44;
%  f3 = 1320;  A3 = 0.12;
%  
% Calculate Waveform
 fs = 22050; % sample frequency (Hz)
 d = 4.0; % duration (s)
 n = fs * d; % number of samples
 tS = (1:n) / fs; % sound data preparation
 % s = sin(2 * pi * f1 * s); % pure tone
 s = zeros(1,n);
 for c = 1:length(A)
     if fS(c) > 0
       s = s + A(c).*sin(2*pi*fS(c)*tS);
     end
 end
 s = s./max(s);
%sound(s, fs); % Generate sound
%pause(d + 0.5); % waiting for sound end
% Save wav file to disk
% filename = 'wav_S3000-3003.wav';
% audiowrite(filename,s,fs);

%%


  
% GRAPHICS ==============================================================
  fs = 14;

figure(1)
   pos = [0.02 0.05 0.30 0.62];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

subplot(3,1,1)   
   xP = ts; yP = signal;
   plot(xP,yP,'b','lineWidth',2);
   title(txt_s,'FontWeight','normal')
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
   title(txt_I,'FontWeight','normal')
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
   title(txt2,'FontWeight','normal')
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

figure(3)   % computed signal
   pos = [0.35 0.5 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = tS; yP = s;
   plot(xP,yP,'b','lineWidth',2);
   title('Computed Signal')
   xlabel('t_S')
   ylabel('s')
   grid on
   set(gca,'fontsize',fs)
   xlim([1 1.04])

toc

% Play manipualted signal
%   sound(abs(hI), Fs);





% FUNCTIONS ===========================================================

function integral = simpson1d(f,a,b)

% [1D] integration - Simpson's 1/3 rule
%      f function    a = lower bound    b = upper bound
%      Must have odd number of data points
%      Simpson's coefficients   1 4 2 4 ... 2 4 1

numS = length(f);               % number of data points

if mod(numS,2) == 1
sc = 2*ones(numS,1);
sc(2:2:numS-1) = 4;
sc(1) = 1; sc(numS) = 1;

h = (b-a)/(numS-1);

integral = (h/3) * f * sc;

else
    
integral = 'Length of function must be an ODD number' 
end

end




