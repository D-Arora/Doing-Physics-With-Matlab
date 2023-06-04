% maths_ft_01.m

% Fourier Transform by direct integration: decomposing a time series into
%  its frequency components

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

% SELECT SIGNAL FUNCTION h(t)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  flagF = 15;

%  1   Gausian function
%  2   Exponential function
%  3   Sinusoidal function
%  4   Superposition of sinusoidal functions
%  5   Square Wave
%  6   Sawtooth function
%  7   Single square pulse
%  8   Damped sinusoidal function
%  9   ECG
% 10   Beats
% 11   Beats: audio file
% 12   Audio file: 440 Hz signal
% 13   Audio file: Guitar 220 Hz
% 14   Audio file: Clarient 220 Hz
% 15   Audio file: Voice 220 Hz
% 16   Audio file: train whistle
% 17   Digital Filtering  flagFF = 0 (off) / flagFF = 1 (on)
       flagFF = 0;
% ===================================================================

switch flagF

case 1  % Gaussian function
     txt = 'Gaussian signal function h(t)';
     nT = 2001; nF = 2001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT-1; end
     if mod(nF,2) == 0; nF = nF-+1; end
     tMin = -3; tMax = 3;
     fMax = 2; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:1:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:0.5:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h = exp(-10.*t.^2);  % must be a column vector
    
    % [n,m] = size(h);
    % if n > 1; h = h'; end

case 2  % Exponential function
     txt = 'Exponential signal function h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end
     tMin = 0; tMax = 2;
     fMax = 10; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.5:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:5:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h = 2.*exp(-3*t);  % must be a column vector
    
     [n,m] = size(h);
     if n > 1; h = h'; end

  %   HT = 2./(3 + 1i.*2.*pi.*f);
 

case 3  % Sinusoidal function
     txt = 'Sinuoidal signal function h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     f0 = 10; A = 1; phi = pi/6;
     T0 = 1/f0;
     tMin = 0; tMax = 4*T0; fMax = 40; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.1:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:5:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h = A.*sin(2*pi*f0*t + phi);  % must be a column vector
    
     [n,m] = size(h);
     if n > 1; h = h'; end      
         
 
 case 4  % superposition of sinusoidal functions
     txt = 'Sinuoidals (superposition) h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     f0 = 20;
     tMin = 0; tMax = 0.2; fMax = 80; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.1:tMax;
     XLIMSF = [fMin fMax]; XTICKSF = fMin:20:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h = 1.*sin(2*pi*1*f0*t) + 0.5.*sin(2*pi*2*f0*t) + 0.75.*sin(2*pi*3*f0*t);   % must be a column vector
    
     [n,m] = size(h);
     if n > 1; h = h'; end 
              

case 5  % Square wave
     txt = 'Square wave singal h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     f0 = 10; A = 1; phi = 0;
     T0 = 1/f0;
     tMin = 0; tMax = 3*T0; fMax = 120; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.05:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:30:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

      h = A.*sin(2*pi*f0*t + phi);
         h(h > 0) = 1;
         h(h<0) = -1; 
    
     [n,m] = size(h);
     if n > 1; h = h'; end 
                  
    
case 6   % sawtooth wave
     txt = 'Sawtooth signal function h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     T0 = 5; 
     tMin = 0; tMax = 10*T0; fMax = 1.00; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:10:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:0.2:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

      h = 1.*rem(t,T0);
          h(h == 0 | f == 10) = 5;
          h = h - mean(h);
    
     [n,m] = size(h);
     if n > 1; h = h'; end 
                  
       
 case 7  % Single square pulse
     txt = 'Single square pulse h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     tMin = -2.0; tMax = 2.0; fMax = 5; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.5:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:1:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h(t > -1) = 0.5; 
     h(t > 1)  = 0;
     

     [n,m] = size(h);
     if n > 1; h = h'; end 

                  
case 8  % Damped oscillator 
     txt = 'Damped oscillator (t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     f0 = 10; A = 1; 
     T0 = 1/f0;  tau = 0.2;
     tMin = 0; tMax = 4*T0; fMax = 40; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.1:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:5:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h = A.*exp(-t/tau).*sin(2*pi*f0*t);  % must be a column vector
    
     [n,m] = size(h);
     if n > 1; h = h'; end      
         
            
case 9   % ecg
     txt = 'ECG signal h(t)';
     load('ecg.mat','ecg')
     h = ecg(1:end)';
     nT = length(h); nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT-1; end
     if mod(nF,2) == 0; nF = nF-1; end

     tMin = 0; tMax = 3600/200; fMax = 8; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:2:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:1:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

      h = h(1:nT);
      h = h - mean(h);
    
     [n,m] = size(h);
     if n > 1; h = h'; end      
         

case 10 % Beats
     txt = 'BEATS  h(t)';
     nT = 5001; nF = 5001;    % must be odd numbers
     if mod(nT,2) == 0; nT = nT+1; end
     if mod(nF,2) == 0; nF = nF+1; end

     A1 = 1; A2 = 1;
     f1 = 1000; f2 = 1008;
     phi1 = 0; phi2 = 0;
     tMin = 0; tMax = 2; fMax = 1200; fMin = -fMax;
     XLIMS = [tMin  tMax]; XTICKS = tMin:0.5:tMax; 
     XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax; 
     [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     h = A1.*sin(2*pi*f1*t + phi1) + A2.* sin(2*pi*f2*t+ phi2);  % must be a column vector
    
     [n,m] = size(h);
     if n > 1; h = h'; end  


 case 11  % BEATS - sound file
       txt = 'BEATS  h(t)'; 
       [signal,Fs] = audioread('wav_S1000_1008.wav');
       sound(signal,Fs);
       dt = 1/Fs;
       h = signal;
       nT = length(signal);
       if mod(nT,2) == 0; nT = nT-1; end
       h = signal(1:nT);
       tMax = dt*(nT-1); tMin = 0;

       nF = 7999; if mod(nF,2) == 0; nF = nF+1; end
       fMax = 1200; fMin = -fMax;
       XLIMS = [tMin  tMax]; XTICKS = tMin:0.5:tMax; 
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax; 
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
       [n,m] = size(h);
       if n > 1; h = h'; end  


case 12  % audio 440 Hz signal
       txt = 'audio signal 440 Hz   h(t)'; 
       [signal, Fs] = audioread('audio440.wav');
       sound(signal, Fs);
       dt = 1/Fs;
       h = signal;
       nT = length(h);
       if mod(nT,2) == 0; nT = nT-1; end
       h = h(1:nT);
       [n,m] = size(h);
       if n > 1; h = h'; end 
       h = h - mean(h);
       nF = 7999; 
       tMax = (nT-1)*dt;
       tMin = 0;
       fMax = 1200; fMin = -fMax;
       XLIMS = [tMin  0.01]; XTICKS = tMin:0.01/5:0.01; 
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax;
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

case 13  % audio Guitar 220 Hz
       txt = 'Guitar 220 Hz signal  h(t)'; 
       [signal, Fs] = audioread('audioGuitar1.wav');
       sound(signal, Fs);
       dt = 1/Fs;
       h = signal;
       nT = length(h);
       if mod(nT,2) == 0; nT = nT-1; end
       h = h(1:nT);
       [n,m] = size(h);
       if n > 1; h = h'; end 
       h = h - mean(h);
       nF = 7999; 
       tMax = (nT-1)*dt;
       tMin = 0;
       fMax = 800; fMin = -fMax;
       XLIMS = [tMin  3]; XTICKS = tMin:0.5:3; 
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax;
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
       
  
 case 14  % audio Clarinet 220 Hz
       txt = 'Clarinet 220 Hz signal  h(t)'; 
       [signal, Fs] = audioread('audioClarinet1.wav');
       sound(signal, Fs);
       dt = 1/Fs;
       h = signal;
       nT = length(h);
       if mod(nT,2) == 0; nT = nT-1; end
       h = h(1:nT);
       [n,m] = size(h);
       if n > 1; h = h'; end 
       h = h - mean(h);
       nF = 7999; 
       tMax = (nT-1)*dt;
       tMin = 0;
       fMax = 800; fMin = -fMax;
       XLIMS = [tMin  3]; XTICKS = tMin:0.5:3; 
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax;
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);


 case 15  % audio Voice 220 Hz
       txt = 'Voice 220 Hz signal  h(t)'; 
       [signal, Fs] = audioread('audioVoice1.wav');
       sound(signal, Fs);
       dt = 1/Fs;
       h = signal(1e5:2.5e5);
       nT = length(h);
       if mod(nT,2) == 0; nT = nT-1; end
       h = h(1:nT);
       [n,m] = size(h);
       if n > 1; h = h'; end 
       h = h - mean(h);
       nF = 7999; 
       tMax = (nT-1)*dt;
       tMin = 0;
       fMax = 800; fMin = -fMax;
       XLIMS = [tMin  2.5]; XTICKS = tMin:0.5:2.5; 
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax;
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);                 
       
       
       % turn filtering on  flagFF = 17 otherwise flagFF = 0
       flagFF = 15;


case 16  % audio Train Whistle
       txt = 'Train Whistle   h(t)'; 
       [signal, Fs] = audioread('Train.wav');
       sound(signal, Fs);
       dt = 1/Fs;
       h = signal;
       nT = length(h);
       if mod(nT,2) == 0; nT = nT-1; end
       h = h(1:nT);
       [n,m] = size(h);
       if n > 1; h = h'; end 
       h = h - mean(h);
       nF = 7999; 
       tMax = (nT-1)*dt;
       tMin = 0;
       fMax = 1200; fMin = -fMax;
       XLIMS = [tMin  1.5]; XTICKS = tMin:0.5:1.5; 
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax;
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);  
       
case 17  % Filtering
     % Filter off 0 / filter on 15
       flagFF = 15;
       txt = 'Filtering OFF';
       nT = 5001; nF = 5001;    % must be odd numbers
       if mod(nT,2) == 0; nT = nT+1; end
       if mod(nF,2) == 0; nF = nF+1; end

       f0 = 1000;
       tMin = 0; tMax = 0.02; fMax = 1200; fMin = -fMax;
       XLIMS = [tMin  tMax]; XTICKS = tMin:0.005:tMax;
       XLIMSF = [fMin fMax]; XTICKSF = fMin:200:fMax; 
       [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

       h = 1.*sin(2*pi*1*f0*t) + 0.8.*sin(2*pi*2*100*t);   % must be a column vector
    
     [n,m] = size(h);
     if n > 1; h = h'; end 

end


% FOURIER TRANSFORM CALCULATIONS ======================================
   H = zeros(1,nF); hI = zeros(1,nT);
%   h(1) = 0; h(end) = 0;
   
% Fourier Transform  H(f)
   for c = 1:nF
     g = h.* exp(1i*2*pi*f(c)*t);
     H(c) = simpson1d(g,tMin,tMax);
   end
 
% #17 Filtering
  if flagFF == 17
      H(f<500) = 0;
      txt = 'Filtering ON';
      XLIMSF = [0 fMax]; XTICKSF = 0:200:fMax; 
  end

% #15 Filtering  
   if flagFF == 15
      H(f<300) = 0;
      txt = 'Filtering ON';
      XLIMSF = [0 fMax]; XTICKSF = 0:200:fMax; 
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

   fprintf('PT = %4.4f  \n \n',PT);
   fprintf('PF = %4.4f  \n \n',PF);
  
   

% GRAPHICS ==============================================================
  fs = 14;

figure(1)
   pos = [0.02 0.05 0.30 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
subplot(2,1,1)   
   xP = t; yP = h;
   plot(xP,yP,'b','lineWidth',2);
   title(txt,'FontWeight','normal')
   xlabel('t')
   ylabel('h(t)')
   grid on
   xlim(XLIMS)
   xticks(XTICKS)
   set(gca,'fontsize',fs)
 subplot(2,1,2)  
   xP = t; yP = real(hI);
   plot(xP,yP,'r','lineWidth',2);
   title('Inverse F.T. h_I(t)','FontWeight','normal')
   xlabel('t')
   ylabel('h_I(t)')
   grid on
   xlim(XLIMS)
   xticks(XTICKS)
   set(gca,'fontsize',fs)


figure(2)
   pos = [0.35 0.05 0.30 0.33];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
subplot(2,1,1)   
   xP = f; yP = abs(H);
   plot(xP,yP,'b','lineWidth',2);
   title('F.T. |H(f)|','FontWeight','normal')
   xlabel('f')
   ylabel('|H|')
   grid on
   set(gca,'fontsize',fs)
   xlim(XLIMSF); xticks(XTICKSF)

subplot(2,1,2)  
   xP = f; yP = psd./max(psd);
   plot(xP,yP,'r','lineWidth',2);
   title('psd','FontWeight','normal')
   xlabel('f')
   ylabel('psd')
   grid on
   set(gca,'fontsize',fs)
   xlim(XLIMSF); xticks(XTICKSF)
   xlim([0, fMax])
 %  xlim([990 1010])
 %  xticks([990 1000 1010])
toc

% FUNCTIONS ===========================================================

function [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF)
   t = linspace(tMin,tMax,nT);
   f = linspace(fMin,fMax,nF);
end

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
    
integral = 'Length of function must be an ODD number'; 
end

end
