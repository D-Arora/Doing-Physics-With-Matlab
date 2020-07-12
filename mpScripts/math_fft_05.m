% math_fft_01.m

% Fast Fourier Transform: decomposing a time series into
%  frequency components

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180712 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

%clear 
close all
%clc

% INPUTS / FUNCTION y(t) ==============================================

f = 10; T = 1/f;
nT = 150;
tMax = 5*T;
xT = linspace(0,tMax,nT);
yT = sin(2*pi*f*xT) ;%+ sin(2*pi*100*xT) + 0.5*sin(2*pi*200*xT);
%yT(yT > 0) = 1;
%yT(yT<0) = -1;

TS = xT(2) - xT(1);
fS = 1/TS;


nF = 1024;
yF = fft(yT,nF);
%nVals = 0:nF-1;
xF = 0:length(yF)-1;

  
  


% Generate the plot, title and labels.
figure(1);
   pos = [0.05 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = xT; yP = yT;
   plot(xP,yP,'lineWidth',2);
   title('Function y(t)')
   xlabel('Time t (s)')
   ylabel('y(t)')
   grid on

figure(2);
   pos = [0.31 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = xF; yP = abs(yF);
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (no shift)');
   xlabel('FFT Sample Domain');
   ylabel('FFT y(t)');
   grid on
 
figure(3);
   pos = [0.57 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = xF./nF; yP = abs(yF);
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (no shift)');
   xlabel('Normalized frequency');
   ylabel('FFT y(t)');
   grid on   
   
figure(4);
   pos = [0.05 0.38 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = fS.*(-nF/2:nF/2-1)./nF;
   yP = abs(fftshift(yF));
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (FFT shift)');
   xlabel('frequency f_F  [Hz]');
   ylabel('FFT y(t)');
   grid on 
   
% figure(5);
%    pos = [0.31 0.38 0.25 0.25];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    xP = fS.*(-nF/2:nF/2-1)./nF;
%    yP = real(fftshift(yF));
%    xP = xP(1+nF/2:nF);
%    yP = yP(1+nF/2:nF);
%    plot(xP,yP,'lineWidth',2);
%    title('Double Sided FFT (FFT shift)');
%    xlabel('frequency f_F  [Hz]');
%    ylabel('FFT y(t)');
%    grid on    

figure(5);
   pos = [0.31 0.38 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = fS.*(-nF/2:nF/2-1)./nF;
   yP = fftshift(yF);
   yP = conj(yP).*yP ./ (nF*length(yT));
   xP = xP(1+nF/2:nF);
   yP = yP(1+nF/2:nF);
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (FFT shift)');
   xlabel('frequency f_F  [Hz]');
   ylabel('Power Spectrum');
   grid on   
   
figure(6)
   pos = [0.57 0.38 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
%    xP = fS.*(-nF/2:nF/2-1)./nF;
%    yP = real(fftshift(yF));
%    yP = conj(yP).*yP ./ (nF*length(yT));
%    xP = xP(1+nF/2:nF);
%    yP = yP(1+nF/2:nF);
   yP = 10.*log10(yP);
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (FFT shift)');
   xlabel('frequency f_F  [Hz]');
   ylabel('Normalized Power Spectrum');
   grid on     
   