% math_fft_06.m

% Fourier Transform by direct integration: decomposing a time series into
%  frequency components

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180712 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

clear 
close all
clc

tic
% INPUTS / FUNCTION h(t) ==============================================
% Time: sample points (must be an odd integer)  [2001]
% Frequency: sample points (must be an odd integer)   [2001]
   nT = 2001;
   nF = 2001;

% Function h(t)   flagF
% 1   Gausian function
   flagF = 1;
   
   switch flagF
       case 1
         tMin = -3;
         tMax = 3;
         fMin = -2;
         fMax = 2;
         
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = exp(-pi.*t.^2);
       
       case 2
           
   end
% f0 = 20;
% T0 = 1/f0;
% 
% nT = 1801;
% nF = 5801;



 %h = sin(2*pi*f0*t) + 0.5.*sin(2*pi*2*f0*t) + 0.25.* sin(2*pi*4*f0*t) + 0.8.*sin(2*pi*6*f0*t);
%  h = sin(2*pi*f0*t);
%  h(h > 0) = 1;
%  h(h<0) = -1;

%h = 2.*exp(-3*t);



H = zeros(1,nF); hI = zeros(1,nT);
HT = zeros(1,nF);

% Fourier Transform
for c = 1:nF
   g = h.* exp(1i*2*pi*f(c)*t);
   H(c) = simpson1d(g,tMin,tMax);
%   HT(c) = 2./(3 + 1i.*2*pi*f(c));
   
end
  H = H ./ tMax;

% INVERSE Fourier Transform    
for c = 1:nT
   g = H.* exp(-1i*2*pi*t(c)*f);
   hI(c) = simpson1d(g,fMin,fMax);
end  
  
  hI = 2.*hI ./ tMax;
  
% Total power
   P = conj(H).*H;
   p_total = simpson1d(h.^2,tMin,tMax)
   p_avg = p_total / tMax
   P_total = 2.*simpson1d(conj(H).*H,fMin,fMax)
  
   




% Generate the plot, title and labels.
figure(1);
   pos = [0.05 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = t; yP = h./max(h);
   plot(xP,yP,'lineWidth',3);
   title('Function y(t)')
   xlabel('Time t (s)')
   ylabel('y(t)')
   grid on
   
   hold on
   xP = t; yP = real(hI)./max(real(hI));
   plot(xP,yP,'r','lineWidth',1,'lineStyle','--');

figure(2);
   pos = [0.31 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = f; yP = abs(H)./max(abs(H));
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (no shift)');
   xlabel('FFT Sample Domain');
   ylabel('FFT y(t)');
   grid on
   hold on
   yP = abs(HT);
   plot(xP,yP,'r')
   
   
figure(3);
   pos = [0.57 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = f; yP = P./max(P);
   plot(xP,yP,'lineWidth',2);
   title('Double Sided FFT (no shift)');
   xlabel('Normalized frequency');
   ylabel('power P');
   grid on   
   
% figure(4);
%    pos = [0.05 0.38 0.25 0.25];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    xP = fS.*(-nF/2:nF/2-1)./nF;
%    yP = abs(fftshift(yF));
%    plot(xP,yP,'lineWidth',2);
%    title('Double Sided FFT (FFT shift)');
%    xlabel('frequency f_F  [Hz]');
%    ylabel('FFT y(t)');
%    grid on 
%    
% % figure(5);
% %    pos = [0.31 0.38 0.25 0.25];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w');
% %    xP = fS.*(-nF/2:nF/2-1)./nF;
% %    yP = real(fftshift(yF));
% %    xP = xP(1+nF/2:nF);
% %    yP = yP(1+nF/2:nF);
% %    plot(xP,yP,'lineWidth',2);
% %    title('Double Sided FFT (FFT shift)');
% %    xlabel('frequency f_F  [Hz]');
% %    ylabel('FFT y(t)');
% %    grid on    
% 
% figure(5);
%    pos = [0.31 0.38 0.25 0.25];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    xP = fS.*(-nF/2:nF/2-1)./nF;
%    yP = fftshift(yF);
%    yP = conj(yP).*yP ./ (nF*length(yT));
%    xP = xP(1+nF/2:nF);
%    yP = yP(1+nF/2:nF);
%    plot(xP,yP,'lineWidth',2);
%    title('Double Sided FFT (FFT shift)');
%    xlabel('frequency f_F  [Hz]');
%    ylabel('Power Spectrum');
%    grid on   
%    
% figure(6)
%    pos = [0.57 0.38 0.25 0.25];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
% %    xP = fS.*(-nF/2:nF/2-1)./nF;
% %    yP = real(fftshift(yF));
% %    yP = conj(yP).*yP ./ (nF*length(yT));
% %    xP = xP(1+nF/2:nF);
% %    yP = yP(1+nF/2:nF);
%    yP = 10.*log10(yP);
%    plot(xP,yP,'lineWidth',2);
%    title('Double Sided FFT (FFT shift)');
%    xlabel('frequency f_F  [Hz]');
%    ylabel('Normalized Power Spectrum');
%    grid on     

toc

function [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF)
   t = linspace(tMin,tMax,nT);
   f = linspace(fMin,fMax,nF);
end

