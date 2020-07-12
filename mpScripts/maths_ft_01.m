% maths_ft_01.m

% Fourier Transform by direct integration: decomposing a time series into
%  its frequency components

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180712 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

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
% 2   Exponential function
% 3   Sinusoidal function
% 4   Superposition of sinusoidal functions
% 5   Square Wave
% 6   Sawtooth function
% 7   Single square pulse
% 8   Damped sinusoidal function
% 9   ecg
%10   Beats
% 11  Train Whistle
% 12  Digital Filtering

flagF = 4;
   
   switch flagF
       case 1         % Gaussian function
         tMin = -3;
         tMax = 3;
         fMax = 2;
         fMin = -fMax;
         XTICKS = -2:0.5:2;
         XLIMS = [0 fMax];
         
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = exp(-pi.*t.^2);
       
       case 2        % Exponential function
         tMin = 0;
         tMax = 3;
         fMax = 30;
         fMin = - fMax;
         XTICKS = 0:5:fMax;
         XLIMS = [fMin fMax];
         
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = 2.*exp(-3*t);
         HT = 2./(3 + 1i.*2.*pi.*f);
         
       case 3          % sinusoidal function
          f0 = 10;
          A = 1;
          phi = pi/6;
          T0 = 1/f0;
          tMin = 0;
          tMax = 4*T0;
          fMax = 50;
          fMin = - fMax;
          XTICKS = -fMax:10:fMax;
          XLIMS = [-fMax fMax];
          
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = A.*sin(2*pi*f0*t + phi);
         
       case 4          % superposition of sinusoidal functions
          nF = 2049;
          f0 = 1;
          A = 1:4;
          T0 = 1/f0;
          tMin = 0;
          tMax = 4*T0;
          fMax = 5;
          fMin = -fMax;
          XTICKS = -fMax:1:fMax;
          XLIMS = [-fMax fMax];
          
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         %h = A(1).*sin(2*pi*f0*t) + A(2).*sin(2*pi*2*f0*t) + A(3).*sin(2*pi*3*f0*t) + A(4).*sin(2*pi*4*f0*t);
          h = 1.*sin(2*pi*1*t) + 0.5.*sin(2*pi*2*t) + 0.25.*sin(2*pi*3*t);   
  
      
         case 5       % square wave
          f0 = 10;
          A = 1;
          phi = 0;
          T0 = 1/f0;
          tMin = 0;
          tMax = 3*T0;
          fMax = 120;
          fMin = - fMax;
          XTICKS = -fMax:30:fMax;
          XLIMS = [-fMax fMax];
          
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = A.*sin(2*pi*f0*t + phi);
         h(h > 0) = 1;
         h(h<0) = -1; 
         
         
       case 6       % sawtooth wave
          T0 = 10;          
          tMin = 0;
          tMax = 5*T0;
          fMax = 1.2;
          fMin = - fMax;
          XTICKS = -fMax:0.2:fMax;
          XLIMS = [-fMax fMax];
           
          [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);   
          
          h = 1.*rem(t,T0);
          h(h == 0 | f == 10) = 5;
          h = h - mean(h);
          
       
       case 7         % single square pulse
          T0 = 1.0;          
          tMin = 0;
          tMax = 5*T0;
          fMax = 10;
          fMin = - fMax;
          XTICKS = -fMax:2:fMax;
          XLIMS = [-fMax fMax];
           
          [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);   
          
          h = zeros(1,nT);
          h(round(nT/3) : round(2*nT/3)) = 1;
          h = h - 0.5;
          
       
       case 8      % Damped oscillator  
         f0 = 10;
          A = 1;
          tau = 0.1;
          T0 = 1/f0;
          tMin = 0;
          tMax = 4*T0;
          fMax = 50;
          fMin = - fMax;
          XTICKS = -fMax:10:fMax;
          XLIMS = [-fMax fMax];
          
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = A.*exp(-t/tau).*sin(2*pi*f0*t);  
          
       case 9       % ecg
          tMin = 0;
          tMax = 3600/200;
          fMax = 8; %4*0.005;
          fMin = -fMax;
          XTICKS = -fMax:fMax/2:fMax;
          XLIMS = [-fMax fMax];
          load('ecg.mat','ecg')
          h = ecg(1:end-1)';
          h = h - mean(h);
          nT = length(h); nF = 7999;%length(h);
          [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
          


      case 10       % Beats
          tMin = 0;
          tMax = 4;
          fMax = 1100;
          fMin = -fMax;
          XTICKS = -fMax:fMax/2:fMax;
          XLIMS = [-fMax fMax];
          [a, b] = audioread('wav_S1000_1008.wav');
          sound(a,b);
          h = a(1:end-1)';
          h = h - mean(h);
          nT = length(h); nF = 7999;%length(h);
          [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);

     case 11       % Beats
          tMin = 0;
          tMax = 4;
          fMax = 1100;
          fMin = -fMax;
          XTICKS = -fMax:fMax/2:fMax;
          XLIMS = [-fMax fMax];
          [a, b] = audioread('Train.wav');
          sound(a,b);
          h = a(1:end-1)';
          h = h - mean(h);
          nT = length(h); nF = 7999;%length(h);
          [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
     
   
    case 12          % sinusoidal function
          f0 = 1000;
          A = 1;
          phi = 0;
          T0 = 1/f0;
          tMin = 0;
          tMax = 10*T0;
          fMax = 1800;
          fMin = - fMax;
          XTICKS = -fMax:600:fMax;
          XLIMS = [-fMax fMax];
          
         [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF);
         
         h = A.*sin(2*pi*f0*t + phi) + 1.* sin(2*pi*100*t);
   
   
   end


% FOURIER TRANSFORM CALCULATIONS ======================================
   H = zeros(1,nF); hI = zeros(1,nT);
   h(1) =0; h(end) = 0;
   
  
   
% Fourier Transform  H(f)
   for c = 1:nF
     g = h.* exp(1i*2*pi*f(c)*t);
     H(c) = simpson1d(g,tMin,tMax);
   end
 
 % Filtering  (uncomment for filtering effects)
 % H(f<500) = 0;

% INVERSE Fourier Transform  hI(t)  
for c = 1:nT
   g = H.* exp(-1i*2*pi*t(c)*f);
   
   hI(c) = simpson1d(g,fMin,fMax);
end  
 
  
% One-sided power spectral density PSD  Ph(f)
   Ph = 2.*conj(H).*H;

% Total power  PT (time domain) and PF (frequenct domain)   
   PT = simpson1d(h.^2,tMin,tMax);
   PF = simpson1d(Ph,fMin,fMax)./2;

   fprintf('PT = %4.4f  \n \n',PT);
   fprintf('PF = %4.4f  \n \n',PF);
  
   

% GRAPHICS ==============================================================
   fs = 14;

figure(1)   % h(t) and hI(t)
   pos = [0.02 0.05 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = t; yP = h;
   plot(xP,yP,'lineWidth',3);
   title('Function h(t)')
   xlabel('t')
   ylabel('h(t)')
   grid on
   
   hold on
   %xP = t(1:20:nT); yP = real(hI(1:20:nT));
   xP = t; yP = real(hI);
   plot(xP,yP,'r','lineWidth',1);
   set(gca,'fontsize',fs)
   legend('h','hI')

figure(2);
   pos = [0.31 0.05 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = f; yP = abs(H);
   plot(xP,yP,'lineWidth',2);
   title('Fourier Tranform | H(f) |');
   xlabel('f');
   ylabel('| H(f) |');
   grid on
   set(gca,'fontsize',fs)
   set(gca,'xtick', XTICKS)
   set(gca,'xlim', XLIMS)
   yP = angle(H)./pi;
   
%    yyaxis right
%    plot(xP,yP,'lineWidth',0.2)
%    ylabel('phase / \pi  [rad]');
%    
   if flagF == 2
      hold on
      xP = f(1:20:nF); yP = abs(HT(1:20:nF));
      plot(xP,yP,'+r','lineWidth',0.2);
   end
   
figure(3);
   pos = [0.60 0.05 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = f; yP = Ph;
   plot(xP,yP,'lineWidth',2);
   title('One-sided PSD');
   ylabel('PSD   P_h(f)');
   xlabel('f');
   grid on   
   set(gca,'fontsize',fs)
   xlim([0 fMax])


figure(4)   % real(H) and imag(H)
   pos = [0.02 0.40 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   xP = f; yP = real(H);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   yP = imag(H);
   plot(xP,yP,'r','lineWidth',2);
   legend('Re','Im');
   title('Real and Imaginary parts of H(f)')
   xlabel('f')
   ylabel('H)')
   grid on   
   set(gca,'fontsize',fs)
   set(gca,'xtick', XTICKS)
   set(gca,'xlim', XLIMS)
   
toc

% FUNCTIONS ===========================================================

function [t, f] = domains(tMin,tMax,nT,fMin,fMax,nF)
   t = linspace(tMin,tMax,nT);
   f = linspace(fMin,fMax,nF);
end

