% math_fft_01.m

% Fast Fourier Transform: decomposing a time series into
%  frequency components

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180712 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

clear 
close all
clc

% INPUTS / FUNCTION y(t) ==============================================

% Number of grid points (power of 2)
   N = 126;
% Time steps
   f = 200;
   T = 1/f;
   tMax = 2*T;
   t = linspace(0,tMax,N);
% Function 
   y = 5.*sin(2*pi*f*t);
% Sampling period
   TS = t(2)-t(1);
% Sampling frequency
  fS = 1/TS;
  FTf = fS .* (0:N/2)/N;
% Compute fft
  F = fft(y,1024);
  
  Fp = F(1:N/2+1).*TS;
  
  


% Generate the plot, title and labels.
figure(1);
   pos = [0.1 0.1 0.3 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   plot(t,y);
   title('Sine Wave Signal');
   xlabel('Time (s)');
   grid on

figure(2);
   pos = [0.4 0.1 0.3 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   plot(FTf,real(Fp));
   title('Power Spectrum of a Sine Wave');
   xlabel('Frequency (Hz)');
  % xlim([0 500]);
   ylabel('Power');
   grid on
   