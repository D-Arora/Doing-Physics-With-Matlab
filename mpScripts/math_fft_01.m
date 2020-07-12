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

% Sampling frequency
  fS = 150;
  
t = 0:1/fS:1; % Time vector of 1 second
f = 10; % Create a sine wave of f Hz.
x = sin(2*pi*t*f);


% 
nfft = 1024; % Length of FFT
% Take fft, padding with zeros so that length(X)
X = fft(x,nfft);
% FFT is symmetric, throw away second half
%Power Spectrum of a Sine Wave
% FFT is symmetric, throw away second half
X = X(1:nfft/2);
% Take the magnitude of fft of x
mx = abs(X);
% Frequency vector
f = (0:nfft/2-1)*fS/nfft;


% Generate the plot, title and labels.
figure(1);
plot(t,x);
title('Sine Wave Signal');
xlabel('Time (s)');

figure(2);
plot(f,mx);
title('Power Spectrum of a Sine Wave');
xlabel('Frequency (Hz)');
ylabel('Power');