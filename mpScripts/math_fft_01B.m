% math_fft_01.m

% Fast Fourier Transform: decomposing a time series into
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

% INPUTS / FUNCTION y(t) ==============================================

N = 256;
t = linspace(0,0.5,N);
h = sin(2*pi*10*t);

Ts = t(2)-t(1);
fs = 10;
Ws = 2*pi/Ts;
H = fft(h);
Hp = H(1:N/2+1).*Ts;
f = fs .* (0:N/2)/N;


% Generate the plot, title and labels.
figure(1);
plot(t,h);
title('Sine Wave Signal');
xlabel('Time (s)');

figure(2);
plot(f,Hp);
title('Power Spectrum of a Sine Wave');
xlabel('Frequency (Hz)');
ylabel('Power');
