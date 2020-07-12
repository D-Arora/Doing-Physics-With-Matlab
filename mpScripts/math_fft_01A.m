% math_fft_01.m

% Fast Fourier Transform: decomposing a time series into
%  frequency components

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180712 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm

clear 
close all
clc

% INPUTS / FUNCTION y(t) ==============================================

N = 128;
t = linspace(0,3,N);
h = 2.*exp(-3*t);

Ts = t(2)-t(1);
fs = 1/Ts;
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
hold on
plot(f,2./(3+1i.*2.*pi.*f),'r');