% bp006.m

% External insulin injection
% Gaussian function

clear
close all
clc


N = 501;

S = zeros(N,1);

S0 = 1;

s = 0.8*60;

t = linspace(0,24*60,N);
ts = 5*60; tD = 3*60;
S = S0.*exp(-(t - (ts+tD)).^2./(2*s^2));


figure(1)

plot(t./60,S)
grid on