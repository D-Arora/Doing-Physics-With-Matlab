% qmElectromBeamSF.m


close all
clc
clear


X = 1e-9;
x = X.*linspace(-0.5,0.5,501);
b = 5e10;
y = 400./(1 + exp(-b.*x));

% x = -0.5e-9;
% y = 0.01;
% 
% b = -log(1/y-1)/x



figure(1)
plot(x,y)


