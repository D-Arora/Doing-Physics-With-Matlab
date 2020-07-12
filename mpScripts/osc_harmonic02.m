%osc_harmonic02.m

% Harmonic Oscillations
% Resonance response curve plot


% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm

clear all
close all
clc


% Inputs
m = 0.506602;   
k = 20;   
b = 2;    

Fmax = 1;

w0 = sqrt(k/m);
f0 = w0 / (2*pi);

nMax = 1000;

fMin = 0.5;
fMax = 1.5*f0;
df = (fMax-fMin)/(nMax-1);

f = fMin : df : fMax;

w = (2*pi) .* f;
A = zeros(nMax,1);

% Calculations

for n = 1 : nMax
   A(n) = Fmax/sqrt((k-m*w(n)^2)^2 + b^2*w(n)^2);
end

%plot
figure(1)
   pos = [0.07 0.05 0.32 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
plot(f/f0,A,'LineWidth',3)
set(gca,'xLim',[fMin fMax]);
ylabel('amplitude  A','FontSize',14)
xlabel('driving frequency  f_D / f_0','FontSize',14)
grid on
tm1 = 'f_0  =  ';
tm2 = num2str(f0,' %2.2f  Hz ');
tm3 = '     b  =  ';
tm4 = num2str(b,' %2.1f ');
tm = [tm1 tm2 tm3 tm4];
title(tm)
set(gca,'fontsize',14);
   