% mathComplex4A.m

% A Visual Approach to the Fourier Transform
%   but what is a Fourier Transform?

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/mathComplex4.htm

% This script was inspiration from the youtube video �But what is the Fourier Transform?
%   A visual Introduction� by 3Blue1Brown
%   https://www.youtube.com/watch?v=spUNpyF58BY
%   Matlab programs were developed to create similar figures
%   and to connect the ideas presented in the video.

% For different simulation you may need to make small chnages to the
% script.

% Matlab Version 2018a / 180917

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% John A Sims 
% email: john.sims@ufabc.edu.br
% Biomedical Engineering Department
% Federal University of ABC   Sao Bernardo Campus Brasil

clear
close all
clc


% INPUT SECTION =======================================================

% Winding frequency  [1 Hz]
fw = 3;
% Simulation time    [1/fw 2/fw 3/fw .....]
tMax = 50;
% Number of calculation  [5000]
N = 5001;
% Signal period  [1]
Ts = 1;
% Default amplitudes  [1]
A = 1; R = A;


% CALCULATION SECTION =================================================
fs = 1/Ts;
w = 2*pi*fs;
t = linspace(0,tMax,N);

% COMPLEX FUNCTIONS 
% Signal function h(t) *************************************************

%  h = R.*sin(2*pi*fs*t);
%  h = R.*sin(2*pi*fs*t) + 0.5;
  h = 1.0.*sin(2*pi*fs*t) + 0.5.*sin(2*pi*2*fs*t) + 0.25.*sin(2*pi*3*fs*t);
%  h = ( 0.8.*sin(2*pi*fs*t) + 0.4.*sin(2*pi*3*fs*t) );
% ***********************************************************************
%nfw = 22;
%fw = li5space(0.2,4.4,nfw);
fwMax = 5;
fwMin = 0;
nfw = 201;
fw = linspace(fwMin,fwMax,nfw);
H = zeros(1,nfw);
%H = complex(zeros(Nfw,1));
for cc = 1 : nfw
  %t = linspace(0,6/fw(cc),N); 
% Complex exponetial function r(t)
   r = exp(-1i*2*pi*fw(cc)*t);
% Complex function s(t)
   u = h .* r;

% spectal component
H(cc) = abs(sum(u));
end


% GRAPHICS SECTION ====================================================
  
figure(1)   % Re & Im vs t 
  pos = [0.05 0.05 0.30 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
   xP = fw; yP = H./max(H);
   plot(xP,yP,'b','linewidth',2)

  set(gca,'fontsize',14)
  xlabel('winding frequency  f_w  [Hz]') 
  ylabel('| H |')
  grid on
  grid minor
  