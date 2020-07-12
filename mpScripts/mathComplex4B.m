% mathComplex4B.m

% A Visual Approach to the Fourier Transform
%   But what is a Fourier Transform?
% Calculation of the frequency spectrum of a signal

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/mathComplex4.htm

% This script was inspirated from the youtube video ”But what is the Fourier Transform?
%   A Visual Introduction” by 3Blue1Brown
%   https://www.youtube.com/watch?v=spUNpyF58BY
%   Matlab programs were developed to create similar figures
%   and to connect the ideas presented in the video.

% For different simulation you may need to make small changes to the
% script.

% Matlab Version 2018a / 180919

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

% Winding frequency fw:  min / max / # elements 
fwMin = 0;
fwMax = 5;
nfw = 201;

% Simulation time t:   tMax / # elements 
tMax = 50;
N = 5001;

% Signal: frequency  / amplitude 
fs = 1;
R = 1;


% CALCULATION SECTION =================================================
% Signal
  Ts = 1/fs;
  ws = 2*pi*fs;
% Simulation time interval
  t = linspace(0,tMax,N);
% Winding frequency
  fw = linspace(fwMin,fwMax,nfw);
% Componenets of frequency spectrum  
  H = zeros(1,nfw);

% Signal function h(t) *************************************************
%  h = R.*sin(2*pi*fs*t);
%  h = R.*sin(2*pi*fs*t) + 0.5;
   h = 1 + 1.0.*sin(2*pi*fs*t) + 0.5.*sin(2*pi*2*fs*t) + 0.25.*sin(2*pi*3*fs*t);
%  h = ( 0.8.*sin(2*pi*fs*t) + 0.4.*sin(2*pi*3*fs*t) );
% ***********************************************************************

% Spectrun
for cc = 1 : nfw
   % Complex exponetial function r(t)
      r = exp(-1i*2*pi*fw(cc)*t);
   % Complex function u(t)
      u = h .* r;
   % Spectal Component
      H(cc) = abs(sum(u));
end


% GRAPHICS SECTION ====================================================
  
figure(1)   
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
  