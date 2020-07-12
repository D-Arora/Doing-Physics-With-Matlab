% mathComplex4A.m

% A Visual Approach to the Fourier Transform
%   but what is a Fourier Transform?

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
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
fw = 0.8;
% Simulation time    [1/fw 2/fw 3/fw .....]
tMax = 4/fw;
% Number of calculation  [5000]
N = 500;
% Signal period  [1]
Ts = 1;
% Default amplitudes  [1]
A = 1; R = A;


% CALCULATION SECTION =================================================
   fs = 1/Ts;
   w = 2*pi*fs;
   t = linspace(0,tMax,N);

% COMPLEX FUNCTIONS 

% Enter Signal function h(t) *******************************************

  h = R.*sin(2*pi*fs*t);
%  h = R.*sin(2*pi*fs*t) + 0.5;
%   h = 1.0.*sin(2*pi*fs*t) + 0.5.*sin(2*pi*2*fs*t) + 0.25.*sin(2*pi*3*fs*t);
%  h = ( 0.8.*sin(2*pi*fs*t) + 0.3.*sin(2*pi*3*fs*t) );
% ***********************************************************************

% Complex exponetial function r(t)
   r = exp(-1i*2*pi*fw*t);
% Complex function s(t)
   u = h .* r;

% spectal component
   H = sum(u);

% GRAPHICS SECTION ====================================================
  
 L = max(abs(u));    % Limits for plots
 
figure(1)   % Re & Im vs t 
  pos = [0.05 0.05 0.30 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = t; yP = real(u);
  plot(xP,yP,'b','linewidth',2)
  
  hold on
  
  xP = t; yP = imag(u);
  plot(xP,yP,'r','linewidth',2)
  
 % ylim([-1.1 1.1])
  legend('Re(u)','Im(u)','location','northoutside','orientation','horizontal')
  grid on
  box on
  set(gca,'fontsize',14)
  xlabel('time  t  [s]') 
  ylabel('u(t)')


figure(2)   % phase plot
  pos = [0.05 0.35 0.30 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  
  xP = real(u); yP = -1.*ones(1,numel(u));
    plot(xP,yP,'b','linewidth',2)
  hold on
  yP = imag(u); xP = -1.*ones(1,numel(u));
    plot(xP,yP,'r','linewidth',2)
  xP = real(u); yP = imag(u);  
    plot(xP,yP,'k','linewidth',2)
  xP = real(mean(u)); yP = imag(mean(u));
    Hplot = plot(xP,yP,'o');
    set(Hplot,'markersize',8,'markerfacecolor',[1 0 1],'markeredgecolor',[1 0 1])
  
  xlim([-1.1*L 1.1*L])
  ylim([-1.1*L 1.4*L])
 % set(gca,'xtick',-L:0.5:L)
 % set(gca,'ytick',-L:0.5:L)
  grid on
  box on
  set(gca,'fontsize',14)
  xlabel('real(u)','color','b') 
  ylabel('imag(u)','color','r') 
  axis equal 
  
  tm1 = 'f_s = ';       tm2 = num2str(fs, '%2.2f'); tm3 = '  Hz';
  tm4 = '    f_w = ';  tm5 = num2str(fw,'%2.2f');
  tm6 = '  Hz';
  tm = [tm1 tm2 tm3 tm4 tm5 tm6];
  title(tm,'fontweight','normal')
  
  tm1 = '  H(f_w)  =  ';
  tm2 = num2str(H,'%3.1f');
 % tm = [tm1 tm2];
 % text(-L, 1.25*L, tm,'fontsize',14)
  tm3 = '      abs( H(f_w) )  =  ';
  tm4 = num2str(abs(H),'%3.1f');
  tm = [tm1 tm2 tm3 tm4];
  text(-1.2*L, 1.2*L, tm,'fontsize',14)
 
  
figure(3)   %  [3D] view of complex function 
  pos = [0.6 0.4 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
 
   xP = t; yP = real(u); zP = imag(u); 
   plot3(xP,yP,zP,'k','linewidth',3);   
  
  hold on
  xP = t; yP = real(u); zP = -R.*ones(1,length(t)); 
  %zP = zeros(1,length(t));        % set off-set to zero
  plot3(xP,yP,zP,'b','linewidth',2);
  
  xP = t; zP = imag(u); yP = R.*ones(1,length(t)); 
  %yP = zeros(1,length(t));        % set off-set to zero
  plot3(xP,yP,zP,'r','linewidth',2);
  
  xP = t(1); zP = imag(u(1)); yP = real(u(1)); 
  Hplot = plot3(xP,yP,zP,'ko');
  set(Hplot,'markersize',6,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0])
  
  xlabel('t','color','k','fontsize',16)
  ylabel('real(u)','color','b','fontsize',16)
  zlabel('imag(u)       ','color','r','fontsize',16,'Rotation',90)

%  zlim([-1.1*L 1.1*L])
%  ylim([-1.1*L 1.1*L])
  set(gca,'xtick',0:1:max(t))
  set(gca,'ytick',-1:0.5:1)
  set(gca,'ztick',-1:0.5:1)
  box on
  grid on
  set(gca,'fontsize',14)
