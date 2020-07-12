% mathComplex4B.m

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

% Winding frequency  [min = 0.5  max = 4 Hz Nfw = 1000]
fwMin = 0.2; fwMax = 4.4; Nfw = 22;
% Winding frequency for plot [1]
fwPlot = 3;
% Simulation time    [ 3 s]
tMax = 50;
% Number of calculation  [5000]
N = 5001;
% Signal period  [1]
T = 1;

% CALCULATION SECTION =================================================
f = 1/T;
w = 2*pi*f;
t = linspace(0,tMax,N);
fw = linspace(fwMin,fwMax,Nfw);

% **********************************************************************
% Signal function h(t)
%  h = A.*sin(2*pi*f*t);
   h = ( 0.8.*sin(2*pi*f*t) + 0.4.*sin(2*pi*3*f*t) );

% ***********************************************************************

%  index for winding frequencty phase plot
nfw = find(fw > fwPlot,1)-1;

%
 s = complex(zeros(Nfw,N)); sumS = complex(zeros(Nfw,1));
for cc = 1 : Nfw
 s(cc,:) = h.* exp(-1i*2*pi*fw(cc)*t);
 sumS(cc) = sum(s(cc,:));
end

% for f = 1:numel(wf)
%     compS(f,:) = s.*exp(-2*pi*1i*wf(f)*t);
%     SumcompS(f) = sum(compS(f,:));

 % GRAPHICS SECTION ====================================================

figure(1)   
  pos = [0.05 0.05 0.30 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = fw; yP = abs(sumS)./max(abs(sumS));
  plot(xP,yP,'b','linewidth',2)
    
  grid on
  box on
  set(gca,'fontsize',14)
  xlabel('winding frequency  f_w  [Hz]') 
  ylabel('sumS)')



figure(2)   % phase plot
  pos = [0.05 0.35 0.25 0.30];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = real(s(nfw,:)); yP = -1.*ones(1,numel(xP));
    plot(xP,yP,'b','linewidth',2)
  hold on
  yP = imag(s(nfw,:)); xP = -1.*ones(1,numel(xP));
    plot(xP,yP,'r','linewidth',2)
   xP = real(s(nfw,:)); yP = imag(s(nfw,:));  
      plot(xP,yP,'k','linewidth',2)
%       
  xlim([-1.1 1.1])
  ylim([-1.1 1.1])
  set(gca,'xtick',-1:0.5:1)
  set(gca,'ytick',-1:0.5:1)
  grid on
  box on
  set(gca,'fontsize',14)
  xlabel('real(s)','color','b') 
  ylabel('imag(s)','color','r') 
  axis square 
  
%   tm1 = 'f = ';       tm2 = num2str(f, '%2.1f'); tm3 = '  Hz';
%   tm4 = '    f_w = ';  tm5 = num2str(fw,'%2.1f');
%   tm6 = '  Hz';
%   tm = [tm1 tm2 tm3 tm4 tm5 tm6];
%   title(tm,'fontweight','normal')
% % 
% figure(3)   %  [3D] view of complex function 
%   pos = [0.6 0.4 0.3 0.3];
%   set(gcf,'Units','normalized');
%   set(gcf,'Position',pos);
%   set(gcf,'color','w');
%  
%    xP = t; yP = real(s); zP = imag(s); 
%    plot3(xP,yP,zP,'k','linewidth',3);   
%   
%   hold on
%   xP = t; yP = real(s); zP = -R.*ones(1,length(t)); 
%   %zP = zeros(1,length(t));        % set off-set to zero
%   plot3(xP,yP,zP,'b','linewidth',2);
%   
%   xP = t; zP = imag(s); yP = R.*ones(1,length(t)); 
%   %yP = zeros(1,length(t));        % set off-set to zero
%   plot3(xP,yP,zP,'r','linewidth',2);
%   
%   xP = t(1); zP = imag(s(1)); yP = real(s(1)); 
%   Hplot = plot3(xP,yP,zP,'ko');
%   set(Hplot,'markersize',6,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0])
%   
%   xlabel('t','color','k','fontsize',16)
%   ylabel('real(s)','color','b','fontsize',16)
%   zlabel('imag(s)       ','color','r','fontsize',16,'Rotation',0)
% 
%   zlim([-1.1 1.1])
%   ylim([-1.1 1.1])
%   set(gca,'xtick',0:1:max(t))
%   set(gca,'ytick',-1:0.5:1)
%   set(gca,'ztick',-1:0.5:1)
%   box on
%   grid on
%   set(gca,'fontsize',14)
