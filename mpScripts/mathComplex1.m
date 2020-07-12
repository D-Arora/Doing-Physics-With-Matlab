% mathComplex1.m

% Visualization of a complex valued function
% The phase of the complex function is displayed using colors

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Matlab Version 2018a / 180912

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% John A Sims 
% email: john.sims@ufabc.edu.br
% Biomedical Engineering Department
% Federal University of ABC   Sao Bernardo Campus Brasil

% Script calls the function  ColorCode.m

clear all
close all
clc

tic

% INPUTS ==============================================================

% Number of grid points  [201]
   N = 201;
% Period   [1]
   T = 1;
  % Number of periods displayed [3]
   nT = 3;
% amplitude  [10]
   R = 10;
   
   
% CALCULATIONS ========================================================
   
w = 2*pi / T;             % angular frequency
   f = 1/T;                  % frequency
   tMax = nT*T;              % time interval
   t = linspace(0,tMax,N);   % t grid
   W = 1i*w;                 % complex angular frequencywave number
   
   u = R .* exp(W.*t);                    % wavefunction

   % Color mapping for phase phi using script ColorCode.m
   m = (380-660) / (2*pi);
   b = 380 - m*pi;
  
     
% GRAPHICS  ============================================================

figure(1)    %---------------------------------------------------------
  pos = [0.1 0.05 0.3 0.25];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = t; yP = real(u);
   plot(xP,yP,'b','linewidth',2);
  hold on
  yP = imag(u);
    plot(xP,yP,'r','linewidth',2)
  yP = abs(u);
    plot(xP,yP,'k','linewidth',2)
  ylabel('u','fontsize',14)
  xlabel('t','fontsize',14)
  legend('real','imaginary','magnitude','location','northoutside','orientation','horizontal')
  box on
  grid on
  set(gca,'fontsize',14)
  set(gca,'ylim',[-R-1 R+1])

figure(2)   % ---------------------------------------------------------
  pos = [0.1 0.4 0.4 0.4];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
 
  xP = t; yP = real(u); zP = imag(u); 
  plot3(xP,yP,zP,'k','linewidth',2);
 
 hold on
  xP = t; yP = real(u); zP = -abs(u).*ones(1,length(t)); 
  %zP = zeros(1,length(t));        % set off-set to zero
  plot3(xP,yP,zP,'b','linewidth',2);
  xP = t; zP = imag(u); yP = abs(u).*ones(1,length(t)); 
  %yP = zeros(1,length(t));        % set off-set to zero
  plot3(xP,yP,zP,'r','linewidth',2);
  
  xlabel('t','color','k','fontsize',16)
  ylabel('real(u)','color','b','fontsize',16)
  zlabel('imag(u)       ','color','r','fontsize',16,'Rotation',0)
  xP = t(1); zP = imag(u(1)); yP = real(u(1)); 
  Hplot = plot3(xP,yP,zP,'ko');
  set(Hplot,'markersize',6,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0])
  box on
  grid on
  set(gca,'fontsize',14)
  
figure(3)   % ---------------------------------------------------------
  pos = [0.42 0.1 0.3 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
   hold on
   
 % thisColorMap = hsv(512);
   for cx = 1:N-1
     phase = angle(u(cx));
     CwL = m * phase + b;
     thisColor = ColorCode(CwL*1e-9);
     xP = t; yP = R.*ones(1,N);
     h_area = area(xP(cx:cx+1),yP(cx:cx+1));
     set(h_area,'FaceColor',thisColor);
     set(h_area,'EdgeColor',thisColor); 
     
     xlabel('t','fontsize',14);
     ylabel('| u |','fontsize',14);
     set(gca,'fontsize',14);
   end
       
figure(4)   % ---------------------------------------------------------
  pos = [0.31 0.4 0.3 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = t; yP = angle(u)./pi;
    plot(xP,yP,'b','linewidth',2);
  box on
  grid on
  xlabel('t','fontsize',14);
  ylabel('phase / \pi','fontsize',14); 
  set(gca,'fontsize',14)
  set(gca,'yLim',[-1.2 1.2])
  set(gca,'ytick',-1:0.5:1)
  
figure(5)   %----------------------------------------------------------
 pos = [0.63 0.4 0.3 0.2];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    set(gcf,'color','w');
    N = 512;
    xP = linspace(-pi,pi,N);
    yP = ones(1,length(xP));
    
    hold on
    
    thisColorMap = hsv(512);
       for cx = 1:N-1
         phase = xP(cx);
         wL = m * phase + b;
         thisColor = ColorCode(wL*1e-9);
         h_area = area(xP(cx:cx+1)/pi,yP(cx:cx+1));
         set(h_area,'FaceColor',thisColor);
         set(h_area,'EdgeColor',thisColor); 
         set(gca,'xLim',[-1 1]);
       end
         xlabel('phase  [ \pi rad]','fontsize',14);
         set(gca,'fontsize',14);
         set(gca,'yTick',[]);
         set(gca,'xtick',-1:0.25:1)


toc
