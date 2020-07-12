% vqm001.m

% Complex valued functions
% Visualization of complex numbers using colors to show the phase
% Plane wave representation

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

clear all
close all
clc

tic

% INPUTS ==============================================================
% Number of grid points
   N = 201;
% Range L
   L = 60;
% Mode number n
   n = 4;
   
   

% SETUP ===============================================================
   wL = 4*L/n;
   k = 2*pi / wL;         % wave number
   x = linspace(-L,L,N);     % x grid
   K = 1i*k;                 % complex wave number
   
% color mapping for phase phi using script ColorCode.m
   m = (380-660) / (2*pi);
   b = 380 - m*pi;
   A = 1/sqrt(2*L);
  
   
% CALCULATIONS ========================================================
     psi = A .* exp(K.*x);                    % wavefunction
      
     An = simpson1d(conj(psi).*psi,-L,L);   % check normalization
     
% GRAPHICS  ============================================================

figure(1)    %---------------------------------------------------------
  pos = [0.1 0.1 0.3 0.25];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = x; yP = real(psi);
   plot(xP,yP,'b','linewidth',2);
  hold on
  yP = imag(psi);
    plot(xP,yP,'r','linewidth',2)
  yP = abs(psi);
    plot(xP,yP,'k','linewidth',2)
  ylabel('u','fontsize',14)
  xlabel('x','fontsize',14)
  legend('real','imaginary','magnitude','location','northoutside','orientation','horizontal')
  box on
  grid on
  set(gca,'fontsize',14)
  set(gca,'xtick',-L:10:L)
  set(gca,'ylim',[-0.12 0.12])

%%  
figure(2)   % ---------------------------------------------------------
  pos = [0.1 0.4 0.5 0.5];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
 
  xP = x; yP = real(psi); zP = imag(psi); 
  plot3(xP,yP,zP,'k','linewidth',2);
 
 hold on
  xP = x; yP = real(psi); zP = -abs(psi).*ones(1,length(x)); 
  plot3(xP,yP,zP,'b','linewidth',2);
  xP = x; zP = imag(psi); yP = abs(psi).*ones(1,length(x)); 
  plot3(xP,yP,zP,'r','linewidth',2);
  xlabel('x','color','k','fontsize',16)
  ylabel('real(u)','color','b','fontsize',16)
  zlabel('imag(u)  ','color','r','fontsize',16,'Rotation',0)
  xP = x(1); zP = imag(psi(1)); yP = real(psi(1)); 
  Hplot = plot3(xP,yP,zP,'ko');
  set(Hplot,'markersize',6,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0])
  set(gca,'xtick',-L:20:L);
  box on
  grid on
  set(gca,'fontsize',14)
%%  
figure(3)   % ---------------------------------------------------------
  pos = [0.42 0.1 0.3 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
   hold on
   
  thisColorMap = hsv(512);
     for cx = 1:N-1
     phase = angle(psi(cx));
     CwL = m * phase + b;
     thisColor = ColorCode(CwL*1e-9);
     xP = x; yP = A.*ones(1,N);
     h_area = area(xP(cx:cx+1),yP(cx:cx+1));
     set(h_area,'FaceColor',thisColor);
     set(h_area,'EdgeColor',thisColor); 
     set(gca,'xLim',[-L L]);
     xlabel('x','fontsize',14);
     ylabel('| u |','fontsize',14);
     set(gca,'fontsize',14);
     set(gca,'xtick',-L:10:L)
     end
     
    
figure(4)   % ---------------------------------------------------------
  pos = [0.31 0.4 0.3 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = x; yP = angle(psi)./pi;
    plot(xP,yP,'b','linewidth',2);
  box on
  grid on
  xlabel('x','fontsize',14);
  ylabel('phase / \pi','fontsize',14); 
  set(gca,'fontsize',14)
  set(gca,'xtick',-L:10:L)

  
%%  
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
