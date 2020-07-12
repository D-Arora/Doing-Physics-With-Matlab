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
   L = 50;
% Mode number n
   n = 4;
   
   

% SETUP ===============================================================
   k = n*pi / (2*L);         % wave number
   x = linspace(-L,L,N);     % x grid
   wL = 2*pi/k;              % wavelength       
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
  pos = [0.1 0.1 0.3 0.2];
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
  legend('real','imaginary','magnitude','location','best')
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
  xP = real(psi); zP = imag(psi); yP = x;
   plot3(xP,yP,zP,'linewidth',2);
  
  ylabel('x')
  xlabel('real(u)','fontsize',14)
  zlabel('imag(u)','fontsize',14)
  
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
