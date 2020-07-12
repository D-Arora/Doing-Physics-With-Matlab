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
% Mode numbers n
   n = [4 4];
   
   

% SETUP ===============================================================
   k = n.* (pi / (2*L));         % wave numbers
   x = linspace(-L,L,N);         % x grid
   wL = 2*pi ./ k;               % wavelength       
   K = 1i .* k;                 % complex wave number
   
% color mapping for phase phi using script ColorCode.m
   m = (380-660) / (2*pi);
   b = 380 - m*pi;
   A = 1/sqrt(2*L);
  
   
% CALCULATIONS ========================================================
     u1 = A .* exp(K(1).*x);                    % wavefunction
     u2 = A .* exp(K(2).*x); 
     An = simpson1d(conj(u1).*u2,-L,L);   % check normalization
     
% GRAPHICS  ============================================================

figure(1)   % ---------------------------------------------------------
  pos = [0.1 0.1 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  thisColorMap = hsv(512);
 
subplot(3,1,1)   
     hold on
     for cx = 1:N-1
     phase = angle(u1(cx));
     CwL = m * phase + b;
     thisColor = ColorCode(CwL*1e-9);
     xP = x; yP = A.*ones(1,N);
     h_area = area(xP(cx:cx+1),yP(cx:cx+1));
     set(h_area,'FaceColor',thisColor);
     set(h_area,'EdgeColor',thisColor); 
     set(gca,'xLim',[-L L]);
     xlabel('x','fontsize',14);
     ylabel('| u |','fontsize',14);
     tm1 = 'mode number n = ';
     tm2 = num2str(n(1),0);
     tm = [tm1 tm2];
     title(tm)
     set(gca,'fontsize',14);
     set(gca,'xtick',-L:10:L)
     end
     
 subplot(3,1,2)
     hold on
     for cx = 1:N-1
     phase = angle(u2(cx));
     CwL = m * phase + b;
     thisColor = ColorCode(CwL*1e-9);
     xP = x; yP = A.*ones(1,N);
     h_area = area(xP(cx:cx+1),yP(cx:cx+1));
     set(h_area,'FaceColor',thisColor);
     set(h_area,'EdgeColor',thisColor); 
     set(gca,'xLim',[-L L]);
     xlabel('x','fontsize',14);
     ylabel('| u |','fontsize',14);
     tm1 = 'mode number n = ';
     tm2 = num2str(n(2),0);
     tm = [tm1 tm2];
     title(tm)
     set(gca,'fontsize',14);
     set(gca,'xtick',-L:10:L)
     end    

subplot(3,1,3)
     hold on
     for cx = 1:N-1
     phase = angle(conj(u1(cx)) .* u2(cx));
     CwL = m * phase + b;
     thisColor = ColorCode(CwL*1e-9);
     xP = x; yP = conj(u1).*u2;
     h_area = area(xP(cx:cx+1),yP(cx:cx+1));
     set(h_area,'FaceColor',thisColor);
     set(h_area,'EdgeColor',thisColor); 
     set(gca,'xLim',[-L L]);
     xlabel('x','fontsize',14);
     ylabel('| u |','fontsize',14);
     tm1 = 'mode number n = ';
     tm2 = num2str(n(2),0);
     tm = [tm1 tm2];
     title(tm)
     set(gca,'fontsize',14);
     set(gca,'xtick',-L:10:L)
     end    

  
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
