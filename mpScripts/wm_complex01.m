% wm_complex01.m

% Animation of a travelling wave
% Visualization of complex numbers using colors to show the phase

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm

clear all
close all
clc

% INPUTS ==============================================================
% Number of grid points
   NX = 201;
% Number of time steps
   NT = 9;
% Wavelength / Number of Wavelengths
   L = 20;
   nL = 3;
   
% Period / Number of periods
   T = 8;
   nT = 5;

% SETUP ===============================================================
   k = 2*pi/L;               % propgation constant
   f = 1/T;                  % frequency
   w = 2*pi*f;               % angular frequency
   xMax = nL*L;              % Max X value
   x = linspace(0,nL*L,NX);  % x grid
   tMax = nT*T;              % Max time interval
   %t = linspace(0,tMax,NT);  % time steps
   t = (0:8) .* (T/8);
   KX = k.*x;                 

% color mapping for phase phi using script ColorCode.m
   m = (380-660) / (2*pi);
   b = 380 - m*pi;
   
   
% CALCULATIONS ========================================================
   psi = zeros(NT,NX);      % wavefunction
   
   for c = 1 : NT
      psi(c,:) = exp(-1i.*(KX - w.*t(c)));
   end

% GRAPHICS  ============================================================
    figure(1)
    pos = [0.1 0.1 0.3 0.2];
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
         
   %%
 figure(2)
   pos = [0.5 0.1 0.3 0.15];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
     
  % h_area = area(x,real(psi(1,:)));
  % set(h_area,'FaceColor',[1 1 0]);
  % set(h_area,'EdgeColor','none');
   
   for c = 9
   thisColorMap = hsv(256);
      for cx = 1:NX-1
         phase = angle(psi(c,cx));
         wL = (m * phase + b)*1e-9;
         thisColor = ColorCode(wL);
         h_area = area(x(cx:cx+1),real(psi(c,cx:cx+1)));
         set(h_area,'FaceColor',thisColor);
         set(h_area,'EdgeColor',thisColor);
         set(gca,'xLim',[0 xMax]);
         set(gca,'yLim',[-1 1]);
         hold on
      end
         tm1 = 't = ';
         tm2 = num2str(t(c),0);
         tm3 = ' / 8';
         tm = [tm1 tm2 tm3];
         title(tm);
         pause(0.1)
         hold off
   end
   
    axis off 
   
 