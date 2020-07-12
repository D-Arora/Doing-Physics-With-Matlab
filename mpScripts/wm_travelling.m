% wm_travelling.m

% TRAVELLING WAVES
% Animation of the motion can be saved as an animated gif
%   Travelling wave / Transverse wave / Longitudinal wave

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% 170517

clear all
close all
clc

% Type of animation
   flagC = 3;
             % flagC = 1   animation of a travelling wave
             % flagC = 2   animation of a transverse wave
             % flagC = 3   animation of a longitudinal wave

% =======================================================================
%    Setup for saving images (im) 
% ======================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_wm_170519C.gif';
% Delay in seconds before displaying the next image  
   delay = 0.20;  
% Frame counter start
   nt = 1; 

switch flagC
    
case 1   % ------------------------------------------------------------
        
% amplitude A  / wavelength L / period T
% no. periods nT / no. wavelengths steps nX
% Number of calcutions Nt Nx
   A = 10;
   L = 20;
   T = 25;
   nT = 3; nX = 3;
   Nx = 200; 
   Nt = 200;
   
% left (-1) or right (+1)
   e = -1;
% angular frequency w / time t / wave number k /displacement s
   w = 2*pi / T;
   k = 2*pi / L;
   tMax = nT * T;
   xMax = nX * L;
   
   t = linspace(0,tMax,Nt);
   x = linspace(0,xMax,Nx);
   
% wave function s(x,t)   
   s = zeros(Nx,Nt);
   
   for c = 1 : Nt
     s(:,c) = A .* sin(k.*x + e*w*t(c));
   end
     v = L / T;
     
    
figure(1)
   pos = [0.1 0.1 0.4 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   t1 = 't  =  ';
   t3 = '   s';
   
for c = 1 : Nt
   hold off
   xP = x; yP = s(:,c);
   plot(xP, yP,'b','lineWidth',2);
   hold on
   xP = x(100); yP = s(100,c);
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
   xP = v * t(c); yP = 0;
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',6,'MarkerEdgeColor',[1 0 1], ...
       'MarkerFaceColor', [1 0 1]);
   
   axis([0 xMax -11 11]);
   ylabel('wave function  s ');
   xlabel('position  x  [ m ]')
   set(gca,'fontSize',14);
   grid on
   t2 = num2str(t(c), '%2.1f \n');
   tm = [t1 t2 t3];
   title(tm);
   pause(0.2)
   
   if f_gif > 0 
       frame = getframe(1);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
     %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
   end
   
end

case 2  % --------------------------------------------------------------
   % amplitude A  / wavelength L / period T
% no. periods nT / no. wavelengths steps nX
% Number of calcutions Nt Nx
   A = 10;
   L = 20;
   T = 25;
   nT = 3; nX = 3;
   Nx = 2000; 
   Nt = 200;
   
% left (-1) or right (+1)
   e = -1;
% angular frequency w / time t / wave number k /displacement s
   w = 2*pi / T;
   k = 2*pi / L;
   tMax = nT * T;
   xMax = nX * L;
   
   t = linspace(0,tMax,Nt);
   x = linspace(0,xMax,Nx);
   
   % wave function s(x,t)   
   s = zeros(Nx,Nt);
   
   for c = 1 : Nt
     s(:,c) = A .* sin(k.*x + e*w*t(c));
   end
    
figure(2) 
   pos = [0.1 0.1 0.4 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   NS = 20; 
for c = 1 : Nt
   hold off
   xP = x(1:NS:end); yP = s((1:NS:end),c);
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[0 0 1], ...
       'MarkerFaceColor', [0 0 1]);
   hold on
   xP = x(1:100:end); yP = s(1:100:end,c);
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
   axis([0 xMax -11 11]);
   axis off
   pause(0.1)
   
   if f_gif > 0 
       frame = getframe(2);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
     %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
   end
end

case 3   % LONGITUDINAL WAVE ---------------------------------------------

% x positions of particles
   x = 0:4:100;
   y = zeros(length(x),1);
   
% amplitude 
   A = 1.2;
% period and wavelength
   L = 20;
   k = 2*pi / L;
   T = 25;
   w = 2*pi / T;
% number of time steps
   Nt = 100;
   t = linspace(0,3*T,Nt);

   
figure(3) 
pos = [0.1 0.1 0.4 0.1];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
for c = 1:Nt
    hold off
    xC = x + A.* sin(k*x - w*t(c));
    xP = xC; yP = y;
    hPlot = plot(xP,yP,'o');
    set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
   hold on
   xP = xP(12); yP = y(12);
    hPlot = plot(xP,yP,'o');
    set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 1], ...
       'MarkerFaceColor', [1 0 1]);
    axis([-10 110 -1 1]);
    
    axis off;
    pause(0.1);
    
    if f_gif > 0 
       frame = getframe(3);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
     %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
   end
    
    
end    
    
end   % switch
   


     
     


   