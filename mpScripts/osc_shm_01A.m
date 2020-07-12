% osc_shm_01A.m

% SIMPLE HARMONIC MOTION
% Simulation for the vertical motion of an object attached to a spring
% Animation of the motion can be saved as an animated gif

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% 170517

clear all
close all
clc


% =======================================================================
%    Setup for saving images (im) 
% ======================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_shm_170518.gif';
% Delay in seconds before displaying the next image  
   delay = 0.10;  
% Frame counter start
   nt = 1; 

% ======================================================================
%    PARAMETERS
% ======================================================================

% amplitude sMax  / period T / no. periods nT / no. time steps N
   sMax = 10;
   T = 100;
   nT = 3;
   N = 200; 
   
% angular frequency w / time t / displacement s
   w = 2*pi / T;
   tMax = nT * T;
   t = linspace(0,tMax,N);
   
   s = sMax .* cos(w*t);
   v = -(sMax*w) .* sin(w*t);
   a = -w^2 .* s;
   
   
% ======================================================================
%   GRAPHICS
% ======================================================================   
   pos = [0.1 0.2 0.6 0.6];
   pos1 = [0.01 0.1 0.2 0.8];
   pos2 = [0.3 0.7 0.6 0.25];
   pos3 = [0.3 0.4 0.6 0.25];
   pos4 = [0.3 0.1 0.6 0.25];
   
figure(1)
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = t; yP = s;
   subplot('Position',pos2) 
   plot(xP, yP,'b','lineWidth',2);
   hold on
   axis([0 tMax -11 11]);
   ylabel('s  [ m ]');
   set(gca,'fontSize',14);
   grid on
   
   xP = t; yP = v;
   subplot('Position',pos3)
   plot(xP, yP,'r','lineWidth',2);
   hold on
   axis([0 tMax -1.1*max(v) 1.1*max(v)]);
   ylabel('v  [ m.s^{-1} ]');
   set(gca,'fontSize',14);
   grid on
   
   xP = t; yP = a;
   subplot('Position',pos4)
   plot(xP, yP,'m','lineWidth',2);
   hold on
   axis([0 tMax -1.1*max(a) 1.1*max(a)]);
   ylabel('a  [ m.s^{-2} ]');
   set(gca,'fontSize',14);
   xlabel('time  t  [ s ]')
   grid on
   
   for c = 1: N
     subplot('Position',pos1) 
     xP = zeros(1,N);
     yP = s;
     hPlot = plot(xP(c),yP(c),'o');
     axis off
     set(hPlot,'MarkerSize',18,'MarkerEdgeColor',[0 0 0], ...
       'MarkerFaceColor', [0 0 0]);
     axis([-1 1 -11 11]);
     
     subplot('Position',pos2)
     xP = t(c); yP = s(c);
     plot(xP, yP,'bo');
     axis([0 tMax -11 11]);
     
     subplot('Position',pos3)
     xP = t(c); yP = v(c);
     plot(xP, yP,'ro');
     axis([0 tMax -1.1*max(v) 1.1*max(v)]);
     
     subplot('Position',pos4)
     xP = t(c); yP = a(c);
     plot(xP, yP,'mo');
     axis([0 tMax -1.1*max(a) 1.1*max(a)]);
     
     pause(0.01);
     
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
   
   



set(gca,'fontSize',14);


   