% osc_shm_01.m

% Simulation for the vertical motion of an object attached to a spring

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% 170517

clear all
close all
clc

   
figure(1)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.45 0.45]);
   set(gcf,'color','w');

% Model 1 -------------------------------------------------------------

% amplitude sMax  / period T / no. periods nT / no. time steps N
   sMax = 10;
   T = 100;
   nT = 3;
   N = 2000;

% angular frequency w / time t / displacement s
   w = 2*pi / T;
   tMax = nT * T;
   t = linspace(0,tMax,N);
   
   s = sMax .* cos(w*t);
   v = -(sMax*w) .* sin(w*t);
   a = -w^2 .* s;   
  
   subplot(3,3,1)
   xP = t; yP = s;
   plot(xP, yP,'b','linewidth',2);
   ylabel('s  [ m ]');
   axis([0 tMax -10 10]);
   set(gca,'yTick',-10:5:10);
   set(gca,'fontSize',14);
   grid on
   
   subplot(3,3,4)
   xP = t; yP = v;
   plot(xP, yP,'r','linewidth',2);
   ylabel('v  [ m.s^{-1} ]');
   axis([0 tMax -2 2]);
   set(gca,'yTick',-2:1:2);
   set(gca,'fontSize',14);
   grid on
      
   subplot(3,3,7)
   xP = t; yP = a;
   plot(xP, yP,'m','linewidth',2);
   ylabel('a  [ m.s^{-2} ]');
   xlabel('t  [ s ]');
   axis([0 tMax -0.4 0.4]);
   set(gca,'yLim',[-0.4 0.4]);
   set(gca,'yTick',-0.4:0.2:0.4);
   grid on
   set(gca,'fontSize',14);

% Model 2 -------------------------------------------------------------

% amplitude sMax  / period T / no. periods nT / no. time steps N
   sMax = 8;
   T = 50;
   
% angular frequency w / time t / displacement s
   w = 2*pi / T;
   s = sMax .* cos(w*t);
   v = -(sMax*w) .* sin(w*t);
   a = -w^2 .* s;   
  
   subplot(3,3,2)
   xP = t; yP = s;
   plot(xP, yP,'b','linewidth',2);
   ylabel('s  [ m ]');
   axis([0 tMax -10 10]);
   set(gca,'yTick',-10:5:10);
   set(gca,'fontSize',14);
   grid on
   
   subplot(3,3,5)
   xP = t; yP = v;
   plot(xP, yP,'r','linewidth',2);
   ylabel('v  [ m.s^{-1} ]');
   axis([0 tMax -2 2]);
    set(gca,'yTick',-2:1:2);
   set(gca,'fontSize',14);
   grid on
   
   subplot(3,3,8)
   xP = t; yP = a;
   plot(xP, yP,'m','linewidth',2);
   ylabel('a  [ m.s^{-2} ]');
   xlabel('t  [ s ]');
   axis([0 tMax -0.4 0.4]);
   set(gca,'yLim',[-0.4 0.4]);
   set(gca,'yTick',-0.4:0.2:0.4);
   set(gca,'fontSize',14);
   grid on
   

% Model 3 -------------------------------------------------------------

% amplitude sMax  / period T / no. periods nT / no. time steps N
   sMax = 5;
   T = 25;
   
% angular frequency w / time t / displacement s
   w = 2*pi / T;
   s = sMax .* cos(w*t);
   v = -(sMax*w) .* sin(w*t);
   a = -w^2 .* s;   
  
   subplot(3,3,3)
   xP = t; yP = s;
   plot(xP, yP,'b','linewidth',2);
   ylabel('s  [ m ]');
   axis([0 tMax -10 10]);
   set(gca,'yTick',-10:5:10);
   set(gca,'fontSize',14);
   grid on
   
   subplot(3,3,6)
   xP = t; yP = v;
   plot(xP, yP,'r','linewidth',2);
   ylabel('v  [ m.s^{-1} ]');
   axis([0 tMax -2 2]);
   set(gca,'yTick',-2:1:2);
   set(gca,'fontSize',14);
   grid on
   
   subplot(3,3,9)
   xP = t; yP = a;
   plot(xP, yP,'m','linewidth',2);
   ylabel('a  [ m.s^{-2} ]');
   xlabel('t  [ s ]');
   axis([0 tMax -0.4 0.4])
   set(gca,'fontSize',14);
   set(gca,'yLim',[-0.4 0.4]);
   set(gca,'yTick',-0.4:0.2:0.4);
   grid on
   