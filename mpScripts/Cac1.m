% Cac1.m

% ac circuit modeling and analysis 
% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180107 


clear all
close all
clc

% ========================================================================
%   INPUTS   default values [ ]
% =======================================================================

% Voltage #1 and #2 [10 8] [1000 1000]   [0 pi/2] 
   V = [10 8];
   f = [1000 1000];
   phi = [0 -pi/4];
% time
   N = 1000;
   T = 1/f(1);
   tMax = 3*T;
   t = linspace(0,tMax,N);
   w = (2*pi).*f;

% =======================================================================  
%   VOLTAGE CALCULATIONS
% =======================================================================   
    v1 = V(1) .* exp(1i*(w(1)*t + phi(1)));
    v2 = V(2) .* exp(1i*(w(2)*t + phi(2)));


% =======================================================================
%   GRAPHICS
% =======================================================================
   fss = 16;
figure(1)   
   pos = [0.1 0.1 0.50 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   pos1 = [0.47 0.15 0.5 0.7];
   ax1 = subplot('Position',pos1);
   
   xP = t.*1e3; yP = real(v1);        
   plot(xP, yP,'b','linewidth',2);   
   hold on
   yP = real(v2);                     
   plot(xP, yP,'r','linewidth',2);   
  
   grid on
   xlabel('time  t  [ ms ]')
   ylabel('voltage  V  [ V ]');
     tm1 = 'Phases  \phi / \pi   ';
     tm2 = '\phi_1  =  ';
     tm3 = num2str(phi(1)/pi,2);
     tm4 = '    \phi_2  =  ';
     tm5 = num2str(phi(2)/pi,2);
     tm = [tm1 tm2 tm3 tm4 tm5];
   title(tm);
   set(gca,'fontsize',fss)
   
   pos2 = [0.00 0.28 0.45 0.45];
   ax2 = subplot('Position',pos2);
   
   A = linspace(0,2*pi,1000);   % plot circle
   x = cos(A); y = sin(A);
   plot(x,y,'k');
   hold on
   x = [-1 1]; y = [0 0];       % plot axes
   plot(x,y,'k');
   y = [-1 1]; x = [0 0];
   plot(x,y,'k');
   
% v1 phasor   
   x = cos(phi(1)); y = sin(phi(1));
   xP = [0 x]; yP = [y y];
   plot(xP, yP,'b','linewidth',1);
   xP = [0 x]; yP = [0 y];
   plot(xP,yP,'b','linewidth',3);
%    xP = [x x]; yP = [0 y];
%    plot(xP, yP,'b','linewidth',1); 
%   xP = [0 x]; yP = [0 0];
%   plot(xP, yP,'b','linewidth',3)

   % v2 phasor  
    x = cos(phi(2)); y = sin(phi(2));
    xP = [0 x]; yP = [y y];
    plot(xP, yP,'r','linewidth',1);
    xP = [0 x]; yP = [0 y];
    plot(xP,yP,'r','linewidth',3);
     
%   xP = [0 x]; yP = [0 0];
%   plot(xP, yP,'r','linewidth',3)
   
   axis square
   axis off
   tm = [tm2 tm3 tm4 tm5];
   text(-1,-1.5,tm,'fontsize',16);