% mpSR03.m
% 190618  Matlab R2018b

% SPECIAL RELATIVITY ANIMATIONS
% Motion is relative;
%  A ball is thrown into the air.
%  Is the path of the ball a straight line or a parabola?
% The answer depends upon the observer

% DOING PHYSICS ONLINE 
%    ../mphome.htm
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ONLINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Documentation
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod72/mod72A.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

close all
clear
clc

% ANIMATION and GIF file ==============================================

% file name for animated gif   
    ag_name = 'ag_SR.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.1;
    frame1 = 0;

% SETUP ===============================================================

% Acceleration due to gravity
  ay = -10;
% Initial velocities
  v1x = 5;
  v1y = 20;
% Number of time steps
  nT = 51;

  
% CALCULATIONS  =======================================================

% Time
  tMax = abs(2*v1y/ay);
  t = linspace(0,tMax, nT);
% Displacements
  x = v1x.*t;
  y = v1y.*t + 0.5.*ay.*t.^2;
  
  
% GRAPHICS ============================================================

 figure(1)
   pos = [0.02 0.05 0.3 0.60];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')

   
 for cc = 1 : nT 
     
 % Mary -----------------------------------------------------    
 subplot('position',[0.08,0.08,0.2,0.50])
     xP = 0; yP = y(cc);
     Hplot = plot(xP,yP,'ok');
     set(Hplot,'markerfacecolor','k','markersize',12)
     hold on
     xlim([-1 1])
     ylim([-1 21])  
     
     set(gca,'fontsize',14)
     H = title('MARY'); 
     set(H,'fontweight','normal','color','r','fontsize',14)
     
     hold off
     set(gca,'xtick',[])
     set(gca,'ytick',[])
     
 % Steve --------------------------------------------------------------
  subplot('position',[0.36,0.08,0.6,0.50])
     xP = x(cc); yP = y(cc);
     Hplot = plot(xP,yP,'ok');
     set(Hplot,'markerfacecolor','k','markersize',10)
     hold on
    
     set(gca,'fontsize',14)
     H = title('STEVE'); 
     set(H,'fontweight','normal','color','b','fontsize',14)
     txt = '  mpSR05.m';
     H = text(15,-2,txt);
     set(H,'fontweight','normal','color','k','fontsize',8)
     set(gca,'xtick',[])
     set(gca,'ytick',[])
     axis equal
     hold off
     xlim([-1 21])
     ylim([-1 21])  
 
 % Steve and Mary --------------------------------------------------------------
 subplot('position',[0.08,0.65,0.8,0.30])
     xP = x(cc); yP = y(cc);
     Hplot = plot(xP,yP,'ok');
     set(Hplot,'markerfacecolor','k','markersize',10)
     hold on
     xP = x(cc); yP = 0;
     Hplot = plot(xP,yP,'sr');
     set(Hplot,'markerfacecolor','r','markersize',30)
     
     xP = 5; yP = -2;
     Hplot = plot(xP,yP,'sb');
     set(Hplot,'markerfacecolor','b','markersize',30)

     set(gca,'fontsize',14)
    
     set(gca,'ytick',[]) 
     axis equal
    
     xlim([-5 21])
     ylim([-5 21]) 
     axis off
     
     pause(0.1)
     hold off
  
         frame1 = frame1 + 1;
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
         % On the first loop, create the file. In subsequent loops, append.
         if frame1 == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
end

   