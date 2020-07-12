% mpSR03.m
% 190618  Matlab R2018b

% SPECIAL RELATIVITY ANIMATIONS
% Einstein's Postulates:
% The speed of light is independent of an observer.
% Water waves propagating as seen by two observers Mary and Steve

% DOING PHYSICS ONLINE 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ONLINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Documentation
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod72/mod72A.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

close all
clear
clc

% ANIMATION and GIF file ==============================================

% file name for animated gif   
    ag_name = 'ag_SR.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.05;
    frame1 = 0;

% SETUP ===============================================================
% Steve: Fixed frame of reference (stationary)
% Mary:  Moving frame of reference (w.r.t Steve's frame)
% Number of time steps
  N = 121;
% Velocity of Mary's frame of reference w.r.t. Steve's frame of reference
  v = 8;
% Speed of light
  c = 13;
% Time in the Steve's frame and time step
  t = linspace(0,18,N);
  dt = t(2) - t(1);
% Position of Mary w.r.t. Steve
  xM = -65 + v.*t;
% limits for X axis
  xMin = -60;
  xMax = 60;
  xLmin = -85;
  xLmax = 85;
% Circle radius
  R = 0;
  
  
% GRAPHICS ============================================================

 figure(1)
   pos = [0.02 0.05 0.3 0.60];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
    
for cc = 1 : N 
  subplot('position',[0.08,0.75,0.8,0.20])
     xlim([xMin xMax])
     ylim([-1 1])
 % Steve's position    
      xP = 0; yP = -0.4;
      Hplot = plot(xP,yP,'sb');
      set(Hplot,'markerfacecolor','b','markersize',30)
      hold on
     
 % Mary's position  
     xP = xM(cc);  yP = 0.4;    
     Hplot = plot(xP,yP,'sr');
     set(Hplot,'markerfacecolor','r','markersize',30)

% Light flash ON
     if xM(cc) > -10
      xP = 0; yP = 0;
      Hplot = plot(xP,yP,'og');
      set(Hplot,'markerfacecolor','g','markersize',20) 
     end

% Light flash OFF
     if xM(cc) > 10 
      xP = 0; yP = 0;
      Hplot = plot(xP,yP,'ow');
      set(Hplot,'markerfacecolor','w','markersize',20) 
     end
     
     H = title('As Mary passes Steve, a stone is dropped into the water'); 
     set(H,'fontweight','normal','fontsize',14)
     axis off
     xlim([xMin xMax])
     ylim([-1 1]) 
     pause(0.1)
     hold off 
 
 % Steve's Frame of reference  ----------------------------------------  
  subplot('position',[0.27,0.44,0.4,0.3])
     xlim([xLmin xLmax])
     ylim([xLmin xLmax])  
     axis square
     box on
 
      xP = 0; yP = 0;
      Hplot = plot(xP,yP,'sb');
      set(Hplot,'markerfacecolor','b','markersize',15)
      axis square
      hold on
      
     if xM(cc) > -10  && R < 125 
         [xP, yP] = circle(R,0,0);
         plot(xP,yP,'g','linewidth',4) 
                         
         xlim([xLmin xLmax])
         ylim([xLmin xLmax]) 
         hold off    
     end 
     
     H = title('Steve: Frame of reference'); 
     set(H,'fontweight','normal','color','b','fontsize',14)
     axis off
     
 % Mary's Frame of reference   ----------------------------------------   
  subplot('position',[0.27,0.08,0.4,0.3])
     xlim([xLmin xLmax])
     ylim([xLmin xLmax])  
     axis square
     box on
 
     xP = 0; yP = 0;
      Hplot = plot(xP,yP,'sr');
      set(Hplot,'markerfacecolor','r','markersize',15)
      axis square
      hold on
     
     if xM(cc) > -10  && R < 125
         [xP, yP] = circle(R,0,0);
         xP = xP - xM(cc); 
         plot(xP,yP,'g','linewidth',4)  
                 
         xlim([xLmin xLmax])
         ylim([xLmin xLmax])  
         hold off
         axis square
         R = R + c*dt;
     end
     
     H = title('Mary: Frame of reference');
     set(H,'fontweight','normal','color','r','fontsize',14)
     
     H = text(50,-90,'mpSR04.m');
     set(H,'fontweight','normal','color','k','fontsize',8)
     
    axis off
    
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
  
  function [xP, yP] = circle(aC,xC,yC)
    theta = linspace(0,2*pi,501);
    xP = xC + aC.*cos(theta);
    yP = yC + aC.*sin(theta);
  end
   