% mpSR06.m
% 190622  Matlab R2018b

% SPECIAL RELATIVITY ANIMATIONS
% RELATIVE VELOCITIES
% Michelson-Morely Experiment
% Boat moving through water
% Observers Mary (M) on boat and Steve (S) on the bank of the river
% Frames of reference: Steve, Mary and Water

% DOING PHYSICS ONLINE 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ONLINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Documentation
%    http://www.physics.usyd.edu.au/teach_res/hsc/spP/mod7/m72B.htm

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

close all
clear
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_SR.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.15;
    frame1 = 0;

% SETUP ===============================================================
% Steve: Fixed frame of reference (stationary)
% Mary:  Moving frame of reference (w.r.t Steve's frame)
% Water: Moving frame of referecne (w.r.t. Steve's frame)

% Number of time steps
  N = 101;
% Velocity Water w.r.t. Steve
  vWS = 2;
% Velocity Mary w.r.t. Water
  vMW = 10;
% Angle Mary w.r.t river bank [deg]
  theta = 160;
% Time 
  tMax = 10;
  t = linspace(0,tMax,N);
  dt = t(2) - t(1);
% Velocity Mary w.r.t. Steve
  vMS1 = vMW + vWS;                % Mary moving downstream
  vMS2 = -(vMW - vWS);             % Mary moving upstream
  vMSx = vMW * cosd(theta) + vWS;  % Mary across river
  vMSy = vMW * sind(theta);
  vMS3 = sqrt(vMSx^2 + vMSy^2);
  
% Position of Mary in boat 
  xM1 = 0    + vMS1 .*t;
  yM1 = zeros(N,1);
  xM2 = 100 + vMS2.*t;
  yM2 = zeros(N,1);
  
  xM3 = 65 + vMSx.*t;
  yM3 = -5 + vMSy.*t;
  indY = find(yM3 > 5,1);
  xM3(indY:N) = xM3(indY);
  yM3(yM3 > 5) = 5;
  
  
% GRAPHICS ============================================================
% Position of Steve / Mary / light pulses
 figure(1)
   pos = [0.35 0.2 0.4 0.7];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
      
   for cc = 1 : N
   
   subplot(3,1,1)
     hold off
     plot([0 100],[5 5],'g','linewidth',3)  
     hold on
     
     plot([0 100],[-5 -5],'g','linewidth',3)  
     xP = xM1(cc); yP = yM1(cc);
     Hplot = plot(xP,yP,'sr');
     set(Hplot,'markerfacecolor','r','markersize',30)
     
     xlim([0 100])
     ylim([-5 5])
     
     txt = ['t = ' num2str(t(cc),'%2.1f') ' s' ];
     Htext = title(txt);
     set(Htext,'fontsize',14)
     
     txt = ['v_{MS} = ' num2str(vMS1,'%2.1f') ' m/s   v_{WS} = ' ...
         num2str(vWS,'%2.1f') ' m/s   v_{MW} = ' num2str(vMW,'%2.1f') ' m/s' ];
     Htext = text(10,10,txt);
     set(Htext,'fontsize',14)
     
     set(gca,'fontsize',14)
     axis equal
     
  subplot(3,1,2)
     hold off
     plot([0 100],[5 5],'g','linewidth',3)  
     hold on
     plot([0 100],[-5 -5],'g','linewidth',3)
     
     xP = xM2(cc); yP = yM2(cc);
     Hplot = plot(xP,yP,'sr');
     set(Hplot,'markerfacecolor','r','markersize',30)
     
     xlim([0 100])
     ylim([-5 5]) 
     
     txt = ['v_{MS} = ' num2str(vMS2,'%2.1f') ' m/s   v_{WS} = ' num2str(vWS,'%2.1f') ' m/s   v_{MW} = ' num2str(vMW,'%2.1f') ' m/s' ];
     Htext = text(10,10,txt);
     set(Htext,'fontsize',14)
     
     set(gca,'fontsize',14)
     axis equal
     
  subplot(3,1,3)
     hold off
     plot([0 100],[5 5],'g','linewidth',3)  
     hold on
     plot([0 100],[-5 -5],'g','linewidth',3)
     
     xP = xM3(cc); yP = yM3(cc);
     Hplot = plot(xP,yP,'sr');
     set(Hplot,'markerfacecolor','r','markersize',30)
     
     xlim([30 70])
     ylim([-5 5])  
     
     txt = ['v_{MS} = ' num2str(vMS3,'%2.1f') ' m/s   v_{WS} = ' num2str(vWS,'%2.1f') ' m/s   v_{MW} = ' num2str(vMW,'%2.1f') ' m/s' ];
     Htext = text(33,0,txt);
     set(Htext,'fontsize',14)
     
     
     set(gca,'fontsize',14)
     xlabel('                                   mpSR06.m','fontsize',8)
     axis equal                      
     
     pause(0.01)
      
% Save animation      
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
      
      
      hold off
   end
   