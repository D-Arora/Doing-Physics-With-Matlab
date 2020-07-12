% mpSR01.m
% 190624  Matlab R2018b

% SPECIAL RELATIVITY ANIMATIONS
% SIMULTANEITY
% STEVE: Stationary frame of reference (BLUE)
% MARY: Moving frame of reference (RED)
% EVENT: Two green lights flash simultaneously in Steve's frame
% Ignores length contract 


% DOING PHYSICS ONLINE 
%    ../mphome.htm
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ONLINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Documentation
%    http://www.physics.usyd.edu.au/teach_res/hsc/spP/mod7/m72C.htm

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
    delay = 0.10;
    frame1 = 0;


% SETUP ===============================================================

% Number of time steps
  N = 101;
% Velocity of Mary w.r.t. Steve
  vMS = 6;
% Velocity of light
  c = 12;

% Time 
  tMax = 40;
  t = linspace(0,tMax,N);
  dt = t(2) - t(1);
  
% Initial position of light flashes of Steve and Mary
% Steve S  /  Mary M /  Left L  /  R right
  xSL = -100; xSR =  100;
  xML = -100; xMR =  100;
  
% Mary's position wrt Steve  
  xMS = -100 + vMS.*t;
  
  flagM = 0;

  
% GRAPHICS ============================================================

figure(1)
   pos = [0.35 0.2 0.4 0.2];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
      
for cc = 1 : N
  
% Steve's position 
    xP = 0;  yP = -5;    
    Hplot = plot(xP,yP,'sb');
    set(Hplot,'markerfacecolor','b','markersize',50) 

    hold on
    
% Mary's position   
    xP = xMS(cc);  yP = 5;    
    Hplot = plot(xP,yP,'sr');
    set(Hplot,'markerfacecolor','r','markersize',50)     
    
% Steve: Position of light flashes 
     if xMS(cc) > 0 && xSL < 0  
        xP = xSL;  yP = -5;    
        Hplot = plot(xP,yP,'og');
        set(Hplot,'markerfacecolor','g','markersize',20)
        xSL = xSL + c*dt;  
     end

    if xMS(cc) > 0 && xSR > 0
        xP = xSR;  yP = -5;    
        Hplot = plot(xP,yP,'og');
        set(Hplot,'markerfacecolor','g','markersize',20)
        xSR = xSR - c*dt; 
    end

% Mary: Position of light flashes 
      if xMS(cc) > 0 && xML < xMS(cc)
        if flagM == 0  
          xP = xML;  yP = 5;    
          Hplot = plot(xP,yP,'og');
          set(Hplot,'markerfacecolor','g','markersize',20)
         xML = xML + c*dt; 
        end
        if xML > xMS(cc); flagM = 1; end 
     end 

     if xMS(cc) > 0  && xMR > xMS(cc)
        xP = xMR;  yP = 5;    
        Hplot = plot(xP,yP,'og');
        set(Hplot,'markerfacecolor','g','markersize',20)
        xMR = xMR - c*dt; 
    end
 
      xlim([-100 130])
      ylim([-7 7])
      
      
      H = text(100,-5,'mpSR01.m'); 
      set(H,'fontweight','normal','color','k','fontsize',8)
      
      axis off
      pause(0.1)
      hold off
      
% Save animation 
    if flagAG == 1
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
 end


    
