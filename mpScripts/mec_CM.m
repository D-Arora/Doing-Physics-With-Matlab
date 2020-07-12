% mec_CM.m

% Animation of an object executing uniform circular motion.
%   The centripetal acceleration vector and the tangential velocity
%   are shown as the object moves around a circle.  

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% VISUAL PHYSICS ONLINE
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod5new/mod52B.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Script calls the function DrawArrow.m (download)
%   DrawArrow(zT,magV,angleV,L,W,LW,col)
%      zT (x,y) coordinates of vector tail
%      magV     length of vector (magnitude)
%      angleV   angle of vector w.r.t. X axis
%      L        length of arrow head
%      W        width of arrow head
%      LW       width of arrow for plotting
%      col      RGB color of arrow  [x x x]  

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190905   Matlab 2018b


clear
close all
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_UCM.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.0;
    frame1 = 0;

    
% SETUP ===============================================================
% Radius of circle
  R = 10;
% Period of motion
  T = 1;
% Number of time steps;
  nT = 400;
% angular frequency
  w = 2*pi/T;
% time interval
  t = linspace(0,T,nT);
% (x,y) coordinates for points on circle
  x = R.*cos(w*t);
  y = R.*sin(w*t);
  
  
% GRAHICS =============================================================
  figure(1)
   pos = [0.1 0.2 0.25 0.35];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   box on
  
for ct = 1 : nT
  
% Plot circle
  xP = x; yP = y;
  plot(xP,yP,'linewidth',2)
  axis equal;
  hold on
 
% Plot Origon
  Hplot = plot(0,0,'ko');
  set(Hplot,'markerfacecolor','k','markersize',8)
  
% Plot centripetal acceleration (color orange)
  zT = x(ct) + 1i*y(ct);
  angleV = atan2(y(ct),x(ct)) + pi;
  magV = 8; L = 1; W = 0.5; LW = 5; col = [1 0.6 0];
  DrawArrow(zT,magV,angleV,L,W,LW,col)

% Plot tangential velocity  (color green)
  xP = x(ct); yP = y(ct);
  plot(xP,yP,'ko')
  zT = x(ct) + 1i*y(ct);
  angleV = atan2(y(ct),x(ct)) + pi/2;
  magV = 5; L = 1; W = 0.5; LW = 4; col = [0 1 0];
  DrawArrow(zT,magV,angleV,L,W,LW,col)
  
% Plot object 
  xP = x(ct); yP = y(ct);
  Hplot = plot(xP,yP,'ko');
  set(Hplot,'markerfacecolor','k','markersize',12);
  
  text(10,-12,'mecCM.m','fontsize',8) 
  axis([-12 12 -12 12])
  axis off
  set(gca,'position',[0.05 0.05 0.9 0.9])
  pause(0.001)
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
  