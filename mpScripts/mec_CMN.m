% mec_CM.m

% Animation of an object executing circular motion with an applied
%   constant torque. The acceleration vectors and the tangential velocity
%   are shown as the object moves around a circle.  

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% VISUAL PHYSICS ONLINE
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod5new/mod52B.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 1909704   Matlab 2018b

% Script calls the function DrawArrow.m (download)
%   DrawArrow(zT,magV,angleV,L,W,LW,col)
%      zT (x,y) coordinates of vector tail
%      magV     length of vector (magnitude)
%      angleV   angle of vector w.r.t. X axis
%      L        length of arrow head
%      W        width of arrow head
%      LW       width of arrow for plotting
%      col      RGB color of arrow  [x x x]  

clear
close all
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_CM.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.0;
    frame1 = 0;

    
% SETUP ==============================================================
% Radius of circle
  R = 10;
% Period of motion
  T = 1;
% Number of time steps;
  nT = 400;
% angular acceleration
  alpha = 1;  
% angular frequency
  w = 2*pi/T;
% time interval
  tMax = 3*T;
  t = linspace(1,tMax,nT);
% (x,y) coordinates for points on circle
  x = R.*cos(w*t);
  y = R.*sin(w*t);
% max angular speed
  wMax = alpha*tMax;

  
% GRAHICS =============================================================
  figure(1)
   pos = [0.1 0.2 0.25 0.35];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   box on
  
for ct = 1 : nT
    
  xP = x; yP = y;
  plot(xP,yP,'linewidth',2)
 
  axis equal;
  hold on

% Origin (0,0) 
  Hplot = plot(0,0,'ko');
  set(Hplot,'markerfacecolor','k','markersize',8)

% (x,y) coordinates of rotating object  
  wN = alpha*t(ct);
  xN = R*cos(wN*t(ct));
  yN = R*sin(wN*t(ct));
  
% tangential velocity (green)
  xP = xN; yP = yN;
  plot(xP,yP,'ko')
  zT = xN + 1i*yN;
  angleV = atan2(yN,xN) + pi/2;
  magV = 10*(wN/wMax); L = 1; W = 0.5; LW = 5; col = [0 1 0];
  DrawArrow(zT,magV,angleV,L,W,LW,col)
  
% centripetal acceleration (orange)
  zT = xN + 1i*yN;
  angleV = atan2(yN,xN) + pi;
  magV = 9*(wN/wMax)^2; L = 1; W = 0.5; LW = 3; col = [1 0.6 0];
  DrawArrow(zT,magV,angleV,L,W,LW,col)
  ac = magV;
  
% tangential acceleration  (orange)
  xP = xN; yP = yN;
  plot(xP,yP,'ko')
  zT = xN + 1i*yN;
  angleV = atan2(yN,xN) + pi/2;
  magV = 5; L = 1; W = 0.5; LW = 3; col = [1 0.6 0];
  DrawArrow(zT,magV,angleV,L,W,LW,col)
  at = magV;

% linear acceleration  (red)
  acc = sqrt(ac^2 + at^2);
  theta = atan2(ac,at);
  magV = acc;
  angleV = theta + angleV;
  L = 1; W = 0.5; LW = 3; col = [1 0 0];
  DrawArrow(zT,magV,angleV,L,W,LW,col)
  
% (x,y) coordinates of object
  xP = xN; yP = yN;
  Hplot = plot(xP,yP,'ko');
  set(Hplot,'markerfacecolor','k','markersize',12);

  text(10,-13,'mecCMN.m','fontsize',8) 
  axis([-14 14 -14 14])
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
  