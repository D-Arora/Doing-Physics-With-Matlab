% oscCMSHM.m

% Animation: Connection between circulation motion and SHM

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
% 190906   Matlab 2018b




clear 
close all
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'agCMSHM.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.1;
    frame1 = 0;

    
% SETUP ===============================================================
% X axis 
  xMin = 15;                
  xMax = 60;                  
  Nx = 800;                     
  x = linspace(xMin,xMax,Nx);   

% Time / Y-axis (wave function y)
% amplitude wave function  / radius of generating circle
   A = 10;
% wavelength / period / angular frequency   
   wL = 15;                     
   T = 20;                     
   w = 2*pi/T;
% initialize wave function
  y = zeros(1,Nx);  
% number of frames for animation / initial time/ time increment
  Nt = 60;
  t = 0;                     
  dt = T/Nt;              


% GRAPHICS ================================================================
figure(1)
% Setup for plot window
 set(gcf,'units','normalized'); 
 set(gcf,'position',[0.1 0.1 0.6 0.30]);
 set(gcf,'Color',[1 1 1]);
 
% Circle
  p = linspace(0,2*pi,500);
  xc = A .* cos(p);
  yc = A .* sin(p);

% Loop for each frame of animation: c is the frame counter
for c = 1 : Nt+1

% Travelling wave
  y = A .* sin(2*pi*t/T - (2*pi/wL).*x);
  plot(x,y,'b','LineWidth',3);
  hold on
 
% Initial point - Y projection of circle
  plot(x(1),y(1),'o','MarkerSize',12,'MarkerFaceColor',[1 0 0])

% Rotating radius vector and circle
  plot([0 A*cos(w*t)],[0 A*sin(w*t)],'lineWidth',3,'Color',[1 0 0]);
  plot(xc,yc,'k','lineWidth',2);

% Y projection    
  plot([0 x(1)], [y(1) y(1)],'k-')
  plot([0 0], [0 y(1)],'k-')
  
  text(55,-12,'oscCMSHM.m','fontsize',8) 
  
  hold off
  axis equal
  axis off
  
  pause(0.1);  
  
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
   
   t = t + dt;               % increment time t
end









