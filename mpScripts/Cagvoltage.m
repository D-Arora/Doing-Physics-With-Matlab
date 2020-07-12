% File name: pulse1_ag.m
% Ian Cooper
% School of Physics, University of Sydney NSW 2006 Australia
% email: ian.cooper@sydney.edu.au
% Demonstration program: how to create an animated gif file.
% Animation of a pulse
%     y = A^3 / (A^2 + (vt - x)^2)
% Animated gif is saved it will loop continuously.

clear all
close all
clc

% SETUP ===================================================================
% all physical quantities in S.I. units unless stated otherwise
% Setup for animated gif
     ag_name = 'ag_Cagvoltage.gif';   % file name for animated gif
%  Delay in seconds before displaying the next image  
      delay = 0.20;  
%  Frame counter start
      nt = 1;
% Save flagS = 1; NOt save flagS = 0;
      flagS = 0;
                             
% Voltage
A = 10;                     % constant A
T = 20;                      % period
w = 2*pi/T;
N = 200;                   % number of frames for animation
t = linspace(0,3*T,N);              % time increment
v  = A .* sin(2*pi*t/T);
%vP = A .* cos(2*pi*t/T);

% GRAPHICS ================================================================
figure(1)
%   Setup for plot window
set(gcf,'units','normalized'); 
set(gcf,'position',[0.1 0.1 0.6 0.30]);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',12);

xP = t; yP = v;
plot(xP,yP,'k','LineWidth',1);

% circle
p = linspace(0,2*pi,500);
xc = A .* cos(p);
yc = A .* sin(p);

   
for c = 1 : 46%N
    hold off
   
   axis equal
   axis off
   
   
   plot([0 A*sin(w*t(c))],[0 -A*cos(w*t(c))],'lineWidth',1,'Color',[0 0 1]);
   hold on
   xP = [0 A*sin(w*t(c))]; yP = [-A*cos(w*t(c)) -A*cos(w*t(c))];
   plot(xP,yP,'r','LineWidth',3)
   
   xP = [12+t(c) 12+t(c)]; yP = [0 v(c)];
   plot(xP,yP,'r','LineWidth',3);
   
   xP = [12 12+t(end)]; yP = [0 0];
   plot(xP,yP,'k','LineWidth',1);
      
   plot(xc,yc,'k','lineWidth',1);
   
   xP = 12+t; yP = v;
   plot(xP,yP,'m','LineWidth',1);
   
   xP = 12+t(1:c); yP = v(1:c);
   plot(xP,yP,'b','LineWidth',3);
   
   xP = [0 0]; yP = [-A A];
   plot(xP,yP,'k','LineWidth',1);
   
   xP = [-A A]; yP = [0 0];
   plot(xP,yP,'k','LineWidth',1);
   
   hold off
   set(gca,'FontSize',12);
   
   axis equal
   axis off
   pause(0.1)
 
if flagS == 1 
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
end  % if








