% chaos11.m

% Dynamics of Linear and Nonlinear Systems.
% PHASE PLANE ANALYSIS: Linear and Non-linear Systems.
% A finite difference method is used to solve the paired coupled D.E.
% The number of calculations needs is enormous for stability N ~ 60e6.
% The simulation time not too large else solution diverges  tMax.
% Need to specify the dimensions of the phase space for each D.E. L. 
% The initial values for the state variables x and y are specified using
     % the ginput: click to select I.C. in the phase space plot 
     % the number of trajectories plotted is given by numT.
     % numT = 1 for only 1 trajectory and the time evolution of x & y
% The D.E. needs to be entered three times: indicated by **** .... ***
% The equations for the nullclines are given for each pair of D.E.
% Sample D.E. is chosen using flagC. 
% A time evolution plot of x and y is displayed for the last I.C.

     
% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180913 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/chaos10.pdf
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/chaos11.pdf
%  

clear 
close all
clc

tic

% INPUTS ==============================================================

% Number of trajectories plotted in phase space portrait
numT = 1;

% Model constants
I = 0;
a = 0.7; b = 0.8; c = 3;



% Fitzhugh-Naguno Model

        N = 100e6;
        tMax = 10;
        L = 3;
        tmC = 'dv/dt = c(v - v^3/3 + r + I)   dr/dt = (-1/c)(v - a + b r) ';
      

% CALCULATIONS ========================================================
 
 ticks = -L:1:L;
 
% State variables
  x = zeros(N,1);
  y = zeros(N,1);
     
% Time domain
  t = linspace(0,tMax,N);
  dt = t(2)-t(1);
  
% Phase space setup: number of vectors nX [16] 
   nX = 20;
% Vector field 
  x1 = linspace(-L,L,nX);
  x2 = x1;
  [xx, yy] = meshgrid(x1,x2);
  
 % ***********************************************
  
       f = c.*(xx - xx.^3/3 + yy + I);   
       g = (-1/c).*(xx - a + b.*yy); 
  
 % **********************************************
  
   fs = f./sqrt(f.^2 + g.^2);    % unit vectors
   gs = g./sqrt(f.^2  +g.^2);
  
  
% GRAPHICS   ==========================================================  

FS = 14;  % fontsize

figure(1)
  pos = [0.35 0.05 0.29 0.39];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  title(tmC,'fontweight','normal')
  hold on
  box on
  
% VECTOR FIELD
   hq = quiver(xx,yy,fs,gs);
   set(hq,'color',[0.2 0.2 0.2],'AutoScaleFactor',0.6);
   set(gca,'fontsize',FS)
   xlim([-L L])
   ylim([-L L])
   set(gca,'xtick',ticks);
   set(gca,'ytick',ticks);
   grid on   
 %  hold on
   xlabel('x'); ylabel('y');

% NULLCLINES
% v nullcline
     v = linspace(-L,L,200);
     r = -(v - v.^3/3 + I);
        xP = v; yP = r;
          plot(xP,yP,'r','linewidth',1.2)
          
% r nullcline 
    r = (a - v)/b;
        xP = v; yP = r;
          plot(xP,yP,'m','linewidth',1.2)
  
        
for cc = 1:numT
   
  [x(1), y(1)] = ginput(1);
  
  % 2nd time step
  % ************************************************
  
        x(2) = x(1) + dt * ( c*(x(1) - x(1)^3/3 + y(1) + I) ); 
        y(2) = y(1) + dt * ( (-1/c)*(x(1) - a + b*y(1)) );   
     
  % ************************************************
  
% Solve coupled D.E. by finite difference method 

  flagS = 0; k = 3;
  while flagS == 0
   
% ************************************************* 
         x(k) = x(k-2) + 2*dt* ( c*(x(k-1) - x(k-1)^3/3 + y(k-1) + I) ); 
         y(k) = y(k-2) + 2*dt* ( (-1/c)*((x(k)+x(k-1))/2 - a + b*y(k-1)) );   
         
       
  
% *************************************************
   
   if abs(x(k)) > 10*L; flagS = 1; end
   if abs(y(k)) > 10*L; flagS = 1; end
   if k > N-10; flagS = 1; end
   k = k+1;
  end
  
   XY = 1:k-1;  % end index for x and y
   xP = x(1); yP = y(1);   
   Hplot = plot(xP,yP,'o');
   set(Hplot,'markersize',6,'markerFaceColor','b','markerEdgeColor','b')
   stepXY = 1e3;
   xP = x(1:stepXY:length(XY));
   yP = y(1:stepXY:length(XY));
   plot(xP,yP,'b','linewidth',2.5)
   
end    
  
   

   
   
figure(2)
  pos = [0.05 0.05 0.29 0.29];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  plot(t(XY),x(XY),'b','linewidth',2)
  hold on
  plot(t(XY),y(XY),'r','linewidth',2)

  tm1 = 'x_{end} = ';
  tm2 = num2str(x(XY(end)),'%3.2f');
  tLx = [tm1 tm2];
  tm1 = 'y_{end} = ';
  tm2 = num2str(y(XY(end)),'%3.2f');
  tLy = [tm1 tm2];
  legend(tLx,tLy,'location','northoutside','orientation','horizontal')
  xlabel('t')
  ylabel('x & y')
  grid on
  set(gca,'fontsize',FS)
  box on 
 
     




toc

