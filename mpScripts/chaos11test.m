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


% Choose coupled D.E. *************************************************
 flagC = 4;
% 1 dx/dt = x*(x^2-1)    dy/dt = y  
% 2 dx/dt = (x+1)*y      dy/dt = y
% 3 dx/dt = x - y        dy/dt = 2x - y - x^2
% 4 dx/dt = (2+x)*(-x+y) dy/dt = (2-x)*(x+y)
% 5 dx/dt = x - 2y -1    dy/dt = 2x - 3y - 3     linear

% Simulation time tMax  /  Boundary for plots LxL
% Number of calculations
  switch flagC
      case 1
        N = 50e5;
        tMax = 5;
        L = 3;
        tmC = 'dx/dt = x(x^2 - 1)   dy/dt = y ';
      case 2
        N = 10e6;
        tMax = 3;
        L = 5;
        tmC = 'dx/dt = (x + 1)y     dy/dt = y';
      case 3
        N = 50e5;   
        tMax = 30;
        L = 3; 
        tmC = 'dx/dt = x - y       dy/dt = 2x - y - x^2';
      case 4
        N = 50e6;
        tMax = 1;
        L = 3;
        tmC = 'dx/dt = (2+x)(-x+y) dy/dt = (2-x)(x+y)';
      case 5
        N = 50e5;
        tMax = 15;
        L = 10;
        tmC = 'dx/dt = x - 2y - 1    dy/dt = 2x - 3y -3';  
 end
% *********************************************************************

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
  switch flagC
      case 1
       f = (xx).*(xx.^2-1);   
       g = yy; 
      case 2
       f = (xx+1).*yy;    
       g = xx.*(yy+3); 
      case 3
       f = xx - yy;
       g = 2.*xx - yy - xx.^2;
      case 4
       f = (2 + xx).*(-xx + yy);
       g = (2 - xx).*( xx + yy);
      case 5
       f = xx - 2.*yy -1;
       g = 2.*xx - 3.*yy -3;
  end 
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
  switch flagC
      case 1
        xP = [0 0]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1.2)
        xP = [-1 -1]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1.2)
        xP = [1 1]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1.2)
        xP = [-L L]; yP = [0 0];
          plot(xP,yP,'m','linewidth',1)
        
      case 2
        xP = [-L L]; yP = [0 0];
          plot(xP,yP,'r','linewidth',1)
        xP = [-1 -1]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1)  
        xP = [0 0]; yP = [-L L];
          plot(xP,yP,'m','linewidth',1)
        xP = [-L L]; yP = [-3 -3];
          plot(xP,yP,'m','linewidth',1)
      case 3
        xP = [-L L]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1) 
        xP = linspace(-L,L, 100); yP = 2.*xP - xP.^2;
          plot(xP,yP,'m','linewidth',1)
      case 4
        xP = [-2 -2]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1)
        xP = [-L L]; yP = [-L L];
          plot(xP,yP,'r','linewidth',1)
        xP = [2 2]; yP = [-L L];
          plot(xP,yP,'m','linewidth',1)
        xP = [-L L]; yP = [L -L];
          plot(xP,yP,'m','linewidth',1) 
      case 5
        xP = linspace(-L,L,100); yP = (xP-1)./2;
         plot(xP,yP,'r','linewidth',1)
        yP = (2.*xP - 3)./3;
         plot(xP,yP,'m','linewidth',1)   
 end
  
   
for cc = 1:numT
   
  [x(1), y(1)] = ginput(1);
  
  % 2ndt time step
  % ************************************************
  switch flagC
      case 1
        x(2) = x(1) + dt* x(1)*(x(1)^2-1); 
        y(2) = y(1) + dt* (y(1)+0);   
      case 2
        x(2) = x(1) + dt* y(1)*(x(1)+1); 
        y(2) = y(1) + dt* x(1)*(2*y(1)+3);
      case 3
        x(2) = x(1) + dt* (x(1) - y(1)); 
        y(2) = y(1) + dt* (2*x(1) - y(1) - x(1)^2);
      case 4
        x(2) = x(1) + dt* (2 + x(1))*(-x(1) + y(1)); 
        y(2) = y(1) + dt* (2 - x(1))*( x(1) + y(1));  
      case 5
        x(2) = x(1) + dt* (x(1) - 2*y(1) -1); 
        y(2) = y(1) + dt* (2*x(1) - 3*y(1) -3);   
   end 
  % ************************************************
  
% Solve coupled D.E. by finite difference method 

  flagS = 0; c = 3;
  while flagS == 0
   
% ************************************************* 
   switch flagC
       case 1
         x(c) = x(c-2) + 2*dt * (x(c-1)+0)*(x(c-1)^2-1); 
         y(c) = y(c-2) + 2*dt *  1*(y(c-1)+0);
       case 2
         x(c) = x(c-2) + 2*dt * (x(c-1)+1)*y(c-1); 
         y(c) = y(c-2) + 2*dt *  x(c-1)*(y(c-1)+3);
       case 3
         x(c) = x(c-2) + dt* (x(c-1) - y(c-1)); 
         y(c) = y(c-2) + dt* (2*x(c-1) - y(c-1) - x(c-1)^2);
       case 4
         x(c) = x(c-2) + dt* (2 + x(c-1))*(-x(c-1) + y(c-1)); 
         y(c) = y(c-2) + dt* (2 - x(c-1))*( x(c-1) + y(c-1));
       case 5
        x(c) = x(c-2) + dt* (x(c-1)   - 2*y(c-1) -1); 
        y(c) = y(c-2) + dt* (2*x(c-1) - 3*y(c-1) -3); 
   end
% *************************************************
   
   if abs(x(c)) > 10*L; flagS = 1; end
   if abs(y(c)) > 10*L; flagS = 1; end
   if c > N-10; flagS = 1; end
   c = c+1;
  end
  
   XY = 1:c-1;  % end index for x and y
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

