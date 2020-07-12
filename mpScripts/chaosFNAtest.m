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
%    ../mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/chaos10.pdf
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/chaos11.pdf
%  

clear 
close all
clc

tic

% INPUTS ==============================================================



I = 0.5;
a = 0.7; b = 0.8; c = 12.5;

dT = [0 80];       % Simulation time interval
  y0 = [-2; -2];     % Initial conditions

% Fitzhugh-Naguno Model

        N = 100e6;
        tMax = 10;
        Lx = [-3 3]; Ly = [-2 3];
        tmC = 'dv/dt = c(v - v^3/3 + r + I)   dr/dt = (-1/c)(v - a + b r) ';
      

% CALCULATIONS ========================================================
K = [I; a; b; c];


  
  
  [t,y] = ode45(@(t,y) FNode(t,y,K), dT,y0);
  
%ode45(@(t,y) thisode(t,y,P1), time, It);
 ticksX = Lx(1):1:Lx(2);
 ticksY = Ly(1):1:Ly(2);
 
 
% State variables
%   x = zeros(N,1);
%   y = zeros(N,1);
     

  
% Phase space setup: number of vectors nX [16] 
   nX = 20;
% Vector field 
  x1 = linspace(Lx(1),Lx(2),nX);
  x2 = linspace(Ly(1),Ly(2),nX);
  [xx, yy] = meshgrid(x1,x2);
  
 % ***********************************************
  
       f = (xx - xx.^3/3 - yy + I);   
       g = (1/c).*(xx + a - b.*yy); 
  
  
  
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
   xlim(Lx)
   ylim(Ly)
   set(gca,'xtick',ticksX);
   set(gca,'ytick',ticksY);
   grid on   
 %  hold on
   xlabel('x'); ylabel('y');

% NULLCLINES
% v nullcline
     v = linspace(Lx(1),Lx(2),200);
     r = (v - v.^3/3 + I);
        xP = v; yP = r;
          plot(xP,yP,'r','linewidth',1.5)
          
% r nullcline 
    r = (v + a)/b;
        xP = v; yP = r;
          plot(xP,yP,'m','linewidth',1.5)
  
xP = y(:,1); yP = y(:,2);
plot(xP,yP,'b','linewidth',2)


xP = y(1,1); yP = y(1,2);
    Hplot = plot(xP,yP,'o');
    set(Hplot,'markersize',8,'markerfacecolor',[0 1 0],'markeredgecolor',[0 1 0])
xP = y(end,1); yP = y(end,2);
Hplot = plot(xP,yP,'o');
    set(Hplot,'markersize',8,'markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
   
   
figure(2)
  pos = [0.05 0.05 0.29 0.29];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = t; yP = y(:,1);
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = y(:,2);
  plot(xP,yP,'r','linewidth',2)

%   tm1 = 'x_{end} = ';
%   tm2 = num2str(x(XY(end)),'%3.2f');
%   tLx = [tm1 tm2];
%   tm1 = 'y_{end} = ';
%   tm2 = num2str(y(XY(end)),'%3.2f');
%   tLy = [tm1 tm2];
%   legend(tLx,tLy,'location','northoutside','orientation','horizontal')
  xlabel('t')
  ylabel('x & y')
  grid on
  set(gca,'fontsize',FS)
  box on 
 





function dydt = FNode(t,y,K)
a = K(2); b = K(3); c = K(4); I = K(1);


dydt = [(y(1) - y(1)^3/3 - y(2) + I); (1/c)*(y(1) + a - b*y(2))];


end


