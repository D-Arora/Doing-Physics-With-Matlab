% cnsFNTW.m

% FITZHUGH-NAGUM0 TRAVELLING WAVES.
% membrane potential v and recovery variable w.
% ode45 function used to solve the pair of coupled differential equations.
% Variable y gives v and w.
% Most of the parameters are specified in the INPUT SECTION.
% For different models, it may be necessary to make changes to the script.
% Results are displayed in as a phase space plot and a time evolution plot.
     
% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180925 / Matlab version R2018b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/cnsFN.pdf
  

clear 
close all
clc

tic

% INPUT SECTION ========================================================

% External current stimulus  [Iext = 1.0]
Iext = 1.5;
% Initial conditions: y0 = [ v(0) = -2.8, w(0) = -1.8] 
y0 = [-1.5; -3/8];
% Simulation time [tMin = 0 tMax = 200]
tSpan = [0 0.5];
% Dimensions X axis [0 30]
Lx = [-3 3]; 
Ly = Lx;
% Phase plot ticks
ticksX = Lx(1):1:Lx(2);
ticksY = Ly(1):1:Ly(2); 
% Phase space setup for vector field: number of vectors nX [20] 
nX = 20;
% Fitzhugh-Nagoma model parameters
K(1) = 10; K(2) = 1;
K(3) = 0.8; K(4) = 1.5; K(5) = 1.25; K(6) = 1;





    

% CALCULATION SECTION ===================================================
% K use to pass variables to function FNode
 %  K = [Iext; a; b; c];
   
% Solve differential equations using ode45
  [t,y] = ode45(@(t,y) FNode(t,y,K), tSpan,y0);

% Vector field 
  x1 = linspace(Lx(1),Lx(2),nX);
  x2 = linspace(Ly(1),Ly(2),nX);
  [xx, yy] = meshgrid(x1,x2);
  f = (xx - xx.^3/3 - yy + Iext);   
  g = (1/c).*(xx + a - b.*yy); 
  fs = f./sqrt(f.^2 + g.^2);    % unit vectors
  gs = g./sqrt(f.^2  +g.^2);
  
% Critical point - output to Command Window
syms p
Sp = vpasolve(p-p^3/3-(p+a)/b + Iext == 0,p,[-3 3]);
Sq = (Sp+a)/b;
Sp = double(Sp); Sq = double(Sq);
disp('Critical point');
fprintf('   v_C =  %2.2f\n', Sp);
disp('   ')
fprintf('   w_C =  %2.2f\n', Sq);


  
% GRAPHICS SECTION=======================================================  

   FS = 14;  % fontsize

% Phase space plot: v vs w   ------------------------------------------
figure(1)
   pos = [0.35 0.05 0.29 0.39];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
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
   xlabel('membrane potential v'); ylabel('recovery variable  w');

% v nullcline
     v = linspace(Lx(1),Lx(2),200);
     w = (v - v.^3/3 + Iext);
        xP = v; yP = w;
          plot(xP,yP,'r','linewidth',1.5)

          % r nullcline 
    w = (v + a)/b;
        xP = v; yP = w;
          plot(xP,yP,'m','linewidth',1.5)
          
% Phase portrait trajectory  
    xP = y(:,1); yP = y(:,2);
      plot(xP,yP,'b','linewidth',2)
    xP = y(1,1); yP = y(1,2);    % initial conditions: start of trajectory
      Hplot = plot(xP,yP,'o');
      set(Hplot,'markersize',8,'markerfacecolor',[0 1 0],'markeredgecolor',[0 1 0])
    xP = y(end,1); yP = y(end,2);   % end of trajectory
      Hplot = plot(xP,yP,'o');
      set(Hplot,'markersize',8,'markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
   
   tm1 = 'I_{ext} = ';
   tm2 = num2str(Iext,'%3.3f');
   tm3 = '  v_C = ';
   tm4 = num2str(Sp,'%3.2f');
   tm5 = '  w_C = ';
   tm6 = num2str(Sq,'%3.2f');
   tm = [tm1 tm2 tm3 tm4 tm5 tm6];   
   hT = title(tm,'FontName','Courier');  

   
% Time evolution of v and w   
figure(2)
  pos = [0.05 0.05 0.29 0.29];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = t; yP = y(:,1);  % v
    plot(xP,yP,'b','linewidth',2)
  hold on
    yP = y(:,2);        % w
  plot(xP,yP,'r','linewidth',2)

  legend('v','w','location','south','orientation','horizontal')
  xlabel('t')
  ylabel('v & w')
  title(tm,'fontName','Courier')
  grid on
  set(gca,'fontsize',FS)
  box on 

  disp('  ')
  toc  
 
% FUNCTIONS ===========================================================

function dydt = FNode(t,y,K)
   a = K(2); b = K(3); c = K(4); Iext = K(1);
   %Iext = 0.003 * t;
   %if t > 50; Iext = 0.5; end
   %if t > 55; Iext =   0; end 
   
   dydt = [(y(1) - y(1)^3/3 - y(2) + Iext); (1/c)*(y(1) + a - b*y(2))];

end


