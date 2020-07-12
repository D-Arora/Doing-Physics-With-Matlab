% cnsFNA.m

% FITZHUGH-NAGUM0 MODEL.
% PHASE PLANE ANALYSIS: membeane potential v and recovery variable w.
% ode45 function used to solve the pair of couped differential equations.
% Variable y gives v and w.
% Most of the parameters are specified in the INPUT SECTION.
% For different model may be necessay to make chnages to the script.
% Results are displayed in as a phase space plot and a time evoltuion plot.
     
% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180913 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/cnsFN.pdf
  

clear 
close all
clc

tic

% INPUTS SECTION ========================================================

% External current stimulus
Iext = 0.31;
% Initial conditions: y0 = [ v(0), w(0)] 
y0 = [-2.8; -1.8];     
% Simulation time [tMin tMax]
tSpanT = [0 200];
% Phase space plot dimensions X axis Lx(1)to Lx(2) / Y axis Ly(1) to Ly(2)
Lx = [-3 3]; Ly = [-2 3];
% Phase plot ticks
ticksX = Lx(1):1:Lx(2);
ticksY = Ly(1):1:Ly(2); 
% Phase space setup for vector field: number of vectors nX [16] 
nX = 20;
% Fitzhugh-Nagoma model parameters
a = 0.7; b = 0.8; c = 12.5;

    

% CALCULATION SECTION ===================================================
% K use to pass variables to function FNode
   K = [Iext; a; b; c];
   
% Solve differental equation using ode45
  [t,y] = ode45(@(t,y) FNode(t,y,K), tSpanT,y0);

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
Sp = vpasolve(p-p^3/3-(p+0.7)/0.8+0.3 == 0,p,[-3 3]);
Sq = (Sp+a)/b;
Sp = double(Sp); Sq = double(Sq);
disp('Critical point');
fprintf('   v_C =  %2.2f\n', Sp);
disp('   ')
fprintf('   w_C =  %2.2f\n', Sq);


  
% GRAPHICS SECTION=======================================================  

FS = 14;  % fontsize

% phase space plot: v vs w   ------------------------------------------
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

% NULLCLINES
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
   tm2 = num2str(Iext,'%3.2f');
   tm3 = '  v_C = ';
   tm4 = num2str(Sp,'%3.2f');
   tm5 = '  w_C = ';
   tm6 = num2str(Sq,'%3.2f');
   tm = [tm1 tm2 tm3 tm4 tm5 tm6];   
   hT = title(tm,'FontName','Courier');  
   
% time evolution of v and w   
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

  legend('v','w','location','north','orientation','horizontal')
  xlabel('t')
  ylabel('v & w')
  title(tm,'fontName','Courier')
  grid on
  set(gca,'fontsize',FS)
  box on 
  
 
% FUNCTIONS ===========================================================

function dydt = FNode(t,y,K)
   a = K(2); b = K(3); c = K(4); I = K(1);

   dydt = [(y(1) - y(1)^3/3 - y(2) + I); (1/c)*(y(1) + a - b*y(2))];

end


