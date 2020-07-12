% cnsIzhikevichA.m

% SPIKING NEURONS IZHIKEVICH MODEL.
% PHASE PLANE ANALYSIS: membrane potential v and recovery variable w.
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
%    ../mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/cnsFN.pdf
  

clear 
close all
clc

tic

% INPUT SECTION ========================================================

% Izhihevich model parameters
a = 0.02; b = 0.2; c = -65; d = 8;
k(1) = 0.04; k(2) = 5; k(3) = 140; k(4) = -1; k(5) = 1;


% External current stimulus  [Iext = 1.0]
Iext = 10;
% Initial conditions: y0 = [ v(0) = -2.8, w(0) = -1.8] 
v0 = -65; u0 = b*v0;
y0 = [v0; u0];
% Simulation time [tMin = 0 tMax = 200]
tSpan = [0 500];
% Phase space plot dimensions: X axis [-3 3] / Y axis [-2 3]
Lx = [-100 50]; Ly = [-10 10];
% Phase plot ticks
ticksX = Lx(1):10:Lx(2);
ticksY = Ly(1):20:Ly(2); 
% Phase space setup for vector field: number of vectors nX [20] 
nX = 20;

    

% CALCULATION SECTION ===================================================
% K use to pass variables to function FNode
   K = [Iext; a; b; c; d; k(1); k(2); k(3); k(4); k(5)];
   
% Solve differential equations using ode45
  [t,y] = ode45(@(t,y) FNode(t,y,K), tSpan,y0);

% Vector field 
  x1 = linspace(Lx(1),Lx(2),nX);
  x2 = linspace(Ly(1),Ly(2),nX);
  [xx, yy] = meshgrid(x1,x2);
  f = k(1)*xx.^2 + k(2)*xx + k(3) + k(4)*yy + k(5)*Iext;   
  g = a.*(b.*xx - yy); 
  fs = f./sqrt(f.^2 + g.^2);    % unit vectors
  gs = g./sqrt(f.^2  +g.^2);
  
% Nullclines
  vX = linspace(Lx(1),Lx(2),200);
  vY = -(k(1).*vX.^2 + k(2).*vX + k(3) * k(5).*Iext)/k(4);
  uY = b.*vX;
  
  % Critical point - output to Command Window
% syms p
% Sp = vpasolve(p-p^3/3-(p+a)/b + Iext == 0,p,[-3 3]);
% Sq = (Sp+a)/b;
% Sp = double(Sp); Sq = double(Sq);
% disp('Critical point');
% fprintf('   v_C =  %2.2f\n', Sp);
% disp('   ')
% fprintf('   w_C =  %2.2f\n', Sq);


  
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
    xP = vX;
    yP = vY;
    plot(xP,yP,'r','linewidth',1.5)
% u nullcline 
    yP = uY;
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
   
%    tm1 = 'I_{ext} = ';
%    tm2 = num2str(Iext,'%3.3f');
%    tm3 = '  v_C = ';
%    tm4 = num2str(Sp,'%3.2f');
%    tm5 = '  w_C = ';
%    tm6 = num2str(Sq,'%3.2f');
%    tm = [tm1 tm2 tm3 tm4 tm5 tm6];   
%    hT = title(tm,'FontName','Courier');  

   
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
%   title(tm,'fontName','Courier')
  grid on
  set(gca,'fontsize',FS)
  box on 

  disp('  ')
  toc  
 
% FUNCTIONS ===========================================================

function dydt = FNode(t,y,K)
   Iext = K(1);
   a = K(2); b = K(3); c = K(4); d = K(5);
   k(1) = K(6); k(2) = K(7); k(3) = K(8); k(4) = K(9); k(5) = K(10);
   
   dydt = [k(1)*y(1)^2+k(2)*y(1)+k(3)+k(4)*y(2)+k(5)*Iext; a*(b*y(1)-y(2))];
   
   if dydt(1) > 30
       dydt(1) = 30;
       dydt(2) = dydt(2) + d;
   end
end


