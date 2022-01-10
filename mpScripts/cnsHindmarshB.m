% cnsHindmarshB.m

% Hindmarsh Model: Bursting Neurons - 2 and 3 coupled ODEs 
%    Solved using ode45 
%    Solution Variables: y(:,1)  membrane potential v
%                        y(:,2)  recovery variable w
%                        y(:,3)  adaptation current z

% Most of the parameters are specified in the INPUT SECTION
% Simulation time: tSpan
% Initial conditions:  y0 = [v(0) w(0) z(0)]
% Current pulse:  K(7)    pulse height  Imax    
%                 K(11)   pulse on  / K(12)  pulse off
%                 K(12) - K(11)  pulse duration
% 2 coupled ODEs constants: K(1) to K(6)
% 3 coupled ODEs constants: K(9)  K(10)
%                           K(8) adaptation variable
%                           K(8) = 0 --> 2 coupled ODEs model 
%                  
% For different models, it may be necessary to make changes to the Script

% OUTPUTS: time evolution plots
%          PHASE PLANE ANALYSIS: phase space plot (v-w trajectory)
%                               v - w vector field
%                               v - w vector field
%                               v and w nullclines
%          Equilibrium (critical) points  vC and wC

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% DOING PHYSICS WITH MATLAB 
%     https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%     https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cnsHindmarshA.m.pdf
% Date: 211216 / Matlab version: R2021b
  
close all
clear
clc
tic

% INPUT SECTION: [default values] / time [ms] / parameters [ a.u.] >>>>>

% External current stimulus  [Iext = 1.0]
  Imax = 1;

% Initial conditions: y0 = [ v(0) = -1.51 w(0) = -10.35  z(0) = 0] 
  y0 = [-1.6180; -12.0902; 0];

% Simulation time [tMin = 0 tMax = 200]
%  tSpan = [0 1500];
  tSpan = linspace(9,250, 5001);

% Phase space plot dimensions: X axis [-2 2] / Y axis [-20 2]
  Lx = [-2 2]; Ly = [-20 2];

% Phase plot ticks
  ticksX = Lx(1):1:Lx(2);
  ticksY = Ly(1):5:Ly(2); 

% Phase space setup for vector field: number of vectors nX [20] 
  nX = 20;

% Model parameters  2 coupled ODEs K(8) = 0;  K(10) = v_C
% K use to pass variables to function FNode
  K(1) = 1; K(2) = -1; K(3) = 3;   
  K(4) = 1; K(5) = -5; K(6) = -1;

  K(7) = Imax;

  K(8) = 0e-4;    % adaptative variable  [1e-3]
                  % K(8) = 0 --> 2 coupled ODE model

  K(9) = 2; K(10) = y0(1);
 
% Current pulse on / off  [K(11) = 25  K(12) = 50]
  K(11) = 50; K(12) = 75;


% CALCULATION SECTION ===================================================

% Solve differential equations using ode45
%  opts = odeset('RelTol',1e-10);
  [t,y] = ode45(@(t,y) FNode(t,y,K), tSpan,y0);

% Vector field 
  x1 = linspace(Lx(1),Lx(2),nX);
  x2 = linspace(Ly(1),Ly(2),nX);
  [xx, yy] = meshgrid(x1,x2);
  f = ( K(1).*yy + K(2).*xx.^3 +K(3).*xx.^2 + Imax );   
  g = (K(4) + K(5).*xx.^2 + K(6).*yy); 
  fs = f./sqrt(f.^2 + g.^2);    % unit vectors
  gs = g./sqrt(f.^2 + g.^2);
  
 % External current stimulus   To chnage IS must alos change Iext i94n
 %                             function 
  tLEN = length(t);
  IS = zeros(tLEN,1);
 
  for c = 1 : tLEN
      if (t(c) > K(11)  &&  t(c) < K(12)); IS(c) = Imax ; end  % pulse
  end   
  
% Critical point - output to Command Window
syms p
Sp = vpasolve( -K(1)*(K(4) + K(5)*p^2)/K(6) + K(2)*p^3 + K(3)*p^2 + Imax == 0,p,Lx);
%Sp = vpasolve( -K(1)*(K(4) + K(5)*p^2)/K(6) + K(2)*p^3 + K(3)*p^2 + 0 == 0,p,Lx);

Sq = -( K(4) + K(5).*Sp.^2 )./K(6);
Sp = double(Sp); Sq = double(Sq);
vC = Sp; wC = Sq;
disp('Critical point');
fprintf('   v_C =  %2.2f\n', vC);
disp('   ')
fprintf('   w_C =  %2.2f\n', wC);

% equilibrium points: no adaptation, Imax = 0;
ep(1,1) = -1.6180; ep(1,2) = -12.0902;
ep(2,1) = -1.0000; ep(2,2) = -4.0000;
ep(3,1) = 0.6180; ep(3,2) = -0.9098;



  
% GRAPHICS SECTION=======================================================


figure(1)   % Phase space plot: v vs w   ---------------------------
   pos = [0.35 0.05 0.27 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 14;  % fontsize
   hold on
   box on
  
% VECTOR FIELD
   hq = quiver(xx,yy,fs,gs);
   set(hq,'color',[0.2 0.2 0.2],'AutoScaleFactor',0.3);
   set(gca,'fontsize',FS)
   xlim(Lx)
   ylim(Ly)
   set(gca,'xtick',ticksX);
   set(gca,'ytick',ticksY);
   grid on   
   xlabel('membrane potential v'); ylabel('recovery variable  w');

% v nullcline
  v = linspace(Lx(1),Lx(2),200);
  vV = -( K(2).*v.^3 + K(3).*v.^2 + Imax )./K(1);
  xP = v; yP = vV;
  plot(xP,yP,'r','linewidth',1.5)

% w nullcline 
  wW = -( K(4) + K(5).*v.^2)/K(6) ;
  xP = v; yP = wW;
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

% Critical points
  %xP = vC; yP = wC; 
  xP = ep(:,1); yP = ep(:,2);
  Hplot = plot(xP,yP,'ok');
  set(Hplot,'markersize',6,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0])
   
  
figure(2)    % Time evolution of v and w ------------------------------
  pos = [0.05 0.05 0.27 0.50];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

subplot(3,1,1)
  xP = t; yP = IS;  % I
  plot(xP,yP,'b','linewidth',2)

  hold on
  yP = y(:,3);  % v
  plot(xP,yP,'r','linewidth',2)

  xlabel('t  [ms]')
  ylabel('I_{ext} &  z')
  title('current stimulus & adaptation','FontWeight','normal')
  grid on
  set(gca,'fontsize',FS)

subplot(3,1,2)  
  xP = t; yP = y(:,1);  % v
  plot(xP,yP,'b','linewidth',2)
  xlabel('t  [ms]')
  ylabel('v')
  title('membrane potential','FontWeight','normal')
  grid on
  set(gca,'fontsize',FS)

subplot(3,1,3)  
  xP = t; yP = y(:,2);  % v
  plot(xP,yP,'b','linewidth',2)
  xlabel('t  [ms]')
  ylabel('w')
  title('recovery variable','FontWeight','normal')
  grid on
  set(gca,'fontsize',FS)



figure(3)   % Phase space plot: v vs w   ---------------------------
   pos = [0.35 0.45 0.27 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 14;  % fontsize
   hold on
   box on  


   xP = y(:,1); yP = y(:,2); zP = y(:,3);
   plot3(xP,yP,zP,'b','LineWidth',2)
   
   xP = y(end,1); yP = y(end,2); zP = y(end,3);
   Hplot = plot3(xP,yP,zP,'ro');
   set(Hplot,'markersize',8,'markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
   xP = y(1,1); yP = y(1,2); zP = y(1,3);
   Hplot = plot3(xP,yP,zP,'go');
   set(Hplot,'markersize',8,'markerfacecolor',[0 1 0],'markeredgecolor',[0 1 0])
   hold on

   view(-72,16)
   grid on
   box on

   ax = gca;
   ax.BoxStyle = 'full'; 
   xlabel('v')
   ylabel('w')
   zlabel('z')
  toc  

 
% FUNCTIONS ===========================================================

function dydt = FNode(t,y,K)
  % y(1) == v; y(2) == w;  y(3) == z
   Imax = 0;
   if t > K(11); Imax = K(7); end
   if t > K(12); Imax =   0; end 

   % Added white noise  not shown on plot - needs fixing !!!
   wN = 0.5*randn(1,1); Imax = Imax + wN;

   
    dydt(1) = K(1)*y(2) + K(2)*y(1)^3 + K(3)*y(1)^2 + Imax  - y(3);
    dydt(2) = K(4) + K(5)*y(1)^2 + K(6)*y(2);
    dydt(3) = K(8)*( K(9)*(y(1) - K(10)) - y(3) );
    dydt = dydt';
end


