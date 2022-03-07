% cnsHindmarshA.m

% Hindmarsh Model - Bursting Neurons - 2 coupled equations.
% PHASE PLANE ANALYSIS: membrane potential v and recovery variable w.
% ode45 function used to solve the pair of coupled differential equations.
% Variable y gives v and w.
% Most of the parameters are specified in the INPUT SECTION.
% For different models, it may be necessary to make changes to the script.
% Results are displayed in as a phase space plot and a time evolution plot.
     
% Ian Cooper
% email: matlabvisualphysics@gmail.com

% 180925 / Matlab version R2018b

% DOING PHYSICS WITH MATLAB 
%  https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    
  

clear 
close all
clc

tic

% INPUT SECTION ========================================================

% External current stimulus  [Iext = 1.0]
Iext = 1;

% Initial conditions: y0 = [ v(0) = -2.8, w(0) = -1.8] 
y0 = [1; 0];

% Simulation time [tMin = 0 tMax = 200]
tSpan = [0 200];

% Phase space plot dimensions: X axis [-3 3] / Y axis [-2 3]
Lx = [-2 2]; Ly = [-20 2];

% Phase plot ticks
ticksX = Lx(1):1:Lx(2);
ticksY = Ly(1):5:Ly(2); 

% Phase space setup for vector field: number of vectors nX [20] 
  nX = 20;

% Model parameters
  K(1) = 1; K(2) = -1; K(3) = 3; % K use to pass variables to function FNode  
  K(4) = 1; K(5) = -5; K(6) = -1;
  K(7) = Iext;

 

% CALCULATION SECTION ===================================================


% Solve differential equations using ode45
  [t,y] = ode45(@(t,y) FNode(t,y,K), tSpan,y0);

% Vector field 
  x1 = linspace(Lx(1),Lx(2),nX);
  x2 = linspace(Ly(1),Ly(2),nX);
  [xx, yy] = meshgrid(x1,x2);
  f = ( K(1).*yy + K(2).*xx.^3 +K(3).*xx.^2 + Iext );   
  g = (K(4) + K(5).*xx.^2 + K(6).*yy); 
  fs = f./sqrt(f.^2 + g.^2);    % unit vectors
  gs = g./sqrt(f.^2 + g.^2);
  
 % External current stimulus   To chnage IS must alos change Iext i94n
 %                             function 
  tLEN = length(t);
  IS = zeros(tLEN,1);
  for c = 1 : tLEN
      if (t(c) > 20  &&  t(c) < 40); IS(c) = Iext; end
  end   
  
% Critical point - output to Command Window
syms p
Sp = vpasolve( -K(1)*(K(4) + K(5)*p^2)/K(6) + K(2)*p^3 + K(3)*p^2 + Iext == 0,p,[-2 2]);
%Sp = vpasolve(p-p^3/3-(p+a)/b + Iext == 0,p,[-3 3]);
Sq = -( K(4) + K(5).*Sp.^2 )./K(6);
Sp = double(Sp); Sq = double(Sq);
vC = Sp; wC = Sq;
disp('Critical point');
fprintf('   v_C =  %2.2f\n', vC);
disp('   ')
fprintf('   w_C =  %2.2f\n', wC);


  
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
     vV = -( K(2).*v.^3 + K(3).*v.^2 + Iext )./K(1);
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
    xP = vC; yP = wC;    
    Hplot = plot(xP,yP,'ok');
      set(Hplot,'markersize',6,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0])
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
  pos = [0.05 0.05 0.29 0.50];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

subplot(3,1,1)
   xP = t; yP = IS;  % v
   plot(xP,yP,'b','linewidth',2)
   xlabel('t  [ms]')
   ylabel('I_{ext}')
   title('current stimulus','FontWeight','normal')
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
  
 
  toc  
 
% FUNCTIONS ===========================================================

function dydt = FNode(t,y,K)
  % y(1) == v; y(2) == w;
   Iext = 0;
   if t > 20; Iext = K(7); end
   if t > 40; Iext =   0; end 
   
 %  dydt = [(y(1) - y(1)^3/3 - y(2) + Iext); (1/c)*(y(1) + a - b*y(2))];
    dydt(1) = K(1)*y(2) + K(2)*y(1)^3 + K(3)*y(1)^2 + Iext;
    dydt(2) = K(4) + K(5)*y(1)^2 + K(6)*y(2);
    dydt = dydt';
end


