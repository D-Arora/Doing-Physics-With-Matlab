% em_motor_1.m
% Animation of DC motor - armature turning in magnetic field
% Ian Cooper
% School of Physics University of Sydney
% email: cooper@physics.usyd.edu.au

clear all; close all; clc

N = 200;               % number of frames  N = 200 default

T = 1;                 % paramaters for torque & curret graphs
t = linspace(0,3*T,N);
xt = -10;
dx = 24/N;

figure(1)  % **********************************************************
set(gcf,'color',[1 1 1]);
set(gcf,'units','normalized'); 
set(gcf,'position',[0.1 0.1 0.6 0.7]);
set(gca,'Xlim',[-20 20]);
set(gca,'xcolor',[1 1 1]);
set(gca,'ycolor',[1 1 1]);

% plot magnet
lw = 4;
R = 7;         % radius of circle
L = 5;
xL(1) = 5; xL(2) = 10;
yL(1) = 5;   yL(2) = 5;

yC = linspace(-yL(1), yL(1),200);
xC = sqrt(R^2 - yC.^2);

plot(xL,yL,'k','linewidth',lw);
hold on
plot(xL,-yL,'k','linewidth',lw);
plot(-xL,yL,'k','linewidth',lw);
plot(-xL,-yL,'k','linewidth',lw);
plot(xC,yC,'k','linewidth',lw);
plot(-xC,yC,'k','linewidth',lw);
text(-9,0,'N','fontsize',16);
text(8,0,'S','fontsize',16);

% plot arrow B field
plot([-3 3],[8 8],'b','linewidth',6);
plot([2 3],[9 7.8],'b','linewidth',6);
plot([2 3],[7 8.2],'b','linewidth',6);
text(4,7,'direction B-field centre N - S poles','fontsize',12)

% plot arrow force
plot([-8 -8],[7 10.5],'r','linewidth',4);
plot([-8.5 -8],[8 10.8],'r','linewidth',4);
plot([-7.5,-8],[8 10.8],'r','linewidth',4);
text(-20,9,'forces on current loop:','fontsize',12)
text(-20,8,'    +/- z direction','fontsize',12)

% plot current vs angle 0 to 6pi  (3 cycles)
xxi(1) = -10; xxi(2) = 16; yxi(1) = -10; yxi(2) = -10;  % x-axis
plot(xxi,yxi,'k');
xyi(1) = -10; xyi(2) = -10; yyi(1) = -7; yyi(2) = -13;  % y-axis
plot(xyi,yyi,'k');

xi = [-10 -8 -8 -4 -4 0 0 4 4 8 8 12 12 14];
yi = [-8 -8 -12 -12 -8 -8 -12 -12 -8 -8 -12 -12 -8 -8];
plot(xi,yi,'g');
text(14,-12,'angle (rev)','fontsize',12);
text(-20,-8,'current(side AB)','fontsize',12);
text(-2,-11,'1');
text(6,-11,'2');
text(14,-11,'3');

% plot xyz axes
plot([-19 -12],[0 0],'k');   % x
plot([-19 -19],[0 5],'k');   % z
plot([-19 -14],[0 3],'k');   % y
text(-11.8,0,'+X');
text(-14,3.5,'+Y');
text(-19.3,5.5,'+Z');

% plot current points - circle & diamond
hi = plot(12,4,'o');
set(hi,'MarkerSize',14,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[ 0 1 1]); 
text(9,2.5,'current +y direction');
hi = plot(12,-1.8,'d');
set(hi,'MarkerSize',14,'MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[ 1 0 1]); 
text(9,-3.3,'current -y direction');

hi = plot(-11,-8,'o');
set(hi,'MarkerSize',12,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[ 0 1 1]); 
hi = plot(-11,-12,'d');
set(hi,'MarkerSize',12,'MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[ 1 0 1]); 

% plot torque graph (tau)
xxi(1) = -10; xxi(2) = 16; yxi(1) = -20; yxi(2) = -20;  % x-axis
plot(xxi,yxi,'k');
xyi(1) = -10; xyi(2) = -10; yyi(1) = -21; yyi(2) = -15;  % y-axis
plot(xyi,yyi,'k');

tt = linspace(0,24,N);
dt = tt(2)-tt(1);
TT = 8;
tau = abs(4.*cos(2*pi*tt./TT))-20;
plot(tt-10,tau,'r')
text(14,-22,'angle (rev)','fontsize',12);
text(-15.2,-16,'net torque','fontsize',12);
text(-2,-21,'1');
text(6,-21,'2');
text(14,-21,'3');

% set up graph window
axis([-20 20 -25 10]);
axis equal
set(gca,'Xlim',[-20 20]);
set(gca,'xcolor',[1 1 1]);
set(gca,'ycolor',[1 1 1]);

M = getframe;
%M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,256,'nodither');  %RGB to indexed images
im(1,1,1,10) = 0;

for c = 1 : N   % *******************************************************
    
    % plot rotating armature
      x(1) = L * cos(-2*pi*t(c)/T);
      y(1) = L * sin(-2*pi*t(c)/T);
      x(2) = - x(1); y(2) = - y(1);
     plot(x,y,'g','linewidth',2*lw);
     
    
    % plot arrow for forces and current y direction
    % left side
       xa1(1) = x(2); xa1(2) = xa1(1); ya1(1) = y(2);
       if x(2) > 0
           ya1(2) = ya1(1) - 2;
           xaa1(1) = xa1(1) + 0.5; xaa1(2) = xa1(2);
           yaa1(1) = ya1(2) + 0.5; yaa1(2) = ya1(2);
           xaa2(1) = xa1(1) - 0.5; xaa2(2) = xa1(2);
           yaa2(1) = ya1(2) + 0.5; yaa2(2) = ya1(2);
           hR = plot(x(2),y(2),'o');
           set(hR,'MarkerSize',16,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[ 0 1 1]); 
       end;
       
       if x(2) <= 0
           ya1(2) = ya1(1) + 2; 
           xaa1(1) = xa1(1) - 0.5; xaa1(2) = xa1(2);
           yaa1(1) = ya1(2) - 0.5; yaa1(2) = ya1(2);
           xaa2(1) = xa1(1) + 0.5; xaa2(2) = xa1(2);
           yaa2(1) = ya1(2) - 0.5; yaa2(2) = ya1(2);
           hL = plot(x(2),y(2),'d');
           set(hL,'MarkerSize',16,'MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[ 1 0 1]); 
       end
      
     % right side
       xb1(1) = x(1); xb1(2) = xb1(1); yb1(1) = y(1);
       if x(1) > 0
           yb1(2) =  yb1(1) - 2;
           xbb1(1) = xb1(1) + 0.5; xbb1(2) = xb1(2);
           ybb1(1) = yb1(2) + 0.5; ybb1(2) = yb1(2);
           xbb2(1) = xb1(1) - 0.5; xbb2(2) = xb1(2);
           ybb2(1)  =  yb1(2) + 0.5;  ybb2(2) = yb1(2);
           hR = plot(x(1),y(1),'o');
           set(hR,'MarkerSize',16,'MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[ 0 1 1]); 
       end;
       
          if x(1) <= 0
           yb1(2) = yb1(1) + 2; 
           xbb1(1) = xb1(1) + 0.5; xbb1(2) = xb1(2);
           ybb1(1) = yb1(2) - 0.5; ybb1(2) = yb1(2);
           xbb2(1) = xb1(1) - 0.5; xbb2(2) = xb1(2);
           ybb2(1) = yb1(2) - 0.5; ybb2(2) = yb1(2);
           hL = plot(x(1),y(1),'d');
           set(hL,'MarkerSize',16,'MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[ 1 0 1]); 
       end  
       
       plot(xa1,ya1,'r','linewidth',4);
       plot(xaa1,yaa1,'r','linewidth',3);
       plot(xaa2,yaa2,'r','linewidth',3);
       
       plot(xb1,yb1,'r','linewidth',4);
       plot(xbb1,ybb1,'r','linewidth',3);
       plot(xbb2,ybb2,'r','linewidth',3);
       
           
    % plot point on current graph
    if t(c)/T < 0.25, yt = -8; end
    if t(c)/T > 0.25, yt = -12; end
    if t(c)/T > 0.75, yt = -8; end 
    if t(c)/T > 1.25, yt = -12; end 
    if t(c)/T > 1.75, yt = -8; end 
    if t(c)/T > 2.25, yt = -12; end 
    if t(c)/T > 2.75, yt = -8; end 
    plot(xt,yt,'go');
    xt = xt + dx;
    
    % plot point on torque graph
    tau1 = abs(4.*cos(2*pi*dt*c/TT))-20;
    plot(-10+24*c/N,tau1,'ro');
    
    % store image
    M = getframe; 
    %M = getframe(gcf);
    im(:,:,1,c) = rgb2ind(M.cdata,map,'nodither');
    
    pause(0.1)
    %pause
   
    % clear plotted lines
    plot(x,y,'w','linewidth',2*lw);
    plot(xa1,ya1,'w','linewidth',5);
    plot(xaa1,yaa1,'w','linewidth',5);
    plot(xaa2,yaa2,'w','linewidth',5);
    plot(xb1,yb1,'w','linewidth',5);
    plot(xbb1,ybb1,'w','linewidth',5);
    plot(xbb2,ybb2,'w','linewidth',5);
    hR = plot(x,y,'o');
    hL = plot(x,y,'d');
    set(hR,'MarkerSize',16,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[ 1 1 1]); 
    set(hL,'MarkerSize',16,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[ 1 1 1]);  
    
    
end

%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously

ag_name = 'animation.gif';
delay = 0.4;
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);


