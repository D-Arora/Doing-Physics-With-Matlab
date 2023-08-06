% chaos23.m



clear
close all
clc

% INPUTS ==========================================================

col = [0 0 1];
%col = [1 0 0];
x0 = 0.1 ;
y0 = 0.1 ;
z0 = 0.1;
t1 = 0;
t2 = 20;

sigma = 10;
beta = 8/3;

rho = 13.926;  %28

% ANIMATION  =======================================================
  flagP = 0;
% flagP = 0  animation not displayed / flagG = 1 animation displayed
 %  flagG = 1;
% f_gif = 0 animated gif NOT saved / f_gif = 1 file saved
   f_gif = 1;
% File name for animated gif   
   ag_name = 'ag_chaos23.gif';
% Delay in seconds before displaying the next image  
   delay = 0.0;  
% Frame counter start
   nt = 1; 

% SETUP =============================================================
tSpan = [t1 t2];
u0 = [x0;y0;z0];

K = [sigma,beta,rho];

[t, SOL] = ode45(@(t,u) FNode(t,u,K), tSpan, u0);

x = SOL(:,1); y = SOL(:,2); z = SOL(:,3);

% CRITICAL POINTS  ================================================
eta = sqrt(beta*(rho-1));
cP1(3) = rho-1;
cP1(2) = eta;
cP1(1) = eta;
cP2(3) = rho-1;
cP2(2) = -eta;
cP2(1) = -eta;
 
% GRAPHICS  ========================================================
  FS = 14;

if flagP == 1  

figure(3)
  pos = [0.68 0.05 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  Hplot3 = plot3(x(1),y(1),z(1),'og');
  set(Hplot3,'markerfacecolor','g','MarkerSize',10)
  hold on
 % Hplot3 = plot3(x(end),y(end),z(end),'o');
 % set(Hplot3,'markerfacecolor',col','MarkerSize',10)
  hold on
 for c = 1 : 10: length(x)
    plot3(x(1:c),y(1:c),z(1:c),'b')
    hold on
    Hplot3 = plot3(x(c),y(c),z(c),'o');
    set(Hplot3,'markerfacecolor',[1 0 0]','MarkerSize',10)
    grid on
    box on
    xlabel('x'); ylabel('y'); zlabel('z')
    xlim([-30 30]); ylim([-30 30]);zlim([0 50])
    set(gca,'FontSize',FS)
    view(25,15)
    pause(0.001)
    hold off

   pause(0.1)
      
      if f_gif > 0 
         frame = getframe(3);
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
 end
end



figure(1)
  pos = [0.05 0.05 0.23 0.5];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
  
subplot(3,1,1)
  plot(t,x,'color',col,'LineWidth',2)
  ylabel('x')
  grid on
  set(gca,'FontSize',FS)
  hold on
subplot(3,1,2)
  plot(t,y,'color',col,'LineWidth',2) 
  ylabel('y')
  grid on
  set(gca,'FontSize',FS)
  hold on
subplot(3,1,3)
  plot(t,z,'color',col,'LineWidth',2) 
  xlabel('t')
  ylabel('z')
  grid on
  set(gca,'FontSize',FS)
  hold on

%%
figure(2)
  pos = [0.38 0.05 0.23 0.23];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  Hplot3 = plot3(x(1),y(1),z(1),'og');
  set(Hplot3,'markerfacecolor','g','MarkerSize',10)
  hold on
  Hplot3 = plot3(x(end),y(end),z(end),'o');
  set(Hplot3,'markerfacecolor',[1 0 0]','MarkerSize',10)
  hold on
  plot3(x,y,z,'color',col)

  Hplot3 = plot3(cP1(1),cP1(2),cP1(3),'om');
  set(Hplot3,'markerfacecolor',[1 0 1]','MarkerSize',8)
  Hplot3 = plot3(cP2(1),cP2(2),cP2(3),'om');
  set(Hplot3,'markerfacecolor',[1 0 1]','MarkerSize',8)
  grid on
  box on
  xlabel('x'); ylabel('y'); zlabel('z')
  set(gca,'FontSize',FS)
  view(25,15)

%%

% FUNCTIONS ==========================================================

function du = FNode(t,u,K)
   sigma = K(1); beta = K(2); rho = K(3);
   x = u(1); y = u(2); z = u(3);
   du(1) = sigma*(y - x);
   du(2) = x*(rho - z) - y;
   du(3) = x*y - beta*z;
   du = du';
end






