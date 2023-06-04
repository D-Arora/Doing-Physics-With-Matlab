% zHOWTO.m


% #1   plot



%% 1  plot
clear; close all; clc

  N = 299; x1 = 0; x2 = 100; 
  x = linspace(x1,x2,N);
  wL = 25;
  y1 = sin(2*pi*x/wL);
  y2 = sin(2*pi*x/wL + pi/4);




% GRAPHICS ===============================================================
figure(1) 
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.1 0.2 0.25]);
  set(gcf,'color','w');
  FS = 14;
  txtX = 'x';
  txtY = 'y';
  txtT = 'sine function';
  leg1 = 'y1'; leg2 = 'y2';
  XLIM = [0 max(x)];
  XTICKS = 0:20:100;
  YLIM = [-1.1 1.8];
  YTICKS = -1:0.5:1;
  
  xP = x; yP = y1;
  plot(xP,yP,'b','LineWidth',2)
  hold on
  yP = y2;
  plot(xP,yP,'r','LineWidth',2)
  xlim(XLIM); ylim(YLIM)
  xticks(XTICKS); yticks(YTICKS);
  grid on; box on; 
  xlabel(txtX); ylabel(txtY)
  title(txtT,'FontWeight','normal')
  legend(leg1,leg2,'Orientation','horizontal','Location','north')
  set(gca,'fontsize',14)


%% #2   ANIMATION --> gif file
  clear; close all; clc

  Nx = 299; x1 = 0; x2 = 100; 
  x = linspace(x1,x2,Nx);
  Nt = 199; t1 = 0; t2 = 90;
  t = linspace(t1,t2,Nt);
  wL = 20;  T = 30;
  k = 2*pi/wL; w = 2*pi/T;

  y = zeros(Nt,Nx);

for c = 1:Nt
    y(c,:) = sin(k*x - w*t(c));
end

% ANIMATION SETUP =======================================================
% (0 no)  (1 yes) for flag1
   flag1 = 1;    
% file name for animated gif   
    ag_name = 'ag001.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.01; 
% Frame to start
    frame1 = 0;

% GRAPHICS ===============================================================
figure(1) 
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.1 0.2 0.25]);
  set(gcf,'color','w');
  FS = 14;
  txtX = 'x';
  txtY = 'y';
  txtT = 'sine function';
 
    XLIM = [0 max(x)];
    XTICKS = 0:20:100;
    YLIM = [-1.1 1.1];
%   YTICKS = -1:0.5:1;
 for c = 1:Nt 
  xP = x; yP = y(c,:);
  plot(xP,yP,'b','LineWidth',2)
  xlim(XLIM); ylim(YLIM)
  xticks(XTICKS);
  grid on; box on; 
  xlabel(txtX); ylabel(txtY)
  title(txtT,'FontWeight','normal')
  set(gca,'fontsize',14)

    if flag1 > 0
         frame1 = frame1 + 1;
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
      % On the first loop, create the file. In subsequent loops, append.
         if frame1 == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
    end 
  pause(0.01)
 end


 %% #3
 clear; clc; close all
 syms a b c
 eq1 = a+b+c -2;
 eq2 = a-b-c;
 eq3 = 2*a+c-c+1;

 [z1 z2 z3] = vpasolve(eq1,eq2,eq3, [a,b,c]);
% z1 = z.a; z2 = z.b; z3 = z.c;
 z1;
 z2;
 z3;

 fprintf('a = %2.2f \n',z1)
 


 