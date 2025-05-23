% p002.m    

% PHOTONICS and WAVES  Refraction at an interface
% Snell's Law
% Critical angle

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 200418 / Matlab version R2020a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/p002.htm
% Download Location
%    http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
  

clear
close all
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1; 
% file name for animated gif   
   ag_name = 'ag_p002A.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.20;
    frame1 = 0;

% Refractive indices    
n1 = 1; n2 = 1.33;
% Angle of incidence
A2 = 0:1:89;
N = length(A2);
x2 = sind(A2);
y2 = cosd(A2);

% Angle refraction
A1 = asind((n2/n1).*sind(A2));
x1 = -sind(A1);
y1 = -cosd(A1);

% Crtitical angle
  Acrit = asind(n1/n2);


% GRAPHICS  ===========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.28 0.40]);
  
  axis off
  
for c = 1:N
  AC = linspace(0,180,360);
  xC = cosd(AC);
  yC = sind(AC);
  xP = xC; yP = yC;
  plot(xP,yP,'k','linewidth',3)
  hold on
  fill(xP,yP,[0.8 0.8 1])

  AC = linspace(180,360,360);
  xC = cosd(AC);
  yC = sind(AC);
  xP = xC; yP = yC;
  plot(xP,yP,'k','linewidth',3)
  hold on
  fill(xP,yP,[1 1 0.7])
  
 % normal
    xP = [0 0]; yP = [-1.2 1.2];
    plot(xP,yP,'k','linewidth',2)
    
 % Incident ray 
   xP = [0 x2(c)]; yP = [0 y2(c)];
   plot(xP,yP,'r','linewidth',2)
 % Reflected ray
   LW = 1;  
   if A2(c) > Acrit; LW = 2; end
   xP = [0 -x2(c)]; yP = [0 y2(c)];
   plot(xP,yP,'r','linewidth',LW)
 % Refracted ray 
   if A2(c) < Acrit
    xP = [0 x1(c)]; yP = [0 y1(c)];
    plot(xP,yP,'r','linewidth',2)
   end
   
   txt = sprintf('n_1 = %2.2f',n1);
   text(0.5,-0.2,txt,'fontsize',12)
   
   txt = sprintf('n_2 = %2.2f',n2);
   text(-0.7,0.2,txt,'fontsize',12)
  
   text(0.2,1.1, 'Incident Ray','fontsize',14,'FontSmoothing','off')
   text(-1,1.1,'Reflected Ray','fontsize',14,'FontSmoothing','off')
   
   if A2(c) > Acrit
     text(-1.1,-1.1, 'NO','fontsize',14)
   else
     text(-1,-1.1, '   Refracted Ray','fontsize',14,'FontSmoothing','off')
   end
   
   title('Snells Law  n_1 sin\theta_1 = n_2 sin\theta_2','FontSmoothing','off')
   
   txt = sprintf('\\theta_C = %2.1f^o',Acrit);
   text(0.06,0.7,txt,'fontsize',12)
   
   txt = sprintf('\\theta_1 = %2.1f^o',A1(c));
   text(-0.5,-0.7,txt,'fontsize',12)
   
   txt = sprintf('\\theta_2 = %2.1f^o',A2(c));
   text(0.05,0.5,txt,'fontsize',12)
   axis equal
   pause(0.1)
   set(gca,'fontsize',12)
   axis off
   set (gca, 'FontSmoothing', 'off')
   
   % Save animation 
    if flagAG == 1
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
end
  
figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.3 0.05 0.30 0.35]);
  set(gcf,'color','w');
  
  plot([0 90],[0 90],'r')
  hold on
  plot([Acrit Acrit],[0 90],'-r')
  
  A = A2(A2<Acrit); LA = length(A)+1;
  plot(A2(1:LA),real(A1(1:LA)),'b','linewidth',2)
    
  xlabel('\theta_2 (incident angle)')
  ylabel('\theta_1 (refracted angle)')
  txt = sprintf('2 \\rightarrow 1: n_2 = %1.2f  >  n_1 = %1.2f',n2,n1);
  title(txt)
  text(90,-10,'p002A.m','fontsize',8)
  txt = sprintf('\\theta_{critical} = %3.2f^o',Acrit);
  text(50,85,txt)
  
  xlim([0 90])
  ylim([0 90])
  grid on
  set(gca,'xtick',0:10:90)
  set(gca,'ytick',0:10:90)
  set(gca,'fontsize',12)
  axis square  
  
  
  