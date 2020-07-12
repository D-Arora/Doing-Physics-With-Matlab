% oscC00NC.m
% 190524  Matlab 2018b

% SIMULATION OF THE MOTION OF N COUPLED OSCILLATORS
% DYNAMICS OF CRYSTAL LATTICES
% All parameters are changed within the Script
% [default values and units]

% DOING PHYSICS ONLINE 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/oscC00NA.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney


clear
close all
clc

% Animated gif  =======================================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flag2
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_osc.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;


% SETUP ===============================================================

% Propagation constant: waves 1 and 2 [units  pi/d: dk <<  k]
   k = 11.35;
   dk = 1.15;

% Time span
   tN = 61;
   tMax = 2;
   
   
% spring costants kS  [10]
   kS = 10;
% mass of atoms or particles m  [0.1]
   m = 0.1;
% separation distance d bewtween atoms or particles  [1]   
   d = 1;
% Parameters for X axis
   xMax = 4*d; xN = 1501;
   
% CALCULATIONS ========================================================
% Propgation constants  /  Wavelengths  wL1 wL2
  k1 = k*pi/d; k2 = (k + dk)*pi/d;
  wL1= 2*pi/k1;  wL2 = 2*pi/k2;
% amplitudes
  A1 = 1; A2 = 1;

% X axis / time span
   x = linspace(0,xMax,xN);
   t = linspace(0,tMax,tN);

% Angular frequencies omega
   omega1 = sqrt(4*kS/m) * abs(sin(k1*d/2));
   omega2 = sqrt(4*kS/m) * abs(sin(k2*d/2));
   domega = abs(omega1 - omega2);
   
% v0
  v0 = d*sqrt(kS/m);   
% Phase velocity
  vP1 = omega1/k1;
  vP2 = omega2/k2;
% Group velocity
  vG = abs( (omega1 - omega2)/(k1-k2) );
 % Displacement
  e1 = zeros(xN,tN);      % wave 1
  e2 = zeros(xN,tN);      % wave 2 
  eEnv = zeros(xN,tN);    % envelope
  
  
for cx = 1:xN
 for ct = 1:tN
       e1(cx,ct) = A1*cos( omega1*t(ct) - k1*x(cx) );
       e2(cx,ct) = A2*cos( omega2*t(ct) - k2*x(cx) );
       eEnv(cx,ct) = (A1+A2)*cos( 0.5*abs((omega1-omega2)*t(ct)) - 0.5*abs((k1-k2)*x(cx)) ) ;     
 end
end
 
 if vP1 < 1e-6
     e1(:,:) = 0;
     vG = 0;
 end

 if vP2 < 1e-6
     e2(:,:) = 0;
     vG = 0;
 end
 
 e = e1 + e2;



% GRAPHICS ============================================================ 

 figure(1)

   pos = [0.52 0.2 0.35 0.6];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   
for ct = 1: tN
   
   subplot('position',[0.1 0.4 0.8 0.5])
   % Sine wave 1
     xP = x(1:end); yP = e1(1:end,ct);
     plot(xP,yP,'r','linewidth',1.5)
     hold on
     xP = t(ct)*vP1; yP = e1(1,1);
     hPlot = plot(xP,yP,'or');
     set(hPlot,'markersize',9,'markerfacecolor','r','markeredgecolor','k')
       
   % Sine wave 2
     xP = x(1:end); yP = e2(1:end,ct);
     plot(xP,yP,'m','linewidth',1.5)
     hold on
     xP = t(ct)*vP2; yP = e2(1,1);
     hPlot = plot(xP,yP,'om');
     set(hPlot,'markersize',9,'markerfacecolor','m','markeredgecolor','k')

  if vG ~= 0   
  % Resultant wave
     xP = x; yP = e(:,ct);
     plot(xP,yP,'b','linewidth',2)
     xP = t(ct)*vG; yP = 0;
     hPlot = plot(xP,yP,'ob');
     set(hPlot,'markersize',9,'markerfacecolor','b')
     
     xP = x; yP = eEnv(:,ct);
     plot(xP,yP,'b','linewidth',1)
  end
     
    ylim([-2.1 2.1])
    xlim([0 xMax])
    grid on
    set(gca,'fontsize',12)
    xlabel('x ')
    ylabel('left   e   right')
    tm = [' t = ' num2str(t(ct),'%2.1f ')];
    title(tm)
    hold off
% =====================================================================  

   subplot('position',[0.1 0.05 0.85 0.25])
   
   xlim([0 100])
   ylim([0 100])
   
   fs = 14;
   txt = ['d = ' num2str(d,'%2.1f') '   m = ' num2str(m,'%2.1f') '   k_S = ' num2str(kS,'%2.1f') ];
   Htext = text(2,99,txt);
   set(Htext,'fontsize',fs)
   
   txt = ['Wave 1  \lambda_1 = ' num2str(wL1,'%2.2f') '   k_1 = ' num2str(k,'%2.2f') '\pi/d' ...
      '   \omega_1 = ' num2str(omega1,'%2.2f') '   v_{P1} = ' num2str(vP1,'%2.2f')];
   Htext = text(2,70,txt);
   set(Htext,'fontsize',fs,'color','r')
   
   txt = ['Wave 2  d\lambda_2 = ' num2str(wL2,'%2.2f') '   k_2 = ' num2str((k+dk),'%2.2f') '\pi/d' ...
       '   \omega_2 = ' num2str(omega2,'%2.2f')   '   v_{P2} = ' num2str(vP2,'%2.2f')];
   Htext = text(2,40,txt);
   set(Htext,'fontsize',fs, 'color','m')
   
   txt = ['Resultant wave    v_G = ' num2str(vG,'%2.2f')];
   Htext = text(2,10,txt);
   set(Htext,'fontsize',fs, 'color','b')
   
   set(gca,'fontsize',fs)
   axis off
  % ======================================================================
  
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
  
   pause(0.01)
   
end


