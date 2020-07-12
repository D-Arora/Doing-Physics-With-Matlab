% mpSR07.m
% 190627  Matlab R2018b

% SPECIAL RELATIVITY ANIMATIONS
% LENGTH CONTRACTION
% Gamma   1/GAMMA
% Moving frame of reference (RED)
% Stationary frame of reference (BLUE) 
% Moving rod as observed by stationary observer

% DOING PHYSICS ONLINE 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ONLINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Documentation
%    http://www.physics.usyd.edu.au/teach_res/hsc/spP/mod7/m72D.htm

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

close all
clear
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_SR.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;


% LENGTH CONTRACTION  =================================================
% Proper length
  L0 = 80;
% Grid points
  nT = 51;
% Time 
  t = linspace(0,200,nT);
% Velocity of moving systems wrt stationary system
  v = [0 0.6 0.8 0.9 0.95 0.999];
% 1/gamma
  gammaT = sqrt(1-v.^2);
% Contracted length
  L = gammaT.*L0;
% Width of line for rod
  W = 10;

figure(1)
 pos = [0.1 0.2 0.3 0.6];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   xlim([-110 110])
   ylim([0 100])
   box on
   grid on
   
for cc = 1 : nT
        
     xP = [-L0/2 L0/2]; yP = [88 88];
     plot(xP,yP,'b','linewidth',W);
        
     hold on
     n = 2;
     x = -90 + v(n).*t(cc);
     xP = [x-L(n), x]; yP = [68, 68];
     plot(xP,yP,'r','linewidth',W);
  
     n = 3;
     x = -90 + v(n).*t(cc);
     xP = [x-L(n), x]; yP = [53, 53];
     plot(xP,yP,'r','linewidth',W);
  
     n = 4;
     x = -90 + v(n).*t(cc);
     xP = [x-L(n), x]; yP = [38, 38];
     plot(xP,yP,'r','linewidth',W);
   
     
     n = 5;
     x = -90 + v(n).*t(cc);
     xP = [x-L(n), x]; yP = [23, 23];
     plot(xP,yP,'r','linewidth',W);
        
     n = 6;
     x = -90 + v(n).*t(cc);
     xP = [x-L(n), x]; yP = [8, 8];
     plot(xP,yP,'r','linewidth',W);
  
     
     % ----------------------------------------------------------------
     n = 1;
     txt = sprintf('  v/c = %2.2f   L_0 = %2.2f  '  ,v(n),L(n) );
     Htext = text(-100,95,txt);
     set(Htext,'EdgeColor','b','linewidth',3,'fontsize',14)
     
     n = 2;
     txt = sprintf('  v/c = %2.2f   L = %2.2f  '  ,v(n), L(n) );
     Htext = text(-100,75,txt);
     set(Htext,'EdgeColor','r','linewidth',3,'fontsize',14)
     
     n = 3;
     txt = sprintf('  v/c = %2.2f   L = %2.2f  '  ,v(n), L(n) );
     Htext = text(-100,60,txt);
     set(Htext,'EdgeColor','r','linewidth',3,'fontsize',14)
     
     n = 4;
     txt = sprintf('  v/c = %2.2f   L = %2.2f  '  ,v(n), L(n) );
     Htext = text(-100,45,txt);
     set(Htext,'EdgeColor','r','linewidth',3,'fontsize',14)
     
    n = 5;
     txt = sprintf('  v/c = %2.2f   L = %2.2f  '  ,v(n), L(n) );
     Htext = text(-100,30,txt);
     set(Htext,'EdgeColor','r','linewidth',3,'fontsize',14) 
   
   n = 6;
     txt = sprintf('  v/c = %2.2f   L = %2.2f  '  ,v(n), L(n) );
     Htext = text(-100,15,txt);
     set(Htext,'EdgeColor','r','linewidth',3,'fontsize',14) 
     
     text(90,0,'mpSR07.m','fontsize',8);
     
     xlim([-110 110])
     ylim([0 100])
     axis off
     pause(0.2)
     
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
    
    hold off
end
   
   
% GAMMA  1/ GAMMA =====================================================

% Velocity
  v = linspace(0,0.99,201);
% Gamma  
  gamma = 1./sqrt(1-v.^2);
  gammaT = 1./gamma;

figure(2)
   pos = [0.5 0.2 0.3 0.2];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   
   xP = v;  yP = gamma;    
   plot(xP,yP,'b','linewidth',2);
   hold on
   yP = gammaT;
   plot(xP,yP,'r','linewidth',2);

  grid on
  set(gca,'fontsize',12)
  legend('\gamma','1/\gamma','location','north','orientation','horizontal')
  xlabel('v / c')
  text(0.85,-1.2,'mpSR07.m','fontsize',8);
     
