% mpSR08.m
% 190628  Matlab R2018b

% SPECIAL RELATIVITY ANIMATIONS
% SNAKE MOVING and observed in two frames of reference
% Gamma   1/GAMMA
% Moving frame of reference (snake RED)
% Stationary frame of reference (BLUE) 


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
    delay = 0.25;
    frame1 = 0;


% SETUP  =================================================
 
% Speed of light  [m/s]
  c = 3e8;
% Proper length of snake [m]
  L0 = 1.00;
% Grid points
  nT = 51;
% Time [ns]
 tMin = -5.76e-9; tMax = 12e-9;
  tS = linspace(tMin,tMax,nT);
% Velocity of snake wrt stationary system
  v = 0.6*c;
% 1/gamma
  gammaT = sqrt(1-(v/c)^2);
% gamma
  gamma = 1/gammaT;
% Contracted length
  L = gammaT.*L0;

% Stationary frame of reference (Steve)  ------------------------------

% Position of snake's tail and head 
  xST = -1.00 + v.*(tS-tMin);
  xSH = xST + L;
% Position of knives K: left L / right R
  xKL = [0 0]; xKR = [1.00 1.00];
  yKL = [2 6]; yKR = [2 6];
  
% Moving frame of reference (snake - Mary)
% Position of knifes


% Position of snake
  xSnake = [0 1]; ySnake = [0 0];
% Position of knifes   
  xM = gamma .* (0 - v.*tS);
%  tM = gammaT .* tS;         ??????
   
% GRAPHICS  ===========================================================
figure(1)
 pos = [0.1 0.2 0.3 0.4];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   xlim([-1.10 2.1])
   ylim([-10 10])
   box on
   grid on
   
for cc = 1 : nT
     
subplot(2,1,1)        
     xP = [xST(cc) xSH(cc)]; yP = [0 0];
     plot(xP,yP,'g','linewidth',10);
     xlim([-1.10 2.10])
     ylim([-10 10])   
     hold on
     
     xP = xKL; yKL = [2 6];
     
   if xST(cc) > 0
        yKL = [-6 -2];
        yKR = [-6 -2];
        H = text(0,3,'0.00 ns');
        set(H,'color','b','edgecolor','b')
        H = text(1,3,'0.00 ns');
        set(H,'color','b','edgecolor','b')
   end
   
     yP = yKL;
     plot(xP,yP,'r','linewidth',10);
     xlim([-1.10 2.10])
     ylim([-10 10])
     
     xP = xKR; yP = yKR;
     plot(xP,yP,'r','linewidth',10)
     xlim([-1.10 2.10])
     ylim([-10 10])
     
     set(gca,'fontsize',12)
     set(gca,'ytick',[])
     txt = sprintf('  t_S = %2.1f   ns  '  ,tS(cc)*1e9);
     title(txt,'color','b') 
     xlabel('x_S   [m]')
     text(-1,8.5,'Steve observations','color','b','fontsize',14);
     hold off
     
 subplot(2,1,2)        
     xP = [xST(cc) xSH(cc)]; yP = [0 0];
     plot(xP,yP,'g','linewidth',10);
     xlim([-1.10 2.10])
     ylim([-10 10])   
    
     xP = xSnake; yP = ySnake;
     plot(xP,yP,'g','linewidth',10);
     xlim([-1.10 2.10])
     ylim([-10 10])   
     hold on
     
     xP = [xM(cc) xM(cc)]; yP = [2 6];
     if xM(cc) < 0; yP = [-6 -2];
        H = text(0,3,'0.00 ns');
        set(H,'color','r','edgecolor','r')
      end
     plot(xP,yP,'r','linewidth',10);
     
     xP = [xM(cc)+L0 xM(cc)+L0]; yP = [2 6];
     if xM(cc)+L0 < 1.25; yP = [-6 -2]; disp(cc);
        H = text(1.25,3,'-2.50 ns');
        set(H,'color','r','edgecolor','r')
     end
     plot(xP,yP,'r','linewidth',10);
     
     xlim([-1.10 2.10])
     ylim([-10 10])   
     
     set(gca,'fontsize',12)
     set(gca,'ytick',[])
%    Not sure how to calculate time in snake's frame ????     
%      tM = gamma * (tS(cc) - v/c^2);
     tM = -8.3111 - 1e9*(xM(cc) - xM(1))/v;
     txt = sprintf('  t_M = %2.1f   ns  '  ,tM);
      title(txt,'color','r') 
     xlabel('x_M   [m]')
     text(-1,8.5,'Snake observations','color','r','fontsize',14);
     
     text(1.5,-15,'mpSR08.m','fontsize',8);
 

      pause(0.8)
     
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
   
   
