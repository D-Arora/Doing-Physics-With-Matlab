%mpSR02.m



clear
close all
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_SR.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.210;
    frame1 = 0;

% Steve: Fixed frame of reference (stationary)
% Mary: moving frame of reference w.r.t Steve

% Speed of light
   c = 3e8;
% Speed Mary w.r.t Steve   
   v = 0.8*c;
% Number of time intevals   
   N = 21;
% Distance from florr to ceiling
  yMax = 10;
% Time for light to travel from floor to ceiling  
  tMax = 2*yMax / c;
% Gamma
  G = 1/sqrt(1-v^2/c^2);  
% Mary's clock (proper time)
  tM = linspace(0,tMax,N);
% Steve's clocks (dilated time)   
  tS = G.*tM;
% y displacement  y = yM = yS
  y = c.*tM;
  y(12:N) = 10 - y(2:11);
  yMax = max(y);
% x displacement: Mary w.r.t Steve 
  x = v.*tS;
  xMax = max(x);
  dS = sqrt(xMax^2 + (2*yMax)^2);


% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.4 0.75];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   
for cc = 1 : N 
  subplot('position',[0.4,0.60,0.2,0.30])
     
 % Mary's frame of reference   
      xP = 0; yP = y(cc);
      Hplot = plot(xP,yP,'og');
      set(Hplot,'markerfacecolor','g','markersize',15)
      hold on
      
      plot([-1 1],[10 10],'k','linewidth',5)
      plot([-1 1],[0 0],'k','linewidth',5)
      set(gca,'xtick',[])
      
      xP = [0 0]; yP = [0 10];
      plot(xP,yP,'r','linewidth',1);
      
      txt = sprintf('Mary: \\Deltat_M = %3.2f ns   y_{Mmax} = %3.2f   \\Deltad_M = %3.2f ', ... 
                max(tM)*1e9, yMax, 2*yMax);
      H = title(txt); 
      set(H,'fontweight','normal','color','r','fontsize',14)
      
%       txt = sprintf('Mary:   \\Deltat_M = %3.1f ns ',tMax*1e9);
%       H = title(txt); 
%       set(H,'fontweight','normal','fontsize',14)
      
      txt = sprintf('t_M = %3.2f ns ',tM(cc)*1e9);
      H = text(-0.6,-3,txt);
      set(H,'fontweight','normal','color','k','fontsize',14)
      
      ylabel('y_M  [ m ]')
      set(gca,'ytick',[0 10])
      xlim([-1 1])
      ylim([-5 12])
      box on
      set(gca,'fontsize',12)
      pause(0.1)
      
      hold off
      
      
   % Steve's frame of reference 
    subplot('position',[0.2,0.1,0.6,0.40])
      xP = x(cc); yP = y(cc);
      Hplot = plot(xP,yP,'og');
      set(Hplot,'markerfacecolor','g','markersize',15)
      hold on
      
      plot([-1 30],[10 10],'k','linewidth',5)
      plot([-1 30],[0 0],'k','linewidth',5)
      
      xP = [0 x(11) x(N)]; yP = [0 10 0];
      plot(xP,yP,'b','linewidth',1)
      
      xP = x(cc); yP = 10;
      Hplot = plot(xP,yP,'sr');
      set(Hplot,'markerfacecolor','r','markersize',13)
           
      txt = sprintf('Steve: \\Deltat_S = %3.2f ns   x_S = %3.2f   y_{Smax} = %3.2f   \\Deltad_S = %3.2f ', ... 
                max(tS)*1e9, xMax, yMax, dS);
      H = title(txt); 
      set(H,'fontweight','normal','color','b','fontsize',14)
      
      txt = sprintf('t_S = %3.2f ns ',tS(cc)*1e9);
      H = text(1,-5,txt);
      set(H,'fontweight','normal','color','k','fontsize',14)
      
      txt = sprintf('t_S = %3.2f ns',tS(cc)*1e9);
      H = text(19,-5,txt);
      set(H,'fontweight','normal','color','k','fontsize',14)
      
      txt = 'EVENT occurs in F of R of Mary';
      H = text(2,18,txt);
      set(H,'fontweight','normal','color','r','fontsize',14)
      
      txt = 'Steve observes the EVENT ';
      H = text(2,15,txt);
      set(H,'fontweight','normal','color','b','fontsize',14)
      
      txt = '  occurring in F of R of Mary';
      H = text(2,13,txt);
      set(H,'fontweight','normal','color','b','fontsize',14)
      
      txt = '  mpSR02.m';
      H = text(25,-15,txt);
      set(H,'fontweight','normal','color','k','fontsize',8)
      
      xlabel('x_S  [ m ]');
      ylabel('y_S  [ m ]')
      set(gca,'ytick',[0 10])
      xlim([0 30])
      ylim([-10 20])
      box on
      set(gca,'fontsize',12)
      axis square
      pause(0.1)
      hold off   
           
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
