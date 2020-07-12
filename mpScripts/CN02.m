%CN02.m

% Phase Diagrams for a Resistor, capcitor and Inductor
% Rotation of the phasors for current and voltage
% The animation can be saved as an animated gif
% ***  Uses the function arrow.m  ***

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% 180119

clear all
close all
clc
  
   N = 36;
   A = linspace(0,2*pi,N);   % arotation angle
   
% Save animation flagA = 0 (no save) or 1 (save)  [default flagA = 0]
   flagA = 0;
% File Name
   ag_name = 'ag_CN02.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.10;  
%  Frame counter start
    nt = 1;

% ======================================================================
%   GRAPHICS
% ======================================================================
f1 = figure(1);
   FS = 14;
   pos = [0.07 0.07 0.27 0.70];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 

   for c = 1:N
       subplot(3,1,1)
       phi = 0;
       axis off
       tm ='Resistor';
       xP = cos(A); yP = sin(A);
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       hold on
       xP = [-1 1]; yP = [0 0];
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       xP = [0 0]; yP = [-1 1];
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       
       arrow([0,0], [cos(A(c)),sin(A(c))],'Width',5, 'EdgeColor',[1 0 0], ...
           'FaceColor',[1 0 0],'TipAngle',20,'length',25);
       arrow([0,0], [0.7*cos(A(c)+phi),0.7*sin(A(c)+phi)],'Width',5, 'EdgeColor',[0 0 1], ...
           'FaceColor',[0 0 1],'TipAngle',30,'length',20)
       
       set(gca,'xLim',[-1 1]);
       set(gca,'yLim',[-1 1]);
       axis square
       hold off
       h_text = text(-0.3,1.2,tm);
       set(h_text,'color',[0 0 0],'fontsize',14);
       set(gca,'visible','off')
    
     
       
       subplot(3,1,2)
       phi = -pi/2; tm ='Capacitor';
       xP = cos(A); yP = sin(A);
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       hold on
       xP = [-1 1]; yP = [0 0];
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       xP = [0 0]; yP = [-1 1];
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       
       arrow([0,0], [cos(A(c)),sin(A(c))],'Width',5, 'EdgeColor',[1 0 0], ...
           'FaceColor',[1 0 0],'TipAngle',20,'length',25);
       arrow([0,0], [0.7*cos(A(c)+phi),0.7*sin(A(c)+phi)],'Width',5, 'EdgeColor',[0 0 1], ...
           'FaceColor',[0 0 1],'TipAngle',30,'length',20);
       
       set(gca,'xLim',[-1 1]);
       set(gca,'yLim',[-1 1]);
       axis square
       hold off
       h_text = text(-0.3,1.2,tm);
       set(h_text,'color',[0 0 0],'fontsize',14);
       set(gca,'visible','off')
       
        subplot(3,1,3)
       phi = pi/2; tm ='Inductor';
       xP = cos(A); yP = sin(A);
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       hold on
       xP = [-1 1]; yP = [0 0];
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       xP = [0 0]; yP = [-1 1];
       plot(xP,yP,'color',[0.6 0.6 0.6],'linewidth',0.5);
       
       arrow([0,0], [cos(A(c)),sin(A(c))],'Width',5, 'EdgeColor',[1 0 0], ...
           'FaceColor',[1 0 0],'TipAngle',20,'length',25);
       arrow([0,0], [0.7*cos(A(c)+phi),0.7*sin(A(c)+phi)],'Width',5, 'EdgeColor',[0 0 1], ...
           'FaceColor',[0 0 1],'TipAngle',30,'length',20)
       
%        xP = [0 cos(A(c))]; yP = [0 sin(A(c))];
%        plot(xP,yP,'r','linewidth',4);
%        xP = [0 0.7*cos(A(c) + phi)]; yP = [0 0.7*sin(A(c) + phi)];
%        plot(xP,yP,'b','linewidth',4);
%        
       set(gca,'xLim',[-1 1]);
       set(gca,'yLim',[-1 1]);
       axis square
       axis on
       h_text = text(-0.8,-1.5,'Current');
       set(h_text,'color',[1 0 0],'fontsize',14);
       h_text = text(0.2,-1.5,'Voltage');
       set(h_text,'color',[0 0 1],'fontsize',14);
       h_text = text(-0.3,1.2,tm);
       set(h_text,'color',[0 0 0],'fontsize',14);
       hold off
       set(gca,'visible','off')
       
       pause(0.001);
       
   if flagA == 1
       frame = getframe(f1);
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