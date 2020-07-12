% op_michA.m

% Michelson Interferometer
%   Plane Wave illumination of mirrors
%   Mirrors precisely aligned at right angles to the beam
%   Detector screen uniformly illuminated
%   Brightness of the screen depends upon the relative phase of the
%     reflections from the two mirrors
%   The intensity on the detector screen is shown as a function
%     the optical path difference between the two reflections

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mo/doc/op_michelson.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/


% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 190722   Matlab 2018b

clear
close all
clc


% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_AA.gif'; 
% Delay in seconds before displaying the next image  
    delay = 1.5;
    frame1 = 0;


% SETUP  ==============================================================
% Animation plot of screen intensity as optical path length varied
% Range for optical path length
  Nmax = 2;
  N = -Nmax:1/8:Nmax;
% Detector screen intensity
  SD = (sin(pi*N)).^2;


% ANIMATION GRAPHICS  =================================================

figure(1)
   pos = [0.1 0.2 0.25 0.35];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   
for cx = 1 : length(N)
  col = [SD(cx) 0 0];
  pos = [-1 -1 2 2];
  H_rect = rectangle('Position',pos,'Curvature',[1 1]);
  set(H_rect,'EdgeColor',col,'Facecolor',col);

  xlim([-1.1 1.1])
  ylim([-1.1 1.1])
  axis square
  axis off

  txt = sprintf('2 {\\Deltad} / {\\lambda} = %3.3f  \n',N(cx));
  title(txt)
  set(gca,'fontsize',14)
  pause(0.05)
 
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



% DETECTOR SCREEN INTENSITY ===========================================

num = 501;
N = linspace(-Nmax,Nmax, num);
SD = (sin(pi*N)).^2;

figure(2)
plot(N,SD,'r','linewidth',2);
grid on
set(gca,'fontsize',14)
xlabel('2 \Deltad / \lambda)')
ylabel('S_D  [a.u.]')
text(1.5,-0.135,'op_michA.m','fontsize',8)
