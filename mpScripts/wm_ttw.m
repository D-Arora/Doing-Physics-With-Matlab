% wm_ttw.m

% 171214
% Animation of a travlling [1D] sinusoidal wave propagating in Z direction
% Animation of the motion can be saved as an animated gif or avi. file



% John A Sims 
% email: john.sims@ufabc.edu.br
% Universidade Federal do ABC 
% http://www.ufabc.edu.br/
% Biomedical Engineering Department,
% Centro de Engenharia, Modelagem e Ciencias Sociais Aplicadas.
% (Centre for Engineering, modelling and applied social sciences)
% Rua Arcturus (Jd Antares)
% Anchieta
% 09606070 - Sao Bernardo do Campo, SP - Brasil

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm


clear all
close all
clc


% =======================================================================
% TRAVELLING WAVE INPUT PARAMETERS
%   S.I. units for all quantities
%   Default values given in brakets
% =======================================================================
   e = 1;       %  left (-1) or right (+1)
   A = 10;       %  amplitude  (10)
   L = 25;       %  wavelength  (25)
   T = 20;       %  period  (20)
   nT = 4;       %  no. periods for animation  (4)
   nX = 3;       %  no. wavelengths displayed (3)
   Nx = 2076;    %  Number of grid X grid points (2076)
   Nt = 200;     %  number of time steps (200)

% Save animation fl agA = 0 or 1 or 2  (defualt flagA = 0) ---------------
   flagA = 0;
   % flagA = 0   Animation NOT saved
   % flagA = 1   Animation saved as an animated gif file
   % flagA = 2   Animation saved as an AVI file
   

% =======================================================================   
% Setup for saving images animated gif file and AVI files
% =======================================================================
switch flagA
case 1
%  ANIMATED GIF:
     ag_name = 'ag_wm_TWLR.gif';
      if flagA == 1
       if e == -1
         ag_name = 'ag_wm_TWLR.gif' ;
       end   % if
       if e == +1
         ag_name = 'ag_wm_TWRL.gif';
       end   %  if
      end
%  Delay in seconds before displaying the next image  
    delay = 0.20;  
%  Frame counter start
    nt = 1;

case 2   
%  AVI: open writer object for avi.
   aviName = 'transWave_LR.avi';
   if flagA == 2
       if e == -1
         aviName = 'avi_TW_LR.avi';
       end   % if
       if e == +1
         aviName = 'avi_TW_RL.avi';
       end   %  if
     aviObj = VideoWriter(aviName);
     aviObj.FrameRate = 10;  % from 30
     aviObj.Quality = 40;    % from 75
     open(aviObj);
   end   % if
end   %  switch
   
% =======================================================================   
%   CALCULATIONS
% =======================================================================
   w = 2*pi / T;    %  angular frequency 
   k = 2*pi / L;    %  wave number k
   tMax = nT * T;   %  max time interval for animation
   zMax = nX * L;   %  max Z distance
   
   t = linspace(0,tMax,Nt);   % grid points for time 
   z = linspace(0,zMax,Nx);   % grid points for Z axis
   
    
%  Wave function s(z,t) ---------------------------------------------
%    Calculate all data and store it in S, rows are x, cols t. 
%    Later, phase relationship between all particles in an instant 
%    of time (i.e spatial domain information) will be taken from , 
%    rows of S and time (temporal domain) for a single particle 
%    taken from columns.
   s = zeros(Nx,Nt); 
   
% Time loop   
   for c = 1 : Nt
     s(:,c) = A .* sin(k.*z + e*w*t(c));
   end

   
% =======================================================================     
%   GRAPHICS 
% =======================================================================  
f1 = figure(1) ;
   pos = [0.1 0.1 0.4 0.8];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   NS = 20; 
   yT1 = zeros(Nt,1);   xT1 = zeros(Nt,1);
   yT2 = zeros(Nt,1);   xT2 = zeros(Nt,1);
   
   
 %  Adjust to reduce/increase the frequency of red circles. 
 %    compared to blue circles.  i.e. Every "NR" is to be red.
   NR = 100; 
   
 %  Define text for three figures based on user choice of e.  
 %    Figure 1 - spatial domain - main title.
 %    Figure 2 - temporal domain - blacktitle (it is the black dot)
 %    Figure 3 - temporal domain - cyantitle (guess why).
   blacktitle= sprintf('Wavefunction of black particles through time.');    
   cyantitle = sprintf('Wavefunction of magenta particles through time.');
   
   switch e
   case 1
     eqnname = 'y(z,t) = 10 sin(k z + \omega  t)';
     maintitle  = sprintf('Phase relation of particles in space (velocity -z direction)');
     blacktitle = sprintf('%s\n y(t) = 10 sin(\\phi_b + \\omega t) ',blacktitle);
     cyantitle  = sprintf('%s\n y(t) = 10 sin(\\phi_m + \\omega t) ',cyantitle); 
   
   case -1
     eqnname = 'y(z,t) = 10 sin(k z - \omega  t)';
     maintitle  = sprintf('Phase relation of particles in space (velocity +z direction)');
     blacktitle = sprintf('%s\n y(t) = 10 sin(\\phi_b - \\omega t) ',blacktitle);
     cyantitle  = sprintf('%s\n   y(t) = 10 sin(\\phi_m - \\omega t) ',cyantitle); 
   
   otherwise %error condition
   error('Error: Value of e must be either 1 or -1!')   
   end
   
   maintitle = sprintf('%s\n %s\n ',maintitle,eqnname);
   %maintitle = sprintf('%s\n %s\n \\fontsize{6}{john.sims@ufabc.edu.br 2017, modified from code by ian.cooper@sydney.edu.au}',maintitle,eqnname);
   blacktitle = sprintf('%s where \\phi_b = k z_b  Initial phase for black = 0',blacktitle);
   cyantitle = sprintf('%swhere \\phi_m = k z_m    Initial phase for magenta \\neq 0',cyantitle);


 nt = 1;
 ct = 0:Nt-1;  % create time step vector
 
for c = ct+1
   
%  wavefunction vs postion plot -----------------------------------------    
   hold off
   xP = z(1:NS:end); yP = s((1:NS:end),c);
   pos1 = [0.15 0.6 0.7 0.3];
   ax1 = subplot('Position',pos1);
   
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',6,'MarkerEdgeColor',[0 0 1]);%, ...
      % 'MarkerFaceColor', [0 0 1]); % How to present blue circles
   hold(ax1, 'on') %Equivalent of hold on when using subplots.

   title(maintitle,'fontweight','normal');
   xlabel('position  z  [ m ]');
   ytit = sprintf('wavefunction \n y(z,t)  [a.u.]');
   ylabel(ytit);
   
   xP = z(1:NR:end); yP = s(1:NR:end,c);
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',6,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
  
   hPlot = plot(xP(1), yP(1),'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[0,0,0], ...
       'MarkerFaceColor', [0,0,0]);
   
    hPlot = plot(L, yP(1),'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[0,0,0], ...
       'MarkerFaceColor', [0,0,0]);
   
   hPlot = plot(xP(2), yP(2),'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor','m', ...
       'MarkerFaceColor', 'm');
   
   hPlot = plot(xP(2)+L, yP(2),'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor','m', ...
       'MarkerFaceColor', 'm');
   
   axis([0 zMax+1 -11 11]);
   hold(ax1, 'off');
   grid on;
   box on;
   set(gca,'xTick',0:5:zMax);
   set(gca,'fontsize',12);
   
% BLACK PARTICLE -------------------------------------------------------
   pos1 = [0.15 0.35 0.7 0.125];
   ax2 = subplot('Position',pos1);   
   hold(ax2, 'on')
   title(blacktitle,'fontweight','normal');
   xlabel('time  t  [ s ]');
   ytit = sprintf('wavefunction \n y(z_B, t)  [a.u.]');
   ylabel(ytit);
   xT1(c) = t(c); 
   yT1(1:c) = s(1,1:c);
   plot(xT1(1:c), yT1(1:c),'color',[0,0,0],'LineWidth',3);
   axis([0 tMax -10 10]);
   set(gca,'xTick',0:10:tMax);
   set(gca,'fontsize',12);
   box on;
   grid on
   hold(ax2, 'off');
   
   
 %  CYAN PARTICLE -------------------------------------------------------
   pos1 = [0.15 0.1 0.7 0.125];
   ax3 = subplot('Position',pos1);
   hold(ax3, 'on')
   title(cyantitle,'fontweight','normal');
   xlabel('time  t  [ s ]');
   ytit = sprintf('wavefunction \n y(z_C, t)  [a.u.]');
   ylabel(ytit);
   xT2(c) = t(c); 
   yT2(1:c) = s(100,1:c);
   plot(xT2(1:c), yT2(1:c),'color','m','LineWidth',3); 
   axis([0 tMax -10 10]);
   set(gca,'xTick',0:10:tMax);
   set(gca,'fontsize',12);
   box on;
   grid on
   hold(ax3, 'off')
   pause(0.1)
   
% Adding the frame to the video object using the 'writeVideo' method
   if flagA == 2
     frame = getframe(f1);
     writeVideo(aviObj,frame);
   end

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

if flagA == 2
  close(aviObj);
end