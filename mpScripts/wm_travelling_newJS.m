% osc_shm_01A.m

% SIMPLE HARMONIC MOTION
% Simulation for the vertical motion of an object attached to a spring
% Animation of the motion can be saved as an animated gif

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170517

%Modified by John A Sims, 
% john.sims@ufabc.edu.br
%Universidade Federal do ABC, 
%http://www.ufabc.edu.br/
%Biomedical Engineering department,
%Centro de Engenharia, Modelagem e CiÃªncias Sociais Aplicadas.Â 
% (Centre for Engineering, modelling and applied social sciences)
%Rua Arcturus (Jd Antares)
%Anchieta
%09606070 - SÃ£o Bernardo do Campo, SP - Brasil
%

clear all
close all
clc


flagC = 2;
% flagC = 1   animation of a trvelling wave
% flagC = 2   animation of a transverse wave
% flagC = 3   animation of a longitudinal wave

% =======================================================================
%    Setup for saving images (im) 
% ======================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_wm_trav.gif';
% Delay in seconds before displaying the next image  
   delay = 0.20;  
% Frame counter start
   nt = 1; 

switch flagC
    
case 1   % ------------------------------------------------------------
        
% amplitude A  / wavelength L / period T
% no. periods nT / no. wavelengths steps nX
% Number of calcutions Nt Nx
   A = 10;
   L = 20;
   T = 25;
   nT = 3; nX = 3;
   Nx = 200; 
   Nt = 200;
   
% left (-1) or right (+1)
   e = 1;
% angular frequency w / time t / wave number k /displacement s
   w = 2*pi / T;
   k = 2*pi / L;
   tMax = nT * T;
   xMax = nX * L;
   
   t = linspace(0,tMax,Nt);
   x = linspace(0,xMax,Nx);
   
% wave function s(x,t)   
   s = zeros(Nx,Nt);
   
   for c = 1 : Nt
     s(:,c) = A .* sin(k.*x + e*w*t(c));
   end
     v = L / T;
     
    
figure(1)
   pos = [0.1 0.1 0.4 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   t1 = 't  =  ';
   t3 = '   s';
   
for c = 1 : Nt
   hold off
   xP = x; yP = s(:,c);
   plot(xP, yP,'b','lineWidth',2);
   hold on
   xP = x(100); yP = s(100,c);
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
   xP = v * t(c); yP = 0;
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',6,'MarkerEdgeColor',[1 0 1], ...
       'MarkerFaceColor', [1 0 1]);
   
   axis([0 xMax -11 11]);
   ylabel('wave function  s ');
   xlabel('position  x  [ m ]')
   set(gca,'fontSize',14);
   grid on
   t2 = num2str(t(c), '%2.1f \n');
   tm = [t1 t2 t3];
   title(tm);
   pause(0.2)
   
   if f_gif > 0 
       frame = getframe(1);
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

%transverse wave.
case 2  % --------------------------------------------------------------
   % amplitude A  / wavelength L / period T
% no. periods nT / no. wavelengths steps nX
% Number of calcutions Nt Nx

   A = 10;
   L = 20;
   T = 25;
   nT = 3; nX = 3;
   Nx = 2000; 
   Nt = 200;
   
% left (-1) or right (+1)
   e = -1;
% angular frequency w / time t / wave number k /displacement s
   w = 2*pi / T;
   k = 2*pi / L;
   tMax = nT * T;
   xMax = nX * L;
   
   t = linspace(0,tMax,Nt);
   x = linspace(0,xMax,Nx);
   
   % wave function s(x,t)   
   s = zeros(Nx,Nt);
   % calculate all data and store it in S, rows are x, cols t. 
   
   % JSims 2017
   % Later, phase relationship between all particles in an instant of time
   % (i.e spatial domain information) will be taken from rows of S, 
   % and time (temporal domain) for a single particle taken from columns.
   
   for c = 1 : Nt
     s(:,c) = A .* sin(k.*x + e*w*t(c));
   end
    
    f2 = figure(2) ;
   pos = [0.1 0.1 0.4 0.7];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   NS = 20; 
   yT1 = zeros(Nt,1);   xT1 = zeros(Nt,1);
   yT2 = zeros(Nt,1);   xT2 = zeros(Nt,1);
   
   %added by J Sims
   redFreq=100; % adjust to reduce/increase the frequency of red circles. 
            % compared to blue circles.  I.e. Every "redFreq" is to be red.
         
   %define text for three figures based on user choice of e.  
   %Figure 1 - spatial domain - main title.
   %Figure 2 - temporal domain - blacktitle (it is the black dot)
   %Figure 3 - temporal domain - cyantitle (guess why).
   blacktitle= sprintf('\nAmplitude of black particle through time.');    
   cyantitle = sprintf('Amplitude of cyan particle through time.\n');
   
   switch e
       case 1
   eqnname = 'y(z,t) = 10sin(kz + \omega t)';
   maintitle = sprintf('Phase relation of particles in space (velocity -z direction)');
   blacktitle = sprintf('%s\ny(t)= 10sin(\\phi_b + \\omega t) ',blacktitle);
   cyantitle = sprintf('%s y(t) 10sin(\\phi_c + \\omega t) ',cyantitle); 
   aviName = 'transWave_RtoL.avi';
       case -1
   eqnname = 'y(z,t) = 10sin(kz - \omega t)';
   maintitle = sprintf('Phase relation of particles in space (velocity +z direction)');
   blacktitle = sprintf('%s\ny(t)= 10sin(\\phi_b - \\omega t) ',blacktitle);
   cyantitle = sprintf('%s y(t) 10sin(\\phi_c - \\omega t) ',cyantitle); 
   aviName = 'transWave_LtoR.avi';
       otherwise %error condition
   error('Error: Value of e must be either 1 or -1!')   
   end
   
   maintitle = sprintf('%s\n %s\n \\fontsize{6}{john.sims@ufabc.edu.br 2017, modified from code by ian.cooper@sydney.edu.au}',maintitle,eqnname);
   blacktitle = sprintf('%s where \\phi_b = kx_b = initial phase for black = 0',blacktitle);
   cyantitle = sprintf('%swhere \\phi_c = kx_c = initial phase for cyan \\neq 0',cyantitle);

   %open writer object for avi.
   aviObj = VideoWriter(aviName);
   aviObj.FrameRate = 10;  % from 30
   aviObj.Quality = 40;    % from 75
   open(aviObj);
   %%%%%%%%%
   t = 0:Nt-1;%create time vector
   
for c = t+1
   hold off
   xP = x(1:NS:end); yP = s((1:NS:end),c);
   ax1 = subplot(4,1,[1,2]); 
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',6,'MarkerEdgeColor',[0 0 1]);%, ...
      % 'MarkerFaceColor', [0 0 1]); % How to present blue circles
   hold(ax1, 'on') %Equivalent of hold on when using subplots.

   title(maintitle);
   xlabel('z (distance)');
   ytit = sprintf('y(z,t)\n (signal amplitude)');
   ylabel(ytit);
   xP = x(1:redFreq:end); yP = s(1:redFreq:end,c);
   hPlot = plot(xP, yP,'o');
   set(hPlot,'MarkerSize',6,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
   hPlot = plot(xP(1), yP(1),'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[0,0,0], ...
       'MarkerFaceColor', [0,0,0]);
   hPlot = plot(xP(2), yP(2),'o');
   set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[0,0,1], ...
       'MarkerFaceColor', 'c');
   
   axis([0 xMax -11 11]);
   hold(ax1, 'off')

   ax2 = subplot(4,1,3);   
   hold(ax2, 'on')
   title(blacktitle);
   xlabel('t (time)');
   ylabel(ytit);
   xT1(c) = t(c); 
   yT1(1:c) = s(1,1:c);
   plot(xT1(1:c), yT1(1:c),'color',[0,0,0],'LineWidth',3);
   axis([0 Nt -11 11]);
   hold(ax2, 'off')
   
    ax3 = subplot(4,1,4);   
   hold(ax3, 'on')
   title(cyantitle);
   xlabel('t (time)');
   ylabel(ytit);
   xT2(c) = t(c); 
   yT2(1:c) = s(100,1:c);
   plot(xT2(1:c), yT2(1:c),'color','c','LineWidth',3); %orange 255 128 0
   axis([0 Nt -11 11]);
   hold(ax3, 'off')
   pause(0.1)
   
   frame = getframe(f2);
   % Adding the frame to the video object using the 'writeVideo' method
   writeVideo(aviObj,frame);


%  code for animated gif.  Replaced by writing of avi, I think it is more
%  didactic since student/teacher can slow or pause the video to think or
%  comment on it.

%    if f_gif > 0 
%        frame = getframe(gca);
%           % Adding the frame to the video object using the 'writeVideo' method
%         writeVideo(writerObj,frame);
        
        %rest of 'if' is about animated gif JS 11/2017
%        im = frame2im(frame);
%        [imind,cm] = rgb2ind(im,256);
%      %  On the first loop, create the file. In subsequent loops, append.
%        if nt == 1
%           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
%        else
%          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
%        end
%        nt = nt+1;
%    end
   %
end

%close avi object.
close(aviObj);

case 3   % LONGITUDINAL WAVE ---------------------------------------------

% x positions of particles
   x = 0:4:100;
   y = zeros(length(x),1);
   
% amplitude 
   A = 1.2;
% period and wavelength
   L = 20;
   k = 2*pi / L;
   T = 25;
   w = 2*pi / T;
% number of time steps
   Nt = 100;
   t = linspace(0,3*T,Nt);

   
figure(3) 
pos = [0.1 0.1 0.4 0.1];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
for c = 1:Nt
    hold off
    xC = x + A.* sin(k*x - w*t(c));
    xP = xC; yP = y;
    hPlot = plot(xP,yP,'o');
    set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 0], ...
       'MarkerFaceColor', [1 0 0]);
   hold on
   xP = xP(12); yP = y(12);
    hPlot = plot(xP,yP,'o');
    set(hPlot,'MarkerSize',10,'MarkerEdgeColor',[1 0 1], ...
       'MarkerFaceColor', [1 0 1]);
    axis([-10 110 -1 1]);
    
    axis off;
    pause(0.1);
    
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
    
end   % switch
   


     
     


   