% op_002.m

% Animation of a [2D] plane wave

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1001.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181103

close all
clear 
clc
tic

% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 0;
%  Enter file name
     ag_name = 'agOP2.gif';
%  Delay in seconds before displaying the next image  
    delay = 0;  
%  Frame counter start
    nt = 1;

% SETUP ===============================================================
   wL = 20;     % wavelength

   nX = 100;    % spatial grid points
   nT = 50;     % number of time steps
   T = 100;     % period

% CALCULATIONS ========================================================
   x = 1:nX;            % X axis 
   t = linspace(0,T,nT);   % time
   k = 2*pi/wL;            % propagation constant
   w = 2*pi/T;             % angular frequency

   USS = zeros(nX,nX);     % [2D] spatial wavefunction

   US = exp(1i.*k.*x);     % [D] spatial wave function

for cc = 1:nX
   USS(cc,:) = US;
end

% GRAPHICS ============================================================

figure(3)
  FS = 12;
  pos = [0.6 0.05 0.35 0.40];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  axis off
  
for cc = 1 : nT
  xy = real(USS.*exp(-1i.*w.*t(cc)));  
  pcolor(xy)
  shading interp
  axis off
  colormap(winter(1024))
  
  if flagS == 1
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
  
  pause(0.1)
end

toc
