% cemEDR03.m

% ELECTRIC DIPOLE RADIATION
%   GENERATION OF ELECTROMAGNETIC WAVES
%   Dipole in Z direction

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Notes
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemEDR.htm
% SCRIPTS
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb



clear
close all
clc

% Animated gif  =======================================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flag2
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_cemEDR03.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0;
    frame1 = 0;

tic
% Constant C
  CMin = 0.25; CMax = 1.35; dC = 0.30;
  C = CMin:dC:CMax; CLen = length(C);
% time
  wtMin =  0; wtMax = 2*pi; dwt = pi/64;
  wt = wtMin:dwt:wtMax; wtLen = length(wt);
% Displacement
  krMin = 0.05; krMax = 12; dkr = 0.01;
  kr = krMin:dkr:krMax; fr = krMax/11; krLen = length(kr);

  x = zeros(krLen,CLen);
  y = zeros(krLen,CLen);

  figure(1)
  set(gcf,'color',[1 1 1])

  for m = 1 :wtLen
      hold off
  for c = 1:CLen
    s2 = abs(C(c)./( (1./kr).*(sin(kr-wt(m))) - cos(kr - wt(m) )  ));

 for n = 1:length(s2)
 x(:,c) = real(fr.*kr.*sqrt(s2));
 y(:,c) = real(fr.*kr.*sqrt(1-s2));
 end
 end
 for c = 1:CLen
 plot(x(:,c),y(:,c),'b-','linewidth',1)
 hold on
 plot(x(:,c), -y(:,c),'b-','linewidth',1);
 plot(-x(:,c),-y(:,c),'b-','linewidth',1);
 plot(-x(:,c), y(:,c),'b-','linewidth',1);
 plot([-12 12], [0 0],'color',[1 1 1],'linewidth',1.0)
 xlim([-12 12])
 ylim([-12 12])
 axis square
 axis off
 
 end
  pause(0.001)

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
  toc