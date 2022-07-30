% chaosLogisticsEq01.m
% Logistic Equation
% Bifurcation Diagram

% DOING PHYSICS WITH MATLAB
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/???
% Scripts
%   https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%   https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% IAN COOPER
% email    matlabvisualphysics@gmail.com

% Matlab 2021b     220526

clear
close all
clc

tic

% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_LogisticEq1.gif';
% Delay in seconds before displaying the next image  
   delay = 0;  
% Frame counter start
   nt = 1; 


% Initial condition   0 <= x0 <= 1
  x0 = 0.25;

% Control parameter 0 <= r <= 1
  rMin = 0.9609;
  rMax = 0.9615;
  Nr = 901;
  r = linspace(rMin,rMax,Nr);

% Number of iterations
  N = 100;

% Logistic Equation

  x = zeros(N,Nr);
  xEND = zeros(Nr,1);

  for k = 1:Nr
      x(1,k) = x0;
    for c = 1:N-1
      x(c+1,k) = 4*r(k)*x(c,k)*(1-x(c,k));
    end
      xEND(k) = x(c+1,k);
  end

figure(1)
  set(gcf,'units','normalized');
 set(gcf,'Position', [0.05 0.05 0.25 0.25])
 set(gcf,'color','w');
 FS = 14;

 for k = 1:2: Nr
   xP = 1:N; yP = x(:,k);
   plot(xP,yP,'b','linewidth',2)
   txt = sprintf('r = %2.3f   \n',r(k));
   title(txt)
   ylim([0 1.1])
   xlim([50 100])
   yticks(0:0.25:1)
   grid on
   box on
   xlabel('n')
   ylabel('x(n)')
   ytickformat('%1.2f')
   set(gca,'fontsize',FS)
   hold off   
   pause(0.01)

    if f_gif > 0 
          frame = getframe(1);
          im = frame2im(frame);
          [imind,cm] = rgb2ind(im,256);
  %    On the first loop, create the file. In subsequent loops, append
          if nt == 1
            imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
          else
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
          end
       nt = nt+1;
   end
 end


figure(2)
  set(gcf,'units','normalized');
 set(gcf,'Position', [0.45 0.05 0.25 0.25])
 set(gcf,'color','w');
 FS = 14;
 
 for k = 1:Nr
     for c = 60 : N
       xP = r(k); yP = x(c,k);
       plot(xP,yP,'b.')
       hold on
       xlabel('r')
       ylabel('x(n)')
       set(gca,'fontsize',FS)
       grid on
       box on
  %    xlim([0.86 0.9])
    end
 end

 toc
