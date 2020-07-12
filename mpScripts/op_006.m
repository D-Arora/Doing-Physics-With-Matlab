% op_006.m

% POLARIZATION: MALUS' LAW
% Animation of light tranmsitted through a polarizer and an analyzer
%   as the analyzer is rotated through 180 deg

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1007.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181111

clc
close all
clear


% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 1;
%  Enter file name
     ag_name = 'ag_A.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.25;  
%  Frame counter start
    nt = 1;

% =====================================================================
% Number of angle through which analyzer is rotated
   N = 181;
% Angle between transmission axes of polarizer and analyzer
  theta = linspace(0,180,N);
% Color coding of light through optical elements  
   c = abs(cosd(theta)).^2;
% Irradiance of light [W.m^-2] 
   S0 = 100;
   S = S0.*cosd(theta).^2;

% GRAPHICS ============================================================
figure(1)
  pos = [0.05 0.05 0.22 0.30];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
    
  
for cc = 1 : N
  
  rectangle('Position',[1,2,5,5],...
    'Curvature',[1,1], 'FaceColor',[c(cc) 0 0])


  tm1 = '\theta  =  ';
  tm2 = num2str(theta(cc), '%2.0f  deg  \n');
  tm3 = '    S  =  ';
  tm4 = num2str(S(cc), '%2.0f  W.m^{-2}  \n');
  tm = [tm1 tm2 tm3 tm4];
  title(tm,'fontsize',14)  

  axis square;
  axis off
  
  if flagS == 1
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
  
  pause(0.1)

end

figure(2)
  pos = [0.30 0.05 0.22 0.30];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = theta; yP = S;
  plot(xP,yP,'r','linewidth',2)
  
  grid on
  xlim([0 180])
  set(gca,'xtick',0:30:180)
  
  xlabel('\theta [ deg ]')
  ylabel('Irradiance  S  [ W.m^{-2} ]')
  
  set(gca,'fontsize',12)
  
  