% p001.m    

% PHOTONICS and WAVES  Phase and Group Velocity Animations
% Travelling sinusoidal waves in non-dispersive and dispersive media  

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 200418 / Matlab version R2020a

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/p001.htm
% Download Location
%    http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
  

clear
close all
clc

% Non dispersive medium  flagM = 1 /dispersive medium flagM =2
  flagM = 1;

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% Delay in seconds before displaying the next image  
    delay = 0.20;
    frame1 = 0;

% SETUP  ==============================================================
    


% Simulation time
  tMax = 10;
  nT = 80;
  t = linspace(0,tMax,nT);

% Z axis
  nZ = 500;
  zMax = 100;
  z = linspace(0,zMax,nZ);

% Amplitude
  A1 = 1; A2 = 1;
% Wavefunction  
  y1 = zeros(nZ,nT);
  y2 = zeros(nZ,nT);
  yR = zeros(nZ,nT);  

% NON DISPERIVE MEDIUM ================================================    
if flagM == 1
% Wavelength
  wL1 = 20;
% Phase velocity
  v1 = 10; v2 = v1;
% Increase in angular velocity
  dw = 0.8;    
% Angular velocity (omega w) / frequency / period
  f1 = v1 / wL1;  T1 = 1/f1;  w1 = 2*pi*f1;
  w2 = (1+dw)*w1; f2 = w2/(2*pi); T2 = 1/f2;
% Wavenumber
  k1 = w1/v1; k2 = w2/v1; wL2 = 2*pi/k2;
% Group Velocity
  vG = (w2-w1)/(k2-k1);
  
% file name for animated gif   
   ag_name = 'ag_photonicsND.gif'; 
end


% DISPERSIVE MEDIUM ==================================================
if flagM == 2
% Frequency
  f1 = 0.5; f2 = 2;
% Phase velocity
  v1 = 10; v2 = 8;
% Angular velocity (omega w) / Wavenumber
  w1  = 2*pi*f1;   w2  = 2*pi*f2;
  k1  = (w1/v1); k2  = (w2/v2);
  wL1 = 2*pi/k1;   wL2 = 2*pi/k1;
% Group Velocity
  vG = (w2-w1)/(k2-k1);
% file name for animated gif   
  ag_name = 'ag_photonicsD.gif';   
end


% WAVE FUNCTIONS  =====================================================
  for c = 1 : nT
      y1(:,c) = A1.*sin(k1.*z - w1*t(c));
      y2(:,c) = A2.*sin(k2.*z - w2*t(c));
      yR(:,c) = y1(:,c) + y2(:,c);
  end

  
% GRAPHICS  ===========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.3 0.65]);
  set(gcf,'color','w');
  
  for c = 1 : nT
    subplot(3,1,1)
     xP = z; yP = y1(:,c);
     plot(xP,yP,'b','linewidth',2);
     hold on
     xP = v1*t(c); yP = 0;
     Hplot = plot(xP,yP,'bo');
     set(Hplot,'markerfacecolor','b')
     grid on
     ylim([-1.2 1.2])
     ylabel('\psi_1','fontsize',14)
     txt = sprintf('t  =  %2.1f',t(c));
     title(txt)
     set(gca,'fontsize',12)
     hold off
    
     subplot(3,1,2)
     xP = z; yP = y2(:,c);
     plot(xP,yP,'k','linewidth',2);
     hold on
     xP = v2*t(c); yP = 0;
     Hplot = plot(xP,yP,'ko');
     set(Hplot,'markerfacecolor','k')
     grid on
     ylim([-1.2 1.2])
     ylabel('\psi_2','fontsize',14)
     set(gca,'fontsize',12) 
     hold off
     
    subplot(3,1,3)
     xP = z; yP = yR(:,c);
     plot(xP,yP,'r','linewidth',2);
     hold on
     xP = vG*t(c); yP = 0;
     Hplot = plot(xP,yP,'ro');
     set(Hplot,'markerfacecolor','r')
     grid on
     ylim([-2.2 2.2])
     ylabel('\psi_R','fontsize',14)
     xlabel('z')
     set(gca,'fontsize',12)
     text(85,-3,'p001.m','fontsize',8)
     hold off 
     pause(0.25);
     
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

