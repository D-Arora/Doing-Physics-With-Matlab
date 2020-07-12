% oscC002.m
% 190517  Matlab 2018b

% SIMULATION OF THE MOTION OF TWO COUPLED OSCILLATORS
% All parameters are changed within the Script
% [default values and units]

% DOING PHYSICS ONLINE 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/oscC002.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney


clear
close all
clc

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_osc.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;

    
% SETUP =============================================================== 

% Particles and Springs  [m = 1 kg   k = 1 N/m] -----------------------
   N = 4;
   m0 = 1;
   m = m0 .* ones(N,1); 
   k0 = 1;
   k = k0 .* ones(N-1,1);
    
%    k(3) = 0.0;
%    k(2) = 3.8225;
%    m(3) = 0.1;
   
% Time / Time Steps [Nt = 1001] [tMin = 0 s  tMax = 50 s] -------------
   Nt = 1001;   % must be an ODD number
   tMin = 0;
   tMax = 50;
   t = linspace(0,tMax,Nt);
   dt = t(2) - t(1);   

% Particle 2: Natural frequency [Hz] / natural period [s] -------------
     f0 = sqrt((k(1)+k(2))/m(2)) / (2*pi);
     T0 = 1/f0;
    
% Damping constant   [b = 0.1  kg.s^-1] ------------------------------
    b = 0.0;
    
% Driving force: frequency [fD0 = 0.2 Hz] / amplitude [AD = 0.0007 N]
    fD0 = 0.1952;
    AD = 0.0000;
    
    FD = zeros(N,Nt);
    FD(2,:) =  AD*cos(2*pi*fD0*t); 

% Displacements --------------------------------------------------------    
% Equilibrium separation between particles [d = 0.1 m]
     d = 0.010;
% Displacement of each particle [m]
     e = zeros(N,Nt);
% Position of each particle
     x = zeros(N,Nt);
     
% Initial conditions [0.0025 m]
   e(2,1) = 0.0025;
   e(2,2) = 0.0025;
   e(3,1) = -0.0025;
   e(3,2) = -0.0025;

% Calculation constants
  KS = dt^2.*k;   Kb = b*dt;   KD = dt^2;   
     

% FINITE DIFFERENCE METHOD: displacement from equilibrium ==============
 
for nt = 1 : Nt-2
  for c = 2 : N-1
    e(c,nt+2) = 2*e(c,nt+1) - e(c,nt)...
                -(KS(c-1)/m(c)) * (e(c,nt+1) - e(c-1,nt+1)) ...
                -(KS(c)/m(c))   * (e(c,nt+1) - e(c+1,nt+1)) ...   
                - (Kb/m(c))     * (e(c,nt+1) - e(c,nt))     ...
                + (KD/m(c))     * FD(c,nt+1);
  end
end


% Fourier Transform  ==================================================

  fMin = 0.1*f0;
  fMax= 2*f0;
  Nf = 501;
  H = zeros(1,Nf);
  f = linspace(fMin,fMax,Nf);
  
  for c = 1:Nf
    g = e(2,:) .*  exp(1i*2*pi*f(c)*t);
    H(c) = simpson1d(g,tMin,tMax);
  end
  
  PSD = 2.*conj(H).*H;
  PSD = PSD./max(PSD);
 
  [xx, yy] = findpeaks(PSD,f,'MinPeakProminence',0.2);
  
  f_peak = yy;
  
  
 % Frequency equation ================================================
   p = -(k(1) + k(2))/m(2) - (k(2) + k(3))/m(3) ;
   q = (k(1)*k(2) + k(2)*k(3) + k(1)*k(3)) / (m(2)*m(3));
   omega(2) =  sqrt((-p + sqrt(p^2 - 4*1*q))/2);
   omega(1) =  sqrt((-p - sqrt(p^2 - 4*1*q))/2);
   f_natural = omega ./ (2*pi);
   ratio = -k(2) ./ ( m(1).*omega.^2 - k(1) - k(2) );
   A_ratio = -k(2) ./ (m(2).*omega - k(1) - k(2));
   
   
% OUPUTS ==============================================================

  disp('FREE MOTION: Natural modes of vibration  ')
  disp('Analytical A ')
  fprintf('   f_A = %3.4f  Hz   \n',f_natural(1))
  fprintf('   f_A = %3.4f  Hz   \n',f_natural(2))
  disp('Numerical N  (Fourier Transform peaks)  ')
  fprintf('   f_N = %3.4f  Hz   \n',f_peak(1))
  if length(f_peak) > 1
    fprintf('   f_N = %3.4f  Hz   \n',f_peak(2))
  end


% GRAPHICS ============================================================

figure(1)
    pos = [0.05 0.1 0.3 0.7];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
  
  subplot(3,1,1)
    xP = t; yP =  1e3.*e(2,:);
    plot(xP,yP,'b','linewidth',2)
    xlabel('t  [ s ]')
    ylabel('e_2  [ mm ]')
    grid on
    set(gca,'fontsize',12)

  subplot(3,1,2)
    xP = t; yP = 1e3.*e(3,:);
    plot(xP,yP,'b','linewidth',2)
    xlabel('t  [ s ]')
    ylabel('e_3  [ mm ]')
    grid on
    set(gca,'fontsize',12)  
    
 subplot(3,1,3)
    xP = t; yP =  1e3.*e(2,:);
    plot(xP,yP,'b','linewidth',2)
    hold on
    xP = t; yP =  1e3.*e(3,:);
    plot(xP,yP,'r','linewidth',2)
    xlabel('t  [ s ]')
    ylabel('e  [ mm ]')
    legend('e_2','e_3','orientation','horizontal')
    grid on
    set(gca,'fontsize',12)    

 figure(2)
    pos = [0.4 0.1 0.3 0.3];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on    
    xP = f; yP = PSD;
    plot(xP,yP,'b','linewidth',2)
    xlabel('f  [ Hz ]')
    ylabel('PSD [ a.u. ]')
    grid on
    set(gca,'fontsize',12) 
    
    tm = sprintf('f_{peak}  =  %3.4f   Hz \n',f_peak);
    title(tm)
    

% ANIMATION ===========================================================    
if flagAG == 1  
    % May have to adjust limits for different simulations
      Xmax = 30; Xtick = 0:5:Xmax;     
    %   Xmax = 50; Xtick = 0:10:Xmax; 

figure(9)
    pos = [0.4 0.5 0.3 0.1];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
    xlim([0 Xmax])
    set(gca,'xtick',Xtick)
    
for c = 1: 5 : Nt
       xP = 1e3.*(d+e(2,c)); yP = 0;
       hPlot = plot(xP,yP,'ob','linewidth',2);
       set(hPlot,'markersize',15,'markerfacecolor','b')
       hold on
       xP = 1e3.*(2*d+e(3,c)); yP = 0;
       hPlot = plot(xP,yP,'or','linewidth',2);
       set(hPlot,'markersize',15,'markerfacecolor','r')
       xlabel('x  [mm]')
       xlim([0 Xmax])
       set(gca,'xtick',Xtick)
    
      
      ylim([-1 1])
      grid on
      set(gca,'ytick',[])
      set(gca,'fontsize',12) 
      pause(0.1)
      hold off
      xlim([0 Xmax])
      set(gca,'xtick',Xtick)
         
       frame1 = frame1 + 1;
         frame = getframe(9);
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