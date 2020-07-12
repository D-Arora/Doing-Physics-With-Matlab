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
    flagAG = 0;    
% file name for animated gif   
    ag_name = 'ag_osc.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;

    
% SETUP =============================================================== 

% Particles and Springs  [m = 1 kg   k = 1 N/m] -----------------------
   N = 21;
   m0 = 0.1;
   m = m0 .* ones(N,1); 
   k0 = 10;
   k = k0 .* ones(N-1,1);
      
% Time / Time Steps [Nt = 1001] [tMin = 0 s  tMax = 50 s] -------------
   Nt = 1001;   % must be an ODD number
   tMin = 0;
   tMax = 60;
   t = linspace(0,tMax,Nt);
   dt = t(2) - t(1);   

% Particle 2: Natural frequency [Hz] / natural period [s] -------------
%     f0 = sqrt((k(1)+k(2))/m(2)) / (2*pi);
%     T0 = 1/f0;
    
% Damping constant   [b = 0.1  kg.s^-1] ------------------------------
    b = 0.0;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
% Driving force: frequency [fD0 = 0.2 Hz] / amplitude [AD = 0.0007 N]
    fD0 = 0.2503;    
    AD = 0.0003;
% Particle to be excited by driving force   (comment / uncomment) 
    nE = ceil(N/2);   % odd harmonics 
%   nE = ceil(N/4);   % even harmonics
    FD = zeros(N,Nt);
    FD(nE,:) =  AD*cos(2*pi*fD0*t); 

% Displacements --------------------------------------------------------    
% Equilibrium separation between particles [d = 0.1 m]
     d = 0.010;
% Displacement of each particle [m]
     e = zeros(N,Nt);
% Position of each particle
     x = zeros(N,Nt);
     
% Initial conditions [0.0025 m]
%   e(11,1) = 0.00;
%   e(11,2) = 0.00;
 %  e(3,1) = -0.0025;
 %  e(3,2) = -0.0025;

 
% CALCULATIONS ======================================================== 
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

  fMin = 0.01*fD0;
  fMax = 2*fD0;
  Nf = 901;
  H = zeros(1,Nf);
  f = linspace(fMin,fMax,Nf);
  
  for c = 1:Nf
    g = e(nE,:) .*  exp(1i*2*pi*f(c)*t);
    H(c) = simpson1d(g,tMin,tMax);
  end
  
  PSD = 2.*conj(H).*H;
  PSD = PSD./max(PSD);
 
  [xx, yy] = findpeaks(PSD,f,'MinPeakProminence',0.2);
  
  f_peak = yy;
  
  
% Natural Frequencies ================================================
  modeN = 5;
  omega = zeros(modeN,1);
  for c = 1 : modeN
      omega(c) = sqrt( (4.*k0/m0) * sin(c*pi/(2*(N+1)))^2 );
  end
  f_mode = omega./(2*pi)
  
%p = -(k(1) + k(2))/m(2) - (k(2) + k(3))/m(3) ;%
%    q = (k(1)*k(2) + k(2)*k(3) + k(1)*k(3)) / (m(2)*m(3));
%    omega(2) =  sqrt((-p + sqrt(p^2 - 4*1*q))/2);
%    omega(1) =  sqrt((-p - sqrt(p^2 - 4*1*q))/2);
%    f_natural = omega ./ (2*pi);
%    ratio = -k(2) ./ ( m(1).*omega.^2 - k(1) - k(2) );
%    A_ratio = -k(2) ./ (m(2).*omega - k(1) - k(2));
%    
  
%
eMax = zeros(N,1);
eMin = zeros(N,1);
for c = 1 : N
   eMax(c) = max(e(c,:));
   eMin(c) = min(e(c,:));
end   
   
% OUPUTS ==============================================================

%   disp('FREE MOTION: Natural modes of vibration  ')
%   disp('Analytical A ')
%   fprintf('   f_A = %3.4f  Hz   \n',f_natural(1))
%   fprintf('   f_A = %3.4f  Hz   \n',f_natural(2))
%   disp('Numerical N  (Fourier Transform peaks)  ')
%  % fprintf('   f_N = %3.4f  Hz   \n',f_peak(1))
%   if length(f_peak) > 1
%     fprintf('   f_N = %3.4f  Hz   \n',f_peak(2))
%   end


% GRAPHICS ============================================================

% ANIMATION ===========================================================    
if flagAG == 1  
    % May have to adjust limits for different simulations
      Xmax = 1e3*(N+1)*d; Xtick = 0:25:Xmax;     
    %  Xmax = 50; Xtick = 0:10:Xmax; 

figure(9)
    pos = [0.5 0.2 0.3 0.4];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
    xlim([0 Xmax])
    set(gca,'xtick',Xtick)

for c = 1: 5 : Nt
    
 subplot('position',[0.1 0.8 0.8 0.1])   
    for cm = 1 : N
       xP = 1e3.*(cm*d + e(cm,c)); yP = 0;
       hPlot = plot(xP,yP,'ob');
       set(hPlot,'markersize',8,'markerfacecolor','b')
       hold on
    end
      xlabel('x  [mm]')
      xlim([0 Xmax])
      set(gca,'xtick',Xtick)
      ylim([-1 1])
      grid on
      set(gca,'ytick',[])
      set(gca,'fontsize',12) 
      pause(0.01)
      hold off
      xlim([0 Xmax])
      set(gca,'xtick',Xtick)
      
 subplot('position',[0.1 0.2 0.8 0.4])  
  %subplot(2,1,2)   
    for cm = 1 : N
       xP = 1e3.*(cm*d + e(cm,c)); yP = 1e3.*e(cm,c);
       hPlot = plot(xP,yP,'ob');
       set(hPlot,'markersize',8,'markerfacecolor','b')
       hold on
    end
      xlabel('x  [mm]')
      xlim([0 Xmax])
      set(gca,'xtick',Xtick)
      ylim([-1e3*max(eMax) 1e3*max(eMax)])
      grid on
      set(gca,'ytick',[])
      set(gca,'fontsize',12) 
      pause(0.01)
      hold off
      xlim([0 Xmax])
      set(gca,'xtick',Xtick)
  
      
%        frame1 = frame1 + 1;
%          frame = getframe(9);
%          im = frame2im(frame);
%          [imind,cm] = rgb2ind(im,256);
%       % On the first loop, create the file. In subsequent loops, append.
%          if frame1 == 1
%            imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
%          else
%           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
%          end
    
end

end

figure(1)
    pos = [0.05 0.1 0.3 0.5];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
  
  subplot(2,1,1)
    xP = t; yP =  1e3.*e(nE,:);
    plot(xP,yP,'r','linewidth',2)
    xlabel('t  [ s ]')
    ylabel('e_D  [ mm ]')
    grid on
    set(gca,'fontsize',12)

 subplot(2,1,2)
    xP = 0:d:(N-1)*d; yP =  1e3.*eMax;
    plot(xP,yP,'b','linewidth',2)
    hold on
    yP =  1e3.*eMin;
    plot(xP,yP,'b','linewidth',2)
   
    yP = zeros(1,length(xP));
    hPlot = plot(xP,yP,'ok');
    set(hPlot,'markersize',4,'markerfacecolor','k')
    
    xP = (nE-1)*d; yP = 0;
    hPlot = plot(xP,yP,'or');
    set(hPlot,'markersize',5,'markerfacecolor','r')
    
    xlabel('x  [ mm ]')
    ylabel('e  [ mm ]')
    ylim([-1.2e3*max(eMax) 1.2e3*max(eMax)])
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
  %  title(tm)
    

 