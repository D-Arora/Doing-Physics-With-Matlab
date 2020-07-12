% oscC001.m
% 190517  Matlab 2018b

% SIMULATION OF THE MOTION OF SINGLE OSCILLATOR
% All parameters are changed within the Script
% [default values and units]

% DOING PHYSICS ONLINE 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/oscC001A.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney


clear
close all
clc

% Animated gif  =======================================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flag2
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_osc.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;
    
% =====================================================================    
% Particles and Springs  [m = 1 kg   k = 1 N/m]
   N = 1;
   m0 = 1;
   m = m0 .* ones(N,1); 
   k0 = 1;
   k = k0 .* ones(N+1,1);
   
% Time / Time Steps [5001] [0 s  50 s]
   Nt = 1001;   % must be an ODD number
   tMin = 0;
   tMax = 50;
   
   t = linspace(0,tMax,Nt);
   dt = t(2) - t(1);   

% Model Parameters
%   Natutal frequency  [Hz]   /   natural period  [s]
    f0 = sqrt((k(1)+k(2))/m(1)) / (2*pi);
    T0 = 1/f0;
    
%   Damping constant   [b = 0.1  kg.s^-1]
    b = 0.1;
    
%   Driving force: period [TD = 1.6 T) s] / amplitude [AD = 0.005 N]
    AD = 0.000;
    TD = 1*T0;
    F = zeros(N,Nt);
    % Amplitude of driving force  [ A = 0.1  N]
   
    F(1,:) =  AD*cos(2*pi*t/TD); 
   
% Displacement  [d = 0.1 m]
%    Equilibrium separation between particles [m]
     d = 0.1;
% Displacement of each particle [m]
     e = zeros(N,Nt);
% Position of each particle
     x = zeros(N,Nt);
     
% Initial conditions   [amplitude  e(1,1) = 0.002 m]
  e(1,1) = 0.0025;
  e(1,2) = 0.0025;
 % e(1,2) = e(1,1) + 0.5*(k(1)+k(2))*e(1,1)*dt^2;
  
% Displacement from equilibrium [m]   FINITE DIFFERENCE METHOD
  
for nt = 1 : Nt-2
    e(1,nt+2) = 2*e(1,nt+1) - e(1,nt)  ...
        -(k(1)*dt^2/m(1)) * e(1,nt+1) ...
        -(k(2)*dt^2/m(1)) * e(1,nt+1) ... 
        -(dt*b/m(1) )* (e(1,nt+1) - e(1,nt)) ...
        + (dt^2/m(1)) * F(1,nt+1); 
end

% Velocity  [m/s]
   v = zeros(1,Nt);
   v(1,:) = gradient(e(1,:),dt);
   v_N = max(v(1,:));   % max velocity - numerical
   v_A = 2*pi*f0*e(1,1);   % max velcoity - analytical
   
% Acceleration  [m/s^2]
  a = zeros(1,Nt); 
  a(1,:) = ( -k(1)*e(1,:) - k(2)*e(1,:) - b*v(1,:) + F(1,:) ) ./ m(1);
  a_N = max(a(1,:));       % max acc - numerical
  a_A = (2*pi*f0)^2*e(1,1);   % max acc - analytical
  
% Fourier Transform
  fMin = 0.1*f0;
  fMax= 2*f0;
  Nf = 501;
  H = zeros(1,Nf);
  f = linspace(fMin,fMax,Nf);
  
  for c = 1:Nf
    g = e(1,:) .*  exp(1i*2*pi*f(c)*t);
    H(c) = simpson1d(g,tMin,tMax);
  end
  
  PSD = 2.*conj(H).*H;
  PSD = PSD./max(PSD);
 
  [xx, yy] = findpeaks(PSD,f,'MinPeakProminence',0.5);
  
  f_peak = yy;
  
  
% OUPUTS ==============================================================
disp('SHM:  N numerical  / A analytical  ')
disp('Natural frequency ')
fprintf('   f_A = %3.4f  Hz   \n',f0)
fprintf('   f_N = %3.4f  Hz   \n',f_peak)
disp('Natural period   ')
fprintf('   T_A = %3.4f  s   \n',T0)
fprintf('   T_N = %3.4f  s   \n',1/f_peak(1))
disp('Max displacement from equilibrium   ')
fprintf('   e_max = %3.4f  m   \n',max(e(1,:)))
disp('Max velocity   ')
fprintf('  v_A = %3.4f  m/s   \n',v_A)
fprintf('  v_N = %3.4f  m/s   \n',v_N)
disp('Max acceleration   ')
fprintf('  a_A = %3.4f  m/s^2   \n',a_A)
fprintf('  a_N = %3.4f  m/s^2   \n',a_N)

% GRAPHICS ============================================================
  figure(1)
    pos = [0.05 0.1 0.3 0.7];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
  
  subplot(4,1,1)
    xP = t; yP = d + e(1,:);
    plot(xP,yP,'b','linewidth',2)
    xlabel('t  [ s ]')
    ylabel('x  [ m ]')
    grid on
    set(gca,'fontsize',12)

  subplot(4,1,2)
    xP = t; yP = v(1,:);
    plot(xP,yP,'b','linewidth',2)
    xlabel('t  [ s ]')
    ylabel('v  [ m/s ]')
    grid on
    set(gca,'fontsize',12)   
  
  subplot(4,1,3)
    xP = t; yP = a(1,:);
    plot(xP,yP,'b','linewidth',2)
    hold on
    xlabel('t  [ s ]')
    ylabel('a  [ m/s^2 ]')
    grid on
    set(gca,'fontsize',12)  
    
 subplot(4,1,4)
    xP = e(1,:); yP = v(1,:);
    plot(xP,yP,'b','linewidth',2)
    ylabel('v  [ m/s ]')
    xlabel('x  [ m ]')
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

if flagAG == 1    
     % May have to adjust limits for different simulations
      Xmax = 6; Xtick = -Xmax:1:Xmax; 
      if max(F) > 0
        Xmax = 100; Xtick = -Xmax:20:Xmax;  
      end
      
figure(9)
    pos = [0.4 0.5 0.3 0.1];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
    
    for c = 1: 5 : Nt
       xP = 1e3*e(1,c); yP = 0;
       hPlot = plot(xP,yP,'ob','linewidth',2);
       set(hPlot,'markersize',12,'markerfacecolor','b')
       xlabel('e  [mm]')
       xlim([-Xmax Xmax])
       set(gca,'xtick',Xtick)
       ylim([-1 1])
       grid on
       set(gca,'ytick',[])
       set(gca,'fontsize',12) 
       pause(0.1)
       hold off
        
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
    