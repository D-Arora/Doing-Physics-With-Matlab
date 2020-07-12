% wm_super01.m


% Animation of the motion can be saved as an animated gif
%   Superposition Principle and Interference

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170604

clear all
close all
clc

% =======================================================================
%    Setup for saving images (im) 
% ======================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_wm_super03.gif';
% Delay in seconds before displaying the next image  
   delay = 0.25;  
% Frame counter start
   nt = 1; 


% INPUTS ================================================================

% Speed of propagation / Period T / Wavelength L / amplitudes A
   v = 2.0;   
   L1 = 25;     L2 = 50;
   T1 = L1/v;   T2 = L2/v;
   A1 = 5;    A2 = 5;
   
% Time t and position x   
   tMax = 2*T1;
   xMax = 100;

   Nt = 200;
   Nx = 2000;

   t = linspace(0,tMax,Nt);
   x = linspace(0,xMax,Nx);

   w1 = 2*pi/T1;  w2 = 2*pi/T2;
   k1 = 2*pi/L1;  k2 = 2*pi/L2;

   s1 = zeros(Nt,Nx);   s2 = zeros(Nt,Nx);
for c = 1 : Nt
   s1(c,:) = A1 .* sin(k1.* x - w1.*t(c));
   s2(c,:)=  A2 .* sin(k2.* x + w2.*t(c));
end


figure(1)
   pos = [0.1 0.1 0.4 0.5];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gca,'fontsize',14);
   
for c = 1 : 1 :  Nt
    
xP = x;  yP = s1(c,:);  % 111111111111111111111111111111111111111111111
subplot(3,1,1)
  hold off
  plot(xP,yP,'lineWidth',2);
  grid on
  axis([0 xMax -11 11]);
  hold on
  xP = x(1500); yP = s1(c,1500);
  hPlot = plot(xP,yP,'ro');
  set(hPlot,'MarkerFaceColor','r');
  xP = [x(1500) x(1500)];  yP = [-11 11]; 
  plot(xP,yP,'k');
  t1 = 't =  ';
  t2 = num2str(t(c), '%2.1f\n');
  t3 = ' s';
  tm = [t1 t2 t3];
  title(tm);
  ylabel('s_1');
  set(gca,'fontsize',14);
  t0 = 'x = 75  m    ';
  t1 = 's_1 = ';
  t2 = num2str(s1(c,1500),'%2.1f\n');
  tm = [t0 t1 t2];
  text(60,7,tm,'fontsize',14);
  set(gca,'xTick',0:10:100);  
  set(gca,'yTick',-10:5:10);   
   
xP = x;  yP = s2(c,:);  % 22222222222222222222222222222222222222222222
subplot(3,1,2)
  hold off
  plot(xP,yP,'lineWidth',2);
  grid on
  axis([0 xMax -11 11]);
  hold on
  xP = x(1500); yP = s2(c,1500);
  hPlot = plot(xP,yP,'ro');
  set(hPlot,'MarkerFaceColor','r');
  xP = [x(1500) x(1500)];  yP = [-11 11]; 
  plot(xP,yP,'k');
  ylabel('s_2');
  set(gca,'fontsize',14);
  t0 = 'x = 75  m    ';
  t1 = 's_1 = ';
  t2 = num2str(s2(c,1500),'%2.1f\n');
  tm = [t0 t1 t2];
  text(60,7,tm,'fontsize',14);
  set(gca,'xTick',0:10:100);
  set(gca,'yTick',-10:5:10);
  
  xP = x;  yP = s1(c,:) + s2(c,:);   % 333333333333333333333333333333333
  subplot(3,1,3);
    hold off
   plot(xP,yP,'lineWidth',2);
   grid on
   axis([0 xMax -11 11]);
   hold on
   xP = x(1500); yP = s1(c,1500)+s2(c,1500);
  hPlot = plot(xP,yP,'ro');
  set(hPlot,'MarkerFaceColor','r');
   xP = [x(1500) x(1500)];  yP = [-11 11]; 
   plot(xP,yP,'k');
   xlabel('  x  [ m ]');
   ylabel('s_1 +  s_2');
   set(gca,'fontsize',14);
   t0 = 'x = 75  m    ';
   t1 = 's_1 + s_2 = ';
   t2 = num2str(s1(c,1500)+s2(c,1500),'%2.1f\n');
   tm = [t0 t1 t2];
   text(60,7,tm,'fontsize',14);
   set(gca,'xTick',0:10:100);
   set(gca,'yTick',-10:5:10);
  
   pause(0.01)

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

