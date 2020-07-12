% wm_refraction_01.m

% REFRACTION  [2D]
% Animation of the refraction of a travelling wave

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% 170529

clear all
close all
clc
tic

% Setup for saving images (im)  =======================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_wm_refraction0.gif';
% Delay in seconds before displaying the next image  
   delay = 0.25;  
% Frame counter start
   nt = 1; 

% CALCULATIONS   ====================================================== 

N = 200;

wL1 = 20; wL2 = 20;
k1 = 2*pi/wL1; k2 = 2*pi/wL2;
v = 25;
T = wL1/v;
tmax = 4*T;

% interface
x1 = 30; y1 = 0; x2 = 80; y2 = 100;


t0 = atand((y2-y1)/(x2-x1));
t1 = 90 - t0;
t2 = asind((wL2/wL1)*sind(t1));

m1 = tand(t0);
b1 = y1 - m1 * x1;
x = linspace(0,100,N);
y = x;

z = zeros(N,N);

% GRAPHICS ===============================================================

figure(1)
   pos = [0.1 0.1 0.4 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
  % axis off
   set(gca,'fontsize',12);
   t1 = ' t / T =  ';
   
   
for t = 0 : 0.01 : tmax
  for cx = 1 : N
    if cx < N/2
        
    z(:,cx)= sin(k1 * x(cx)-2*pi*t/T);
    else
   % z(:,cx)= sin(k1 * x(cx)-2*pi*t/T);   % for no refraction - no change in speed  
     z(:,cx)= sin(k1 * x(N/2) + k2 * (x(cx)-x(N/2)) - 2*pi*t/T);    
    end
  end
  set(gca,'ytick',[]);
  pcolor(z)
  set(gca,'ytick',[]);
  % axis off
  % axis square
   shading interp %flat
   colormap(winter)
   %drawnow
   t2 = num2str(t/T,'%2.1f \n');
   tm = [t1 t2];
   set(gca,'fontsize',12);
   title(tm);
   pause(0.01)
   set(gca,'ytick',[]);
  % set(gca,'ycolor',[1 1 1]);
   set(gca,'fontsize',12)
   set(gca,'xtick',0:20:200)
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
   clf
  % set(gca,'ytick',[]);
end

toc
