% refraction_02.m

% REFRACTION  [2D]
% Animation of the refraction of a travelling wave

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% 170529
clear all
close all
clc
tic

% Setup for saving images (im)  =======================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_wm_refraction2.gif';
% Delay in seconds before displaying the next image  
   delay = 0.25;  
% Frame counter start
   nt = 1; 

% CALCULATIONS   ====================================================== 

N = 500;

wL1 = 10; wL2 = 5;
k1 = 2*pi/wL1; k2 = 2*pi/wL2;
v = 20;
T = wL1/v;
tmax = 2*T;

% INTERFACES 
   x1 = 30; y1 = 0; x2 = 80; y2 = 100;    % angle of incidence > 0
   

t0 = atand((y2-y1)/(x2-x1));
t1 = 90 - t0;
t2 = asind((wL2/wL1)*sind(t1));
p2 = t1 - t2;

m1 = tand(t0);
b1 = y1 - m1 * x1;

x = linspace(0,100,N);
y = x;

m2 = -tand(p2);

z = zeros(N,N);


% GRAPHICS ===============================================================

figure(1)
   pos = [0.1 0.1 0.25 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   axis off
   set(gca,'fontsize',12)
   t1 = ' t =  ';
   t3 = '  s';
   
for t = 0 : 0.02 : tmax
 for cx = 1 : N
    for cy = 1 : N
    z(cy,cx)= sin(k1 * x(cx)-2*pi*t/T);
        if x(cx) > x1 && y(cy) < m1 * x(cx) + b1
        b2 = y(cy) - m2 * x(cx);
        xc = (b2 - b1)/(m1 - m2);
        yc = m1 * xc + b1;
        dx = x(cx) - xc;
        d1 = xc;
        d2 = dx / cosd(p2);
        z(cy,cx)= sin(k1 * d1 + k2 * d2 - 2*pi*t/T);    
        end
    end
 end   
   pcolor(z)
   axis off
   axis square
   shading flat
   colormap(winter)
   %drawnow
   t2 = num2str(t,'%2.2f \n');
   tm = [t1 t2 t3];
   title(tm);
   set(gca,'fontsize',12)
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

 


toc
