% refraction_2.m


clear all
close all
clc
tic
N = 500;

wL1 = 10; wL2 = 4;
k1 = 2*pi/wL1; k2 = 2*pi/wL2;
v = 20;
T = wL1/v;
tmax = 2*T;

% interface
x1 = 30; y1 = 0; x2 = 80; y2 = 100;


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

figure(1)
set(gf
set(gcf,'color',[1 1 1]);
axis off
c = 0;
for cx = 1:N
z(:,cx)= sin(k1 * x(cx)-2*pi*0/T);
end
pcolor(z)
%contourf(z)
axis off
shading flat
colormap(winter)
axis equal

M = getframe;
[im,map] = rgb2ind(M.cdata,256,'nodither');  %RGB to indexed images
im(1,1,1,10) = 0;

for t = 0 : 0.1/3 : tmax
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
%contourf(z)
axis off
axis equal
shading flat
colormap(winter)
drawnow
c = c+1;

M = getframe;
im(:,:,1,c) = rgb2ind(M.cdata,map,'nodither');

%pause(0.01)
end

 ag_name = 'refraction_3.gif';
   delay = 0;
  % imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);


toc
