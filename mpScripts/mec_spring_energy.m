% mec_spring_energy.m
% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

clear all
close all
clc

xE = [80 80];
yE = [20 80];
xK = [88 88];
yK = yE;
xU = [96 96];
yU = yE;

%xS = [0 5 15 25 35 45 55 60];
xp = linspace(0,60,8);
dx = xp(2)-xp(1);
xS(2) = xp(1)+dx/2;
xS(3) = xp(2)+dx; xS(4) = xp(3)+dx; xS(5) = xp(4)+dx;
xS(5) = xp(4)+dx; xS(6) = xp(5)+dx; xS(7) = xp(6)+dx;
xS(8) = xp(8);

y1 = 40; y2 = 60;
yS = [50 y1 y2 y1 y2 y1 y2 40];
xM = 60; yM = 50;
xS(end) = 60; yS(end) = 50;

ms = 30;
cc = 0;

A = 8;
T = 100;

figure(1)
set(gcf,'color',[1 1 1]);
set(gcf,'Position',[200 200 500 300]);
% plot(xE,yE,'b','linewidth',10);
% hold on
% plot(xK,yK,'r','linewidth',10);
% plot(xU,yU,'m','linewidth',10);
% axis([0 100 0 100])
% plot(xS,yS,'k','linewidth',3)
% plot(xM,yM,'o','markersize',ms,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold off
M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,256,'nodither');  %RGB to indexed images
im(1,1,1,1) = 0;


for t = 0:100
x = -A*cos(2*pi*t/T);
xp = linspace(0,60+x,8);
xS(1) = 0;
xS(2) = xp(1)+dx/2;
xS(3) = xp(2)+dx; xS(4) = xp(3)+dx; xS(5) = xp(4)+dx;
xS(5) = xp(4)+dx; xS(6) = xp(5)+dx; xS(7) = xp(6)+dx;
xS(8) = xp(8);

xs = xS + x;
xs(1) = 0; xM = xs(end);

xs(xs<0) = 0;

yU(2) = yU(1) + 60*(x/A)^2; 
yK(2) = yK(1) + 60*(1-(x/A)^2);

plot(xE,yE,'b','linewidth',10);
hold on
plot(xK,yK,'r','linewidth',10);
plot(xU,yU,'m','linewidth',10);
axis([0 100 0 80])
plot(xs,yS,'k','linewidth',3)
plot(xM,yM,'o','markersize',ms,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot([60 60],[30 70],'k');
text(55.5,27,'x = 0 ','fontsize',14);
text(78,15,'E','fontsize',14);
text(86,14,'E_K','fontsize',14);
text(95,14,'E_P','fontsize',14);
hold off
axis off
pause(0.1)

M = getframe(gcf) ;
cc = cc+1;
im(:,:,1,cc) = rgb2ind(M.cdata,map,'nodither');
end


delay = 0.2;
ag_name = 'ag_mec_spring_energy.gif';   % file name for animated gif
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);
