%  Waves - interference - animation
%  Ian Cooper
%  School of Physics, University of Sydney

clear
close all
clc
tic

nfile = 'interference.m';

% INPUTS --------------------------------------------------------
A = 1;                    % wave amplitude
wL = 4e-3;                   % wavelength
v = 10;                   % speed of wave

xmin = 5e-3;                 % detector dimensions
xmax = 4e-2;
ymin = -3e-2;
ymax = 3e-2;
N = 1000;
Nt = 13*3;
xs = 0;                   % source (slit) positions
ys(1) = 0.6e-2;
ys(2) = -0.6e-2;

% SETUP ---------------------------------------------------------
f = v / wL;              % frequency
w = 2*pi*f;              % angulat frequency
T = 1/f;
k = 2*pi/wL;             % wave number
dt = T/12;
t = 0;
x = linspace(xmin,xmax,N);   % detector positions
y = linspace(ymin,ymax,N);

[xx yy] = meshgrid(x,y);
r1 = sqrt((xx-xs).^2 + (yy-ys(1)).^2);
r2 = sqrt((xx-xs).^2 + (yy-ys(2)).^2);

% CALCULATIONS --------------------------------------------------

psi1 = (A./r1) .* exp((1i*k).*r1);
psi2 = (A./r2) .* exp((1i*k).*r2);

psi = psi1 + psi2;

% figure(1)
% contourf(xx,yy,real(psi));
% shading flat
% colorbar
% axis equal
% axis([xmin xmax ymin ymax])
% colormap jet

figure(2)
set(gcf,'color',[1 1 1]);

contourf(xx,yy,real(psi));
shading flat
colormap(jet)
axis equal
axis off
%caxis([-190, +190]);

M = getframe;
[im,map] = rgb2ind(M.cdata,512,'nodither');  %RGB to indexed images
im(1,1,1,10) = 0;

for c = 1 : Nt
t = t+dt;   
psit = psi .* exp(-1i*w*t);
psit(1,end) = max(max(real(psi)));
psit(1,1) = -max(max(real(psi)));

contourf(xx,yy,real(psit));
%caxis([-200, +200]);
shading flat
axis equal
axis off
%axis([xmin xmax ymin ymax]);
  

drawnow
M = getframe;
im(:,:,1,c) = rgb2ind(M.cdata,map,'nodither');

end

toc


%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously
ag_name = 'ag_interference.gif';
delay = 0;
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

