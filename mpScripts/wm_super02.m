% wm_super02.m


% Animation [2D] single source wave pattern
%   

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% 170605

clear all
close all
clc

tic
% =======================================================================
%    Setup for saving images (im) 
% ======================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 1;
   ag_name = 'ag_wm_super02_02.gif';
% Delay in seconds before displaying the next image  
   delay = 0;  
% Frame counter start
   nt = 1; 

% SETUP =================================================================
Nx = 1000;
Nt = 25;

% Source positions
 % xS = 0; yS = 0;     % single source

  %xS = [-0.1,0.1];  yS = [0,0];
  xS = linspace(-0.5, 0.5, 1500);
  Ns = length(xS);
  yS = zeros(1,Ns); 
  
dx = 0.1;
x = linspace(-2.5,2.5,Nx);
y = linspace(0.5,2.5,Nx);
[xx, yy] = meshgrid(x,y);

A = 1;
f = 100;
v = 25;
T = 1/f;
L = v * T;
w = 2*pi/T;
k = 2*pi/L;

t = linspace(0,1*T,Nt);

sXY = zeros(Nx,Nx);

for c = 1 : Ns
  xxS = xx - xS(c); %xxS(xxS<0.1) = 0.1;
  yyS = yy - yS(c); %yyS(yyS<0.1) = 0.1;
  R = sqrt((xxS).^2+(yyS).^2);
  cosT = (xxS + eps)./(R+eps);
  sinT = (yyS + eps)./(R+eps);
  kx = k .* cosT; ky = k .* sinT;
  
  sXY = sXY + A .* exp(1i*(kx.*xxS + ky.*yyS))./sqrt(R);
  
end
  sXY = sXY ./ max(max(sXY));
  %sXY(sXY<-2) = -2;

% GRAPHICS ===============================================================

figure(1)
  pos = [0.1 0.2 0.35 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   

for c = 5 : Nt
 s = sXY .* exp(-1i*w*t(c));
 xP = x; yP = y; zP = real(s);
% xP(1) = x(1); yP(1) = y(1); zP(1,1) = 2;
% xP(end) = x(end); yP(end) = y(end); zP(end,end) = -2;
 pcolor(xP,yP,zP);
 %surf(xP,yP,zP);
 colormap winter(16);
 shading flat; 
 axis equal
 set(gca,'fontsize',14);
 xlabel('x  [ m ]','fontsize',14);
 ylabel('y  [ m ]','fontsize',14);
 set(gca,'zLim',[-1 1]);
 set(gca,'xLim', [-2,2]);
 set(gca,'yLim', [0.5,2.5]);
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
%%
figure(2)
  pos = [0.5 0.2 0.35 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gca,'fontsize',14);
  
  %subplot(1,2,1)
  xP = x; yP = y; zP = 1-(s.*conj(s)); zP =  (s.*conj(s)).^(-0.1); 
  pcolor(xP,yP,zP);
  colormap autumn
  shading flat; 
  axis equal
  set(gca,'yLim',[0.5,2.5]);
  set(gca,'fontsize',14);
  xlabel('x  [ m ]','fontsize',14);
  ylabel('y  [ m ]','fontsize',14);
  hTitle = title('intensity');
  set(hTitle,'FontWeight','normal');
  
%   subplot(1,2,2)
%   xP = x; yP = y; zP =  (s.*conj(s)).^(-0.1); %-log(s.*conj(s));
%   pcolor(xP,yP,zP);
%   colormap autumn
%   shading flat; 
%   axis equal
%   set(gca,'yLim',[0.5,2.5]);
%   set(gca,'fontsize',14);
%   xlabel('x  [ m ]','fontsize',14);
%   ylabel('y  [ m ]','fontsize',14);
%   hTitle = title('intensity (log scaling)');
%   set(hTitle,'FontWeight','normal');
%  % set(gca,'zLim',[0 2]);
%   
 
 toc