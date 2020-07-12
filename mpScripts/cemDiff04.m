% cemDiff03.m
% 18 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% DIVERGENCE and CURL OF A VECTOR FIELD
   

clear all
close all
clc

% INPUTS  ================================================================
% Number of grid points (integer)
   N = 101;
% Range for X and Y values  [minX maxX minY maxY minZ maxZ]
  XY = [-1 1 -1 1 -1 1];
 
  
% indices for slice plots in figures 1 and 3 
  dx = 1:20:N; dy  = dx;   dz = 76; 
   
% CALCULATIONS  ==========================================================   
% X and Y ranges
   minX = XY(1);  maxX = XY(2);
   minY  = XY(3); maxY = XY(4);
   minZ  = XY(5); maxZ = XY(6);
   x = linspace(minX, maxX,N);
   y = linspace(minY, maxY,N);
   z = linspace(minZ, maxZ,N);
   [xx, yy, zz] = meshgrid(x,y,z);

% Define vector field V
%    rr = sqrt(xx.^2 + yy.^2);
%    tt = atan2(yy,xx);
%    Vxx = -rr.*sin(tt);
%    Vyy = rr.*cos(tt);
%    Vzz = 0.*zz.^2;
%    Vxx = xx .* yy;  Vyy = yy.^2;  Vzz = 0 .* Vxx; 
     Vxx = yy.^2;   Vyy = 2.*xx.*yy;  Vzz = 2.*yy.*zz;
     
% DIVERGENCE
divV = divergence(xx, yy, zz, Vxx, Vyy, Vzz);

[curlVxx, curlVyy, curlVzz] = curl(xx, yy, zz, Vxx, Vyy, Vzz);


% GRAPHICS ===============================================================

figure(1) % VECTOR FIELD V------------------------------------------------
    set(gcf,'units','normalized','position',[0.05 0.2 0.3 0.4]);
    p1 = xx(dx,dy,dz); p2 = yy(dx,dy,dz); p3 = zz(dx,dy,dz); 
    p4 = Vxx(dx,dy,dz); p5 = Vyy(dx,dy,dz); p6 = Vzz(dx,dy,dz); 
    h = quiver3(p1, p2, p3, p4, p5, p6);
    set(h,'color',[0 0 1],'linewidth', 2);
    axis tight
    set(gca,'xLim',[minX, maxX]);
    set(gca,'yLim',[minY, maxY]);
    % set(gca,'zLim',[minZ, maxZ]);
    title('Vector Field V');
    xlabel('x'); ylabel('y'); zlabel('z');
    set(gca,'fontsize',14)
    rotate3d 
    box on
    axis tight
    
figure(2);  %  divergence ------------------------------------------------
   set(gcf,'units','normalized','position',[0.38 0.2 0.3 0.4]);
   %h = slice(xx,yy,zz,divV,[minX maxX],maxY,minZ);
   h = slice(xx,yy,zz,divV,0,0,0.5);
   colormap jet
   shading interp
   daspect([1 1 1])
   axis tight
   camlight
   set([h(1),h(2)],'ambientstrength',0.6)
   colorbar;
   rotate3d;
   xlabel('x'); ylabel('y'); zlabel('z');
   title('divergence');
   set(gca,'fontsize',14)  
   axis tight
    
 figure(3) % curl --------------------------------------------------------
    set(gcf,'units','normalized','position',[0.7 0.2 0.30 0.4]);
    p1 = xx(dx,dy,dz); p2 = yy(dx,dy,dz); p3 = zz(dx,dy,dz); 
    p4 = curlVxx(dx,dy,dz); p5 = curlVyy(dx,dy,dz); p6 = curlVzz(dx,dy,dz); 
    
    h = quiver3(p1, p2, p3, p4, p5, p6);
    set(h,'color',[0 0 1],'linewidth', 2);
    axis tight
     set(gca,'xLim',[minX, maxX]);
     set(gca,'yLim',[minY, maxY]);
    % set(gca,'zLim',[minZ, maxZ]);
    title('curl');
    xlabel('x'); ylabel('y'); zlabel('z');
    set(gca,'fontsize',14)
    rotate3d 
    box on   