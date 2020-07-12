% cemDiff05.m
% 18 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% CALCULATION OF THE LAPLACIAN
%     [2D] scalar function  flag == 1;
%     [3D] vector function  flag == 2;
   

clear all
close all
clc

tic

% =========================================================================
% INPUTS

% Number of grid points (integer)    N = 101 default value
    N = 201;
% Range for X, Y and Z values  [minX maxX minY maxY minZ maxZ]
    % XY = [-50 50 -50 50 -50 50 ];
     XY = [-1 1 -1 1 -1 1];
 % flag = 1 for scalar field f or flag = 2 for vector field (V(Vx Vy Vz)
    flag = 2;
 
% ========================================================================
if flag == 1;  % [2D] scalar function f
    
% CALCULATIONS  ----------------------------------------------------------   
% X and Y ranges
minX = XY(1); maxX = XY(2); minY  = XY(3); maxY = XY(4);
x = linspace(minX, maxX,N);
y = linspace(minY, maxY,N);
[xx, yy] = meshgrid(x,y);
hx = x(2)-x(1); hy = y(2)-y(1);

% Define scalar field f
       wLx = 50; wLy = 40;
       kx = 2*pi/wLx; ky = 2*pi/wLy;
       f = sin(kx.*xx) .* cos(ky.*yy);
%      f = xx.^2 + yy.^2;

% Laplacian
       del2f = del2(f,hx,hy)*4;
       

% GRAPHICS ---------------------------------------------------------------

figure(1); 
   set(gcf,'units','normalized','position',[0.2 0.2 0.5 0.4]);
   subplot(1,2,1);   % scalar function f
   zz = f;
   surf(xx,yy,zz);
   shading interp;
   colormap(summer)
   set(gca,'fontsize',14)
   axis square
   box on
   view(-19,55);
   title('[2D] scalar function  f(x,y)');
   xlabel('x'); ylabel('y'); zlabel('f'); 
  
   
   subplot(1,2,2)   %  Laplacian del2f  of scalar fucntion f 
   surf(xx,yy,del2f);
   shading interp;
   colormap(autumn)
   set(gca,'fontsize',14)
   axis square
   box on
   view(-19,55);
   title('[2D] Laplacian of f(x,y)');
   xlabel('x'); ylabel('y'); zlabel('del2f');
   colormap(autumn)
   set(gca,'fontsize',14);
   
else
    
% =======================================================================
% Laplacian of the vector field V    flag == 2;
% =======================================================================

% X and Y ranges
    minX  = XY(1); maxX = XY(2);
    minY  = XY(3); maxY = XY(4);
    minZ  = XY(5); maxZ = XY(6);

    x = linspace(minX, maxX,N);
    y = linspace(minY, maxY,N);
    z = linspace(minY, maxY,N);

    [xx, yy, zz] = meshgrid(x,y,z);
    hx = x(2)-x(1); hy = y(2)-y(1); hz = z(2)-z(1);

% Define vector field V
    Vxx = xx.^2;
    Vyy = 3 .* xx .* zz.^2;
    Vzz = -2 .* xx .* zz;
    
% Calculate Laplacian: X Y Z components / magnitude of vector    
    del2Vx = del2(Vxx,hx,hy,hz)*6;
    del2Vy = del2(Vyy,hx,hy,hz)*6;
    del2Vz = del2(Vzz,hx,hy,hz)*6;

    del2Vmag = sqrt(del2Vx.^2 + del2Vy.^2 + del2Vz.^2);
       
% GRAPHICS ---------------------------------------------------------------

 figure(3) % curl --------------------------------------------------------
    set(gcf,'units','normalized','position',[0.7 0.2 0.30 0.4]);
    % indices dx dy dz for slice plots in figures 1 and 3 
    dx = 1:10:N; dy  = dx;   dz = 76;
   
    % plot variable: p1 p2 p3 p4 p5 p6
    p1 = xx(dx,dy,dz); p2 = yy(dx,dy,dz); p3 = zz(dx,dy,dz); 
    p4 = del2Vx(dx,dy,dz); p5 = del2Vy(dx,dy,dz); p6 = del2Vz(dx,dy,dz); 
    
    h = quiver3(p1, p2, p3, p4, p5, p6);
    set(h,'color',[0 0 1],'linewidth', 2, 'MaxHeadSize',1);
    axis tight
    set(gca,'xLim',[minX, maxX]);
    set(gca,'yLim',[minY, maxY]);
    % set(gca,'zLim',[minZ, maxZ]);
    title('Laplacian of a vector field');
    xlabel('x'); ylabel('y'); zlabel('del2');
    set(gca,'fontsize',14)
    rotate3d 
    box on 
    view(0,90);
    
figure(4);  %  ---------------------------------------------------------
   set(gcf,'units','normalized','position',[0.38 0.2 0.3 0.4]);
   % INPUTS:  Hx Hy Hz slice places
   Hx = 0; Hy = 0; Hz = 0.5;
   h = slice(xx,yy,zz,del2Vmag,Hx,Hy,Hz);
   colormap jet
   shading interp
   daspect([1 1 1])
   axis tight
   camlight
   set([h(1),h(2)],'ambientstrength',0.6)
   colorbar;
   rotate3d;
   xlabel('x'); ylabel('y'); zlabel('z');
   title('magnitude of Laplacian of vector field');
   set(gca,'fontsize',14)  
   axis tight    

end

toc
