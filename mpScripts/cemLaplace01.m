% cemLaplace01.m
% 18 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Solution of Laplace's equation


clear all
close all
clc
tic

% INPUTS  ================================================================
% Number of XY grid points (integer)
     Nx = 7;  Ny = 7;  
% Range for X and Y values  [minX maxX minY maxY minZ maxZ]
     XY = [0 2 0 1];
% tolerance for ending iterations
     tol = 0.9;
% number of iteriations 
     maxN = 50;
% limits for electric field plot
   minXL = 0   ; maxXL = 2;
   minYL = 0; maxYL = 1;

% Boundary values for the potential   
    V0 = 10;
    V = zeros(Ny,Nx);
    V(:,1) = V0;
    V(:,end) = V0;
    V(1,:) = V0;
    V(end,:) = V0;
    
% SETUP ==================================================================
    
    minX  = XY(1); maxX = XY(2);
    minY  = XY(3); maxY = XY(4);
    x = linspace(minX, maxX,Nx);
    y = linspace(minY, maxY,Ny);
    [xx, yy] = meshgrid(x,y);
    hx = x(2) - x(1); hy = y(2) - y(1);
    Kx = hx^2/(2*(hx^2+hy^2)); Ky = hy^2/(2*(hx^2+hy^2));
    
  
    
    % CALCULATIONS ===========================================================  

% dSum difference in sum of squares  /  n  number of iterations
dSum = 1; n = 0;  cc = 0;

while n < maxN
%while  dSum > tol
    sum1 =  sum(sum(V.^2));
        
    for ny = 2: Ny-1
    
        for nx = 2: Nx-1
            V(ny,nx) = Ky * (V(ny,nx+1) + V(ny,nx-1)) + Kx * (V(ny+1,nx) + V(ny-1,nx));
       % end
    
        %    V(ny,nx) = Ky * (V(ny,nx+1) + V(ny,nx-1)) + Kx * (V(ny+1,nx) + V(ny-1,nx));
        end
      %  V(ny,Nx) = Ky * (V(ny,Nx) + V(ny,Nx-2)) + Kx * (V(ny+1,Nx) + V(ny-1,Nx)); 
    end
   
   sum2 =  sum(sum(V.^2));
   dSum = abs(sum2 - sum1);
   n = n + 1; cc = cc+1;
end
  % electric field  
  [Exx, Eyy] = gradient(V,hx,hy);
  Exx = -Exx;  Eyy = -Eyy; 
  E = sqrt(Exx.^2 + Eyy.^2);
  
% GRAPHICS ===============================================================
   
figure(1) % VECTOR FIELD V  -----------------------------------------
    set(gcf,'units','normalized','position',[0.05 0.5 0.4 0.35]);
    surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),V(1:5:end,1:5:end));
    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('V  [ V ]');
    set(gca,'fontsize',14)
    rotate3d 
    box on
    axis tight
    colorbar
    view(43,54);
    
figure(2)
       set(gcf,'units','normalized','position',[0.5 0.5 0.4 0.35]);
       c = 1; yStep = zeros(1,8);
    for n = 1 : round(Ny/15): round(Ny/2)
        xP = xx(n,:); yP = V(n,:);
        plot(xP,yP,'b','linewidth',1.2);
        hold on
        c = c+1;
        yStep(c) = yy(c,1);
        
    end
       tt = num2str(yStep,2);
       text(0.2,95,'y  =');
       text(0.4,95,tt,'fontsize',13)
       xlabel('x  [m]'); ylabel('V  [ V ]');
       set(gca,'fontsize',16)    
    
 figure(3)
     set(gcf,'units','normalized','position',[0.05 0.05 0.4 0.35]);
     index1 = 1 : Nx; index2 = 1 : Ny;
     p1 = xx(index1, index2); p2 = yy(index1, index2);
     p3 = Exx(index1, index2); p4 = Eyy(index1, index2); 
     h = quiver(p1,p2,p3,p4);
     set(h,'color',[0 0 1],'linewidth',2)
     xlabel('x  [m]'); ylabel('y  [m]');
     set(gca,'xLim',[minXL, maxXL]);
     set(gca,'yLim',[minYL, maxYL]);
     
 figure(4)
     set(gcf,'units','normalized','position',[0.5 0.05 0.4 0.35]);
    surf(xx,yy,E);
    shading interp
    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('|E|  [ V/m ]');
    set(gca,'fontsize',14)
    set(gca,'zTick',[0 1000 2000]);
    rotate3d 
    box on
    axis tight
    colorbar
    view(43,54);

    
    cc
    V

toc
