% cemLaplace04.m
% 08 april 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Solution of Laplace's equation
% Two concentric metal squares
% potential / electric field / charge / capactiance calculations
% default values shown between [ ... ]
% Nx, Ny must be odd number for using Simpon's rule
% array iS must have an odd number of elements for using Simpson's rule


clear all
close all
clc
tic

% INPUTS  ================================================================

% Number of XY grid points ( ODD integer ) 
     Nx = 101;   % [101]  
     Ny = 101;   % [101]
% L1 side length of outer square and L2 side length of inner square ]
     L1 = 10;    % [10] 
     L2 =  2;    % [2]
% boundary values for the potential outer square V1 and inner square V2
     V1 = 10;     % [10]
     V2 = 5;      % [ 5]
% tolerance for ending iterations
     tol = 0.01; % [0.01] 
    
        
% SETUP ==================================================================
    
    % Nx and Ny must be odd
      if mod(Nx,2) == 0; Nx = Nx + 1; end;
      if mod(Ny,2) == 0; Ny = Ny + 1; end;
      
    % XY coordinate system
       minX1  = -L1/2; maxX1 = L1/2;
       minY1  = -L1/2; maxY1 = L1/2;
       minX2  = -L2/2; maxX2 = L2/2;
       minY2  = -L2/2; maxY2 = L2/2;
       x = linspace(minX1, maxX1,Nx);
       y = linspace(minY1, maxY1,Ny);
       [xx, yy] = meshgrid(x,y);
       hx = x(2) - x(1); hy = y(2) - y(1);
       Kx = hx^2/(2*(hx^2+hy^2)); Ky = hy^2/(2*(hx^2+hy^2));
    
    % indices for inner square
    % check that iS has an odd number of elements
       ind1 = find(x >= minX2,1);
       ind2 = find(x >= maxX2,1);
       iS = ind1:ind2;
       if mod(length(iS),2) == 0; iS = ind1:ind2+1; end;
       
    % boundary values
       V = zeros(Ny,Nx);
       V(:,1)   = V1;    % outer square
       V(:,end) = V1;
       V(1,:)   = V1;
       V(end,:) = V1;
       V(iS,iS) = V2;    % inner square

% CALCULATIONS ===========================================================  
    eps0 = 8.854e-12;    % constant   permittivity of free space

% dSum difference in sum of squares  /  n  number of iterations
    dSum = 1; n = 0;

while  dSum > tol
    sum1 =  sum(sum(V.^2));
        
    for ny = 2: Ny-1
            V(iS,iS) = V2;
        for nx = 2: Nx-1
            V(ny,nx) = Ky * (V(ny,nx+1) + V(ny,nx-1)) + Kx * (V(ny+1,nx) + V(ny-1,nx));
            V(iS,iS) = V2;
        end
    end
    
   sum2 =  sum(sum(V.^2));
   dSum = abs(sum2 - sum1);
   n = n+1;
end

% electric field / potential / line integral  
   [Exx, Eyy] = gradient(V,hx,hy);
   Exx = -Exx;  Eyy = -Eyy; 
   E = sqrt(Exx.^2 + Eyy.^2);
   Ey1 = Exx(:,1);
   Ey2 = Exx(iS,ind1);
   Ex21 = Exx(ind1,1:ind1); 
   V21 = simpson1d(Ex21,minX1,minX2);
 
% charge Q / capacitance per unit length / 
   Q2 = -4*eps0*simpson1d(Ey2',minY2,maxY2);
   Q1 =  4*eps0*simpson1d(Ey1',minY1,maxY1);
   Q = abs(abs(Q2)-abs(Q1))/2;
   Cap = Q/(V1-V2);

%  charge density one side of inner square
  sigma = Ey2; 
  
  % theoretical capacitance of two concentric cylinders
  b = 8*maxX1/(2*pi); a = 8*maxX2/(2*pi);
  CapT = 2*pi*eps0/log(b/a);
 
% GRAPHICS ===============================================================
   
figure(1) % VECTOR FIELD V  -----------------------------------------
    set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
    surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),V(1:5:end,1:5:end));
    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('V  [ V ]');
    title('potential','fontweight','normal');
    set(gca,'fontsize',14);
    rotate3d 
    box on
    axis tight
    colorbar
    view(55,49);
    set(gca,'ZTick',[0 5 10]);
    

   figure(2)
       set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32]);
       c = 0; yStep = zeros(1,9);
       
    for n = 1 : round(Ny/15): 1+ceil(Ny/2)
        xP = xx(n,:); yP = V(n,:);
        plot(xP,yP,'b','linewidth',1.2);
        hold on
        c = c+1;
        yStep(c) = yy(n,1);
        
    end
       tt = num2str(yStep,2);
       tm1 = 'y  =  ';
       tm2 = tt;
       tm = [tm1 tm2];
       xlabel('x  [m]'); ylabel('V  [ V ]');
       set(gca,'fontsize',16)  
       h_title = title(tm);
       set(h_title,'fontsize',12,'FontWeight','normal');
      
 %%   
 figure(3)
     set(gcf,'units','normalized','position',[0.65 0.1 0.3 0.32]);
     
     hold on
     sx = -5;
     for sy = -5:5
         h = streamline(xx,yy,Exx,Eyy,sx,sy);
         set(h,'linewidth',1,'color',[1 0 1]);
     end
     
     sx = 5;
     for sy = -5:5
         h = streamline(xx,yy,Exx,Eyy,sx,sy);
         set(h,'linewidth',1,'color',[1 0 1]);
     end
     
     sy = -5;
     for sx = -5:5
         h = streamline(xx,yy,Exx,Eyy,sx,sy);
         set(h,'linewidth',1,'color',[1 0 1]);
     end
     
     sy = 5;
     for sx = -5:5
         h = streamline(xx,yy,Exx,Eyy,sx,sy);
         set(h,'linewidth',1,'color',[1 0 1]);
     end
     
     index1 = 1 : 10: Nx; index2 = 1 : 10 : Ny;
     p1 = xx(index1, index2); p2 = yy(index1, index2);
     p3 = Exx(index1, index2); p4 = Eyy(index1, index2); 
     h = quiver(p1,p2,p3,p4);
     set(h,'color',[0 0 1],'linewidth',2)
     xlabel('x  [m]'); ylabel('y  [m]');
     
     title('electric field','fontweight','normal');
     hold on
     h = rectangle('Position',[minX2,minY2,2*maxX2,2*maxY2]');
     set(h,'Edgecolor',[1 0 0],'lineWidth',2);
     h = rectangle('Position',[minX1,minY1,2*maxX1,2*maxY1]');
     set(h,'Edgecolor',[1 0 0],'linewidth',1);
     
     axis equal
     set(gca,'xLim',[-0.5 + minX1, 0.5 + maxX1]);
     set(gca,'yLim',[-0.5 + minY1, 0.5 + maxY1]);
     set(gca,'xTick',minX1:maxX1);
     set(gca,'yTick',minY1:maxY1);
     set(gca,'fontsize',14)
     box on
     
     
%%
 figure(4)
    set(gcf,'units','normalized','position',[0.02 0.1 0.3 0.32]);
    surf(xx,yy,E);
    shading interp
    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('|E|  [ V/m ]');
    set(gca,'fontsize',14)
    %set(gca,'zTick',[0 1000 2000]);
    rotate3d 
    box on
    axis tight
    colorbar
    view(43,54);
    title('electric field','fontweight','normal');
    
figure(5)
    set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
    contourf(xx,yy,V,16);
    shading interp
    xlabel('x  [m]'); ylabel('y  [m]');
    title('potential','fontweight','normal');
    set(gca,'fontsize',14)
    box on
    colorbar
    axis equal
    
figure(6)
    set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
    contourf(xx,yy,E,32);
    shading interp
    xlabel('x  [m]'); ylabel('y  [m]');
    title('electric field','fontweight','normal');
    set(gca,'fontsize',14)
    box on
    colorbar
    axis equal    

figure(7)
    set(gcf,'units','normalized','position',[0.5 0.5 0.3 0.32]);
    plot(x(iS),sigma,'b','linewidth',2);
    xlabel('y  [m]'); ylabel('charge density  [a.u.]');
    title('Inner Square x = -1: Charge density','fontweight','normal');
    set(gca,'fontsize',14)
    box on
    grid on
    
   Q1
   Q2
   Q
   Cap
   CapT
   V21
   tol
    
   toc
