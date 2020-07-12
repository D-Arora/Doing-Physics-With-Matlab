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
     Nx = 101;  Ny = 101;  
% Range for X and Y values  [minX maxX minY maxY minZ maxZ]
     XY = [0 4 -1 1];
% Boundary potential
    
% tolerance for ending iterations
tol = 0.001;
% limits for electric field plot
   minXL = 0; maxXL = XY(2);
   minYL = 0; maxYL = XY(4);
   
% boundary values for the potential
    V0 = 100; V1 = -5; V2 = -10;
    V = zeros(Ny,Nx);
    V(:,1) = V0;
   % V(:,end) = 0;
   % V(1,:) = 0;
   % V(end,:) = 0;

    % linearly increasing potentials on boundaries
       %enter INPUT details in setup
     
       
% SETUP ==================================================================
   
    minX  = XY(1); maxX = XY(2);
    minY  = XY(3); maxY = XY(4);
    x = linspace(minX, maxX,Nx);
    y = linspace(minY, maxY,Ny);
    [xx, yy] = meshgrid(x,y);
    hx = x(2) - x(1); hy = y(2) - y(1);
    Kx = hx^2/(2*(hx^2+hy^2)); Ky = hy^2/(2*(hx^2+hy^2));
    
   % V(:,1) = V0 .* cos(2*pi*y/(4*maxY));
   
% linearly increasing potentials on boundaries
%        m1 = 20/maxY; m2 = 20/maxX; m3 = 30/maxY; m4 = 10/maxX;
%        b1 = 0;       b2 = 20;      b3 = 10;   b4 = 0;
%        
%        V(:,1)   = m1 .* y + b1;
%        V(:,end) = m3 .* y + b3;
%        V(end,:)   = m2 .* x + b2;
%        V(1,:) = m4 .* x + b4;
       
% CALCULATIONS ===========================================================  

% dSum difference in sum of squares  /  n  number of iterations
dSum = 1; n = 0;

while  dSum > tol
    sum1 =  sum(sum(V.^2));
        
    for ny = 2: Ny-1
    
        for nx = 2: Nx-1
            V(ny,nx) = Ky * (V(ny,nx+1) + V(ny,nx-1)) + Kx * (V(ny+1,nx) + V(ny-1,nx));
        end
       % uncomment the next line of code if you need to calculate this column of values 
      % V(ny,Nx) = Ky * (V(ny,Nx) + V(ny,Nx-2)) + Kx * (V(ny+1,Nx) + V(ny-1,Nx)); 
    end
   
   sum2 =  sum(sum(V.^2));
   dSum = abs(sum2 - sum1);
   n = n+1;
end
  % electric field  
  [Exx, Eyy] = gradient(V,hx,hy);
  Exx = -Exx;  Eyy = -Eyy; 
  E = sqrt(Exx.^2 + Eyy.^2);
  
  
 % exact solution for boundary conditions
     % left boundary V0, other boundaries V = 0
    % Vexact = (2*V0/pi) .* atan(sin(pi*y(51)/maxX)./sinh(pi*x/maxX));
   %  Vexact = (2*V0/pi) .* atan(1./sinh(pi*x/maxX));
     
     
 % linearly increasing potentials on boundaries
%        m1 = 20/maxY; m2 = 20/maxX; m3 = 30/maxY; m4 = 10/maxX;
%        b1 = 0;       b2 = 20;      b3 = 10;   b4 = 0;
%        
%        V(:,1)   = m1 .* y + b1;
%        V(:,end) = m3 .* y + b3;
%        V(end,:)   = m2 .* x + b2;
%        V(1,:) = m4 .* x + b4;
 
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
   % set(gca,'ZTick',[0 5 10]);
    
%%
   figure(2)
       set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32]);
       c = 0; yStep = zeros(1,6);
       
    for n = 1:10:51
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
      
       
      % hold on
      % plot(x,Vexact,'r'); 
      
 %%   
 figure(3)
     set(gcf,'units','normalized','position',[0.65 0.1 0.3 0.32]);
     
     hold on
     sx = x(5);
         for sy = -1 : 0.1 : 1
           h = streamline(xx,yy,Exx,Eyy,sx,sy);
           set(h,'linewidth',2,'color',[1 0 1]);
         end
     
     
     index1 = 10 : 1: Nx; index2 = 10 : 1 : Ny;
     p1 = xx(index1, index2); p2 = yy(index1, index2);
     p3 = Exx(index1, index2); p4 = Eyy(index1, index2); 
     h = quiver(p1,p2,p3,p4);
     set(h,'color',[0 0 1],'linewidth',2)
     xlabel('x  [m]'); ylabel('y  [m]');
     
     title('electric field','fontweight','normal');
     hold on
     axis equal
     set(gca,'xLim',[minX,maxX]);
     set(gca,'yLim',[minY,maxY]);
     set(gca,'xTick',minX:1:maxX);
     set(gca,'yTick',minY:1:maxY);
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


toc
