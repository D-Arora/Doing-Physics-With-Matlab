% cemLaplace01.m
% 18 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% Solution of Laplace's equation



clear all
close all
clc
tic

% INPUTS  ================================================================

% Number of XY grid points (integer)
     Nx = 101;  
% Range for X and Y values  [minX maxX minY maxY minZ maxZ]
     XY = [0 5];
% Boundary potential
     V0 = 100;
% tolerance for ending iterations
    tol = 0.1;

% SETUP ==================================================================

    V = zeros(Nx,1);
    minX = XY(1);  maxX = XY(2);
    x = linspace(minX, maxX,Nx);
    hx = x(2) - x(1); 
    V(1) = V0;  V(Nx) = 0;

% CALCULATIONS ===========================================================

% dSum difference in sum of squares  /  n  number of iterations
dSum = 100; n = 0;
while dSum > tol
    sum1 =  sum(sum(V.^2));
    
    for nx = 2: Nx-1
           V(nx) = 0.5*(V(nx+1) + V(nx-1));
    end
   
    sum2 =  sum(sum(V.^2));
    dSum = abs(sum2 - sum1);
    n = n+1;
end

% exact solution
VT = -20.*x + 100; 

% electric field
E = -gradient(V,hx);
    
% GRAPHICS ===============================================================
   
figure(1) 
    set(gcf,'units','normalized','position',[0.05 0.2 0.3 0.4]);
    subplot(2,1,1)
    plot(x,V,'b','lineWidth',2);
    hold on
    plot(x,VT,'r','lineWidth',2);
    %set(gca,'xLim',[minX, maxX]);
    %set(gca,'yLim',[minY, maxY]);
    xlabel('x  [ m ]'); ylabel('V  [ V ]'); 
    set(gca,'fontsize',14)
    
    subplot(2,1,2)
    plot(x,E,'b','lineWidth',2)
    set(gca,'yLim',[0, 25]);
    xlabel('x  [ m ]'); ylabel('E  [ V / m ]'); 
    set(gca,'fontsize',14)
    
 % =======================================================================  
    toc
    