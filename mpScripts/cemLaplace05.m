% cemLaplace05.m
% 8 april 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% Solution of Poisson's equation for different cases
% potential / electric field 
% Default values shown between [ ... ]
% May have to change code for plots with different cases


clear all
close all
clc
tic

% INPUTS  ================================================================

% Number of XY grid points ( ODD integer ) 
     Nx = 101;   % [101]  
     Ny = 101;   % [101]

% Lx X length of region / Ly  Y length of region
     Lx = 10;    % [10] 
     Ly = 10;    % [10]

% tolerance for ending iterations
     tol = 0.001; % [0.01] 

% Goto SETUP section for defining: boundary conditions / initial V / charge density
     flag = 5;
               % flag == 1: uniform charge density distribution
               % flag == 2; rho = eps0 x y
               % flag == 3: constant voltage at origin
               % flag == 4: charge Q at origin
               % flag == 5: two points held at a constant voltage
               % flag == 6: central square region held at a constant
               %            voltage
               % flag == 7: centre square region - insulator
               %            constant charge density 
               % flag == 8; Linear variation of boundary conditions -->
               %            uniform electric field
               % flag == 9; same as 8, except gradient for boundary
               %            conditions   NEUMANN CONDITIONS
               
% SETUP ==================================================================
    
% Nx and Ny must be odd
    if mod(Nx,2) == 0; Nx = Nx + 1; end;
    if mod(Ny,2) == 0; Ny = Ny + 1; end;
      
% XY coordinate system
    minX1  = -Lx/2; maxX1 = Lx/2;
    minY1  = -Ly/2; maxY1 = Ly/2;
    x = linspace(minX1, maxX1,Nx);
    y = linspace(minY1, maxY1,Ny);
    [xx, yy] = meshgrid(x,y);
    hx = x(2) - x(1); hy = y(2) - y(1);
    Kx = hx^2/(2*(hx^2+hy^2)); Ky = hy^2/(2*(hx^2+hy^2));
    
% DEFINE: boundary conditions / initial potential / charge density
     eps0 = 8.854e-12;    % constant   permittivity of free space
     K = hx^2 * hy^2 /(2*(hx^2 + hy^2)*eps0);
     
     V = zeros(Ny,Nx);  
     rho_E = zeros(Ny,Nx);
     
     switch flag
         
         case 1   % uniform charge density distribution
           rho = eps0 .* ones(Ny,Nx);
           rho_E = K .* rho;
          
         case 2   % rho = eps0 x y
           rho = eps0 .* xx .* yy;
           rho_E = K .* rho; 
           
         case 3    % constant voltage at origin
            V1 = 100;
            indx = ceil(Nx/2); indy = ceil(Ny/2);
            V(indy,indx) = V1;
                       
         case 4     % charge Q at origin
            indx = ceil(Nx/2); indy = ceil(Ny/2);
            Q = 1.0e-9;
            rho = Q / (hx*hy);
            rho_E(indy,indx) = K .* rho; 
            
         case 5
             V1 = 10; V2 = -10;
            % indx1 = ceil(Nx/3); indx2 = ceil(2*Nx/3);
             indy1 = ceil(Nx/2); indy2 = ceil(Nx/2);
             indx1 = 41; indx2 = 61;
             V(indy1,indx1) = V1; V(indy2,indx2) = V2;
         
         case 6
             V1 = 0;  V2 = 1.4;
             iS = 41 : 61;
             V(:,1)   = V1;    % outer square
             V(:,end) = V1;
             V(1,:)   = V1;
             V(end,:) = V1;
             V(iS,iS) = V2;    % inner square
            
         case 7
              iS = 41 : 61;
               rho = eps0 .* ones(Ny,Nx);
               rho_E(iS,iS) = K .* rho(iS,iS);
               
         case 8
             V1 = 10; V2 = 5;
             V(:,1) = V1;
             V(:,Nx) = V2;
             m = (V2-V1)/Lx;
             b = V(1,1) - m * x(1);
             V(1,:) = m .* x + b;
             V(Ny,:) = V(1,:);
             
         case 9
             V1 = 10; V2 = 5;
             V(:,1) = V1;
             V(:,Nx) = V2;
             m = (V2-V1)/Lx;
             b = V(1,1) - m * x(1);
             V(1,:) = m .* x + b;
             V(Ny,:) = V(1,:);
     end  
    
% CALCULATIONS ===========================================================  
    
% dSum difference in sum of squares  /  n  number of iterations
    dSum = 1; n = 0;

while  dSum > tol
    sum1 =  sum(sum(V.^2));
        
    for ny = 2: Ny-1
            
        for nx = 2: Nx-1
            if flag == 3; V(indy,indx) = V1; end; 
            if flag == 5; V(indy1,indx1) = V1; V(indy2,indx2) = V2; end; 
            if flag == 6; V(iS,iS) = V2; end;
            if flag == 9; V(ny,Nx) = V(ny,Nx-1) + hx * y(ny)^2; end
            
            V(ny,nx) = Ky * (V(ny,nx+1) + V(ny,nx-1)) + Kx * (V(ny+1,nx) + V(ny-1,nx)) + rho_E(ny,nx);
        
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
 %  Ey1 = Exx(:,1);
 %  Ey2 = Exx(iS,ind1);
 %  Ex21 = Exx(ind1,1:ind1); 
 %  V21 = simpson1d(Ex21,minX1,minX2);

% COMMAND WINDOW DISPLAY =================================================
    % case 3: constant volatage at origin
       if flag == 3
          Q = 4*pi*eps0*hy*V(indy+1,indx)
       end


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
    h =  colorbar;
    h.Label.String = 'V   [ V ]';
    view(55,49);
    %set(gca,'ZTick',[0 5 10]);
    
%%
   figure(2)
       set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32]);
       c = 0; yStep = zeros(1,6);
       
   % for n = 1 : round(Ny/15): 1+ceil(Ny/2)
     for n = [1 11 21 31 41 51]
        xP = xx(n,:); yP = V(n,:);
        plot(xP,yP,'linewidth',2);
        hold on
        c = c+1;
        yStep(c) = yy(n,1);
     end
     
       %tt = num2str(yStep,2);
       %tm1 = 'y  =  ';
       %tm2 = tt;
       %tm = [tm1 tm2];
       tm = 'Potential profiles for different y values';
       xlabel('x  [m]'); ylabel('V  [ V ]');
       set(gca,'fontsize',16)  
       h_title = title(tm);
       set(h_title,'fontsize',12,'FontWeight','normal');
       set(gca,'xTick',-5:5);
       tm1 = num2str(yStep(1),3);
       tm2 = num2str(yStep(2),3); tm3 = num2str(yStep(3),3);
       tm4 = num2str(yStep(4),3);  tm5 = num2str(yStep(5),3);
       tm6 = num2str(yStep(6),3);
       h = legend(tm1,tm2,tm3,tm4,tm5,tm6,'Orientation','horizontal');
       set(h,'Location','north');
       
 %%   
 figure(3)
     set(gcf,'units','normalized','position',[0.65 0.1 0.3 0.32]);
     
     hold on
     for sx = -5 : 1 : 5;
     for sy = -5 : 1 : 5;
        % if sy ~= 0; 
         h = streamline(xx,yy,Exx,Eyy,sx,sy);
         set(h,'linewidth',1,'color',[1 0 1]);
        % end
     end
     end
     

     index1 = 1 : 10 : Nx-1; index2 = 1 : 10 : Ny-1;
    
     %index1 = [21 26 31 36 41 46 61 66 71 76 81 86 91]; index2 = index1;
     p1 = xx(index1, index2); p2 = yy(index1, index2);
     p3 = Exx(index1, index2); p4 = Eyy(index1, index2); 
     h = quiver(p1,p2,p3,p4,'autoscalefactor',1);
     set(h,'color',[0 0 1],'linewidth',2)
     xlabel('x  [m]'); ylabel('y  [m]');
     
      title('electric field','fontweight','normal');

      h = rectangle('Position',[minX1,minY1,2*maxX1,2*maxY1]');
      set(h,'Edgecolor',[1 0.7 0],'linewidth',2);
      
      if flag == 7
         h = rectangle('Position',[-1,-1,2,2]');
         set(h,'Edgecolor',[0 1 0],'FaceColor', [0 1 0],'linewidth',2);
      end
      
     axis equal
     set(gca,'xLim',[-0.5 + minX1, 0.5 + maxX1]);
     set(gca,'yLim',[-0.5 + minY1, 0.5 + maxY1]);
     set(gca,'xTick',minX1:2:maxX1);
     set(gca,'yTick',minY1:2:maxY1);
     set(gca,'fontsize',14)
     box on
     
     if flag == 5
        h = plot(x(41),y(51),'ro');
        set(h,'MarkerFaceColor',[1 0 0],'MarkerSize',8);
        h = plot(x(61),y(51),'ko');
       set(h,'MarkerFaceColor',[0 0 0],'MarkerSize',8);
     end
     
    if flag == 7
     figure(9)
       set(gcf,'units','normalized','position',[0.4 0.4 0.3 0.32]);
       index1 = 41 : 2 : 61; index2 = 41 : 2 : 61;
       p1 = xx(index1, index2); p2 = yy(index1, index2);
       p3 = Exx(index1, index2); p4 = Eyy(index1, index2); 
       h = quiver(p1,p2,p3,p4);
       set(h,'color',[0 0 1],'linewidth',3)
       xlabel('x  [m]'); ylabel('y  [m]');
       h = rectangle('Position',[-1,-1,2,2]');
       set(h,'Edgecolor',[0 1 0],'linewidth',3);
       h = title('electric field inside insulator','fontweight','normal');
       set(h,'fontsize',14);
       axis equal
       set(gca,'xTick',-1:0.5:1);
       set(gca,'yTick',-1:0.5:1);
       set(gca,'fontsize',14);
       set(gca,'xLim',[-1.5, 1.5]);
       set(gca,'yLim',[-1.5, 1.5]);
    end
    
     
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
    h =  colorbar;
    h.Label.String = '| E |   [ V/m ]';
    view(43,54);
    title('electric field','fontweight','normal');
    
%%    
figure(5)
    set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
    contourf(xx,yy,V,16);
    %pcolor(xx,yy,V);
    shading interp
    xlabel('x  [m]'); ylabel('y  [m]');
    title('potential','fontweight','normal');
    set(gca,'fontsize',14)
    box on
    h =  colorbar;
    h.Label.String = 'V   [ V ]';
    %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
    %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
    axis square
    box on
    
%%    
figure(6)
    set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
    contourf(xx,yy,E,20);
    shading interp
    xlabel('x  [m]'); ylabel('y  [m]');
    title('electric field','fontweight','normal');
    set(gca,'fontsize',14)
    box on
    h =  colorbar;
    h.Label.String = '| E |   [ V/m ]';
    axis equal    

%%    
   toc
