% cemVE01.m
% 8 april 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Potential and Electric Field in [2D] region
% max 5 charges


clear all
close all
clc
tic

% INPUTS  ================================================================

% Number of grid point    [N = 10001]
   N = 1001;
   
% Charge  Q = [10, 0, 0, 0, 0] 
   Q = [-20, 0, 0, 0, 0] .* 1e-6;
   
% X & Y components of position of charges  [0, 0, 0, 0, 0]
   xC = [0,  -0.5,  0.5, -0.5, 0];
   yC = [0,  -0.5, -0.5, 0.5, 0];

% 5 random charges   uncomment to run the program for 5 random charges
%   Q = (1 + 9 .* rand(5,1)) .* 1e-6;
%   xC = -2 + 4 .* rand(5,1); 
%   yC = -2 + 4 .* rand(5,1); 
   
% Dimensions of region   [dimensions of region -2 to 2 / minR = 1e-6 / Esat = 1e6 / Vsat = 1e6]
   minX = -2;  
   maxX =  2;
   minY = -2;
   maxY =  2;
   minR = 1e-6;
   minRx = 1e-6;
   minRy = 1e-6;
   Vsat = 1e6;
   Esat = 1e6;
   
% SETUP  =================================================================

eps0 = 8.854e-12;
kC = 1/(4*pi*eps0);
V = zeros(N,N);
Ex = zeros(N,N); Ey = zeros(N,N);

x  = linspace(minX,maxX,N);
y = linspace(minY, maxY,N);

[xG, yG] = meshgrid(x,y);


% CALCULATION: POTENTIAL & ELECTRIC FIELD ================================

for n = 1 : 5
   Rx = xG - xC(n);
   Ry = yG - yC(n);
   
   index = find(abs(Rx)+ abs(Ry) == 0); 
   Rx(index) = minRx;  Ry(index) = minRy;
   
   R = sqrt(Rx.^2 + Ry.^2);
   R(R==0) = minR;
   V = V + kC .* Q(n) ./ (R);
   
   R3 = R.^3;
   Ex = Ex + kC .* Q(n) .* Rx ./ R3;
   Ey = Ey + kC .* Q(n) .* Ry ./ R3;
end

   if max(max(V)) >=  Vsat; V(V > Vsat)  = Vsat; end;
   if min(min(V)) <= -Vsat; V(V < -Vsat) = -Vsat; end;

   E = sqrt(Ex.^2 + Ey.^2);
   if max(max(E)) >=  Esat; E(E >  Esat)  =  Esat; end;
   if min(min(E)) <= -Esat; E(E < -Esat)  = -Esat; end;
   

% GRAPHICS ===============================================================
figure(1)
   set(gcf,'units','normalized','position',[0.02 0.52 0.28 0.32]);
   surf(xG,yG,V);
   colorbar
   shading interp
   axis square
   set(gca,'xLim',[-2, 2]); set(gca,'yLim',[-2, 2]);
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('V  [ V ]');
   title('potential','fontweight','normal');
   box on
   h =  colorbar;
   h.Label.String = 'V   [ V ]';
   rotate3d
   view(30,50);
   set(gca,'fontsize',12)
   
figure(2)
   set(gcf,'units','normalized','position',[0.31 0.52 0.28 0.32]);
   contourf(xG,yG,V,16)
   shading interp
   colorbar
   xlabel('x  [m]'); ylabel('y  [m]');
   title('potential','fontweight','normal');
   h =  colorbar;
   h.Label.String = 'V   [ V ]';
   %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
   %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
   axis square
   set(gca,'fontsize',12)
   box on
   
figure(3)
   set(gcf,'units','normalized','position',[0.60 0.52 0.28 0.32]);
       c = 0; yStep = zeros(1,5);
       
      for n = ceil(N/2) : ceil(N/10): N
   %  for n = [1 11 21 31 41 51]
        xP = xG(n,:); yP = V(n,:);
        plot(xP,yP,'linewidth',2);
        hold on
        c = c+1;
        yStep(c) = yG(n,1);
     end
     
       %tt = num2str(yStep,2);
       %tm1 = 'y  =  ';
       %tm2 = tt;
       %tm = [tm1 tm2];
       tm = 'Potential profiles for different y values';
       xlabel('x  [m]'); ylabel('V  [ V ]');
       set(gca,'fontsize',12)  
       h_title = title(tm);
       set(h_title,'fontsize',12,'FontWeight','normal');
       %set(gca,'xTick',-5:5);
       tm1 = num2str(yStep(1),2);
       tm2 = num2str(yStep(2),2); tm3 = num2str(yStep(3),2);
       tm4 = num2str(yStep(4),2);  tm5 = num2str(yStep(5),2);
       %tm6 = num2str(yStep(6),3);
       h = legend(tm1,tm2,tm3,tm4,tm5,'Orientation','horizontal');
       set(h,'Location','north');
      % 

figure(4)
   set(gcf,'units','normalized','position',[0.02 0.1 0.28 0.32]);
   surf(xG,yG,E);
   colorbar
   shading interp
   axis square
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('E  [ V / m ]');
   title('electric field | E |','fontweight','normal');
   box on
   h =  colorbar;
   h.Label.String = '| E |   [ V/m ]';
   rotate3d
   view(30,50);
   set(gca,'fontsize',12) 
   
figure(5)
   set(gcf,'units','normalized','position',[0.31 0.10 0.28 0.32]);
   contourf(xG,yG,E,16)
   shading interp
   colorbar
   xlabel('x  [m]'); ylabel('y  [m]');
   title('electric field | E | ','fontweight','normal');
   h =  colorbar;
   h.Label.String = '|E|   [ V/m ]';
   %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
   %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
   axis square
   set(gca,'fontsize',12)
   box on   
   
 figure(6)
   set(gcf,'units','normalized','position',[0.60 0.10 0.28 0.32]);
       c = 0; yStep = zeros(1,5);
       
      for n = ceil(N/2) : ceil(N/10): N
   %  for n = [1 11 21 31 41 51]
        xP = xG(n,:); yP = E(n,:);
        plot(xP,yP,'linewidth',2);
        hold on
        c = c+1;
        yStep(c) = yG(n,1);
     end
     
       %tt = num2str(yStep,2);
       %tm1 = 'y  =  ';
       %tm2 = tt;
       %tm = [tm1 tm2];
       tm = '| E | profiles for different y values';
       xlabel('x  [m]'); ylabel('| E |  [ V/m ]');
       set(gca,'fontsize',12)  
       h_title = title(tm);
       set(h_title,'fontsize',12,'FontWeight','normal');
       %set(gca,'xTick',-5:5);
       tm1 = num2str(yStep(1),2);
       tm2 = num2str(yStep(2),2); tm3 = num2str(yStep(3),2);
       tm4 = num2str(yStep(4),2);  tm5 = num2str(yStep(5),2);
       %tm6 = num2str(yStep(6),3);
       h = legend(tm1,tm2,tm3,tm4,tm5,'Orientation','horizontal');
       set(h,'Location','north');  
%%       
figure(7)

     hold on
     index1 = 1 : 5 : N; index2 = 1 : 5 : N;
    
     %index1 = [21 26 31 36 41 46 61 66 71 76 81 86 91]; index2 = index1;
     p1 = xG(index1, index2); p2 = yG(index1, index2);
     p3 = Ex(index1, index2)./(E(index1,index2)); p4 = Ey(index1, index2)./(E(index1,index2));
    % p3 = Ex(index1, index2); p4 = Ey(index1, index2);  
     h = quiver(p1,p2,p3,p4,'autoscalefactor',0.8);
     set(h,'color',[0 0 1],'linewidth',2)
     xlabel('x  [m]'); ylabel('y  [m]');
     
     hold on
     for sx = -2 : 0.25 : 2;
     for sy = -2 : 0.25 : 2;
        % if sy ~= 0; 
         h = streamline(xG,yG,Ex,Ey,sx,sy);
         set(h,'linewidth',1,'color',[1 0 1]);
        % end
     end
     
     end
     set(gca,'xLim',[-2,2]); set(gca,'yLim', [-2, 2]);
     axis([-2 2 -2 2]);
     axis equal
     box on
 
     %%
     figure(99)
     sx = [0 -1 1 0 1]; sy = [-1 -1 0 1 1];
     h = streamline(xG,yG,Ex,Ey,sx,sy);
     
   