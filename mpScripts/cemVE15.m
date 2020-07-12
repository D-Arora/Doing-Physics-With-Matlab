% cemVE14.m
% 28 April 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% ELECTRIC FIELD
%    GAUSS'S LAW: charges enclosed in a cube
%    calculation of total flux through cube

% SI units used unless stated otherwise

clear all
close all
clc

tic

% ========================================================================
% INPUTS  
% ========================================================================

% Number of grid point for each face of cube  [N = 101]
%  must be an ODD number
   N = 101;  
% Charges enclosed inside sphere:  5 charges
  NQ = 5;
  Q  = [3, -1, 1, -5, 2] .* 1e-6;
  % xQ = [1, 2, 3, 5, 8];
  % yQ = [5, 5, 5, 5, 5];
  % zQ = [5, 5, 5, 5, 5];
  a = 0.1; b = 9.9;
  xQ = a + (b - a) .* rand(1,5);
  yQ = a + (b - a) .* rand(1,5);
  zQ = a + (b - a) .* rand(1,5);

% Dimensions of cube
   L = 10;
   
   
% =======================================================================
% SETUP 
% =======================================================================

% constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);

% Grid spacing
   G = linspace(0,L,N);
   [G1, G2] = meshgrid(G,G); 
   
   flux = zeros(1,6);
   
% =======================================================================
% CALCULATION: ELECTRIC FIELD 
% =======================================================================

% Side 1: nx0
    p1 = 0; p2 = G1; p3 = G2;
    E = zeros(N,N);
    for n = 1 : NQ
       Rx = p1 - xQ(n);
       Ry = p2 - yQ(n);
       Rz = p3 - zQ(n);
       R = sqrt(Rx.^2 + Ry.^2+ Rz.^2);
       R3 = R.^3;
       E = E + (kC * Q(n) * Rx) ./R3;
    end
      flux(1) = -simpson2d(E,0,L,0,L);
            
% Side 2: nx1
    p1 = L; p2 = G1; p3 = G2;
    E = zeros(N,N);
    for n = 1 : NQ
       Rx = p1 - xQ(n);
       Ry = p2 - yQ(n);
       Rz = p3 - zQ(n);
       R = sqrt(Rx.^2 + Ry.^2+ Rz.^2);
       R3 = R.^3;
       E = E + (kC * Q(n) * Rx) ./R3;
    end
    flux(2) = simpson2d(E,0,L,0,L);
        
 % Side 3: ny0
    p1 = G1; p2 = 0; p3 = G2;
    E = zeros(N,N);
    for n = 1 : NQ
       Rx = p1 - xQ(n);
       Ry = p2 - yQ(n);
       Rz = p3 - zQ(n);
       R = sqrt(Rx.^2 + Ry.^2+ Rz.^2);
       R3 = R.^3;
       E = E + (kC * Q(n) * Ry) ./R3;
    end   
      flux(3) = -simpson2d(E,0,L,0,L);
      
% Side 4: ny1
    p1 = G1; p2 = L; p3 = G2;
    E = zeros(N,N);
    for n = 1 : NQ
       Rx = p1 - xQ(n);
       Ry = p2 - yQ(n);
       Rz = p3 - zQ(n);
       R = sqrt(Rx.^2 + Ry.^2+ Rz.^2);
       R3 = R.^3;
       E = E + (kC * Q(n) * Ry) ./R3;
    end   
      flux(4) = simpson2d(E,0,L,0,L);
          
 % Side 5: nz0
    p1 = G1; p2 = G2; p3 = 0;
    E = zeros(N,N);
    for n = 1 : NQ
       Rx = p1 - xQ(n);
       Ry = p2 - yQ(n);
       Rz = p3 - zQ(n);
       R = sqrt(Rx.^2 + Ry.^2+ Rz.^2);
       R3 = R.^3;
       E = E + (kC * Q(n) * Rz) ./R3;
    end   
      flux(5) = -simpson2d(E,0,L,0,L);
      
% Side 6: nz1
    p1 = G1; p2 = G2; p3 = L;
    E = zeros(N,N);
    for n = 1 : NQ
       Rx = p1 - xQ(n);
       Ry = p2 - yQ(n);
       Rz = p3 - zQ(n);
       R = sqrt(Rx.^2 + Ry.^2+ Rz.^2);
       R3 = R.^3;
       E = E + (kC * Q(n) * Rz) ./R3;
    end   
      flux(6) = simpson2d(E,0,L,0,L);
          
   eps0Flux = eps0 * sum(flux);

   Qenclosed = eps0Flux
   Qtotal = sum(Q)
   flux
   
% % ======================================================================= 
% % GRAPHICS 
% % ======================================================================= 
   fs = 12;

figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
   plot3([0 0],[0 0],[0 L],'k');
   hold on
   plot3([0 0],[0 L],[L L],'k');
   plot3([0 L],[0 0],[L L],'k');
   
   for n = 1 : NQ
      if Q(n) > 0; col = [1 0 0]; else col = [0 0 0]; end
       h = plot3(xQ(n),yQ(n),zQ(n),'o');
       set(h,'markerFaceColor',col,'markerEdgeColor',col);
       set(h,'markersize',8);
   end
   
   xlabel('x','fontsize',12); ylabel('y','fontsize',12); zlabel('z','fontsize',12);
   
   grid on
   box on
   
   set(gca,'fontsize',12)
   
   
%    zP = V./1e6;
%    contourf(xG,yG,zP,16);
%    
% % titles
%    xlabel('x  [m]'); ylabel('y  [m]'); 
%    title('potential','fontweight','normal');
%    
% % charged conductors 
%    col = [0 0 0];
%    pos = [-RS, -RS, 2*RS, 2*RS];
%    h = rectangle('Position',pos,'Curvature',[1,1]);
%    set(h,'FaceColor',col,'EdgeColor',col);
%       
%    col = [1 0 0]; b = 0.2;
%    pos = [-b+xQ, -b+yQ, 2*b, 2*b];
%    h = rectangle('Position',pos,'Curvature',[1,1]);
%    set(h,'FaceColor',col,'EdgeColor',col);
%       
% % graphics parameters
%    shading interp
%    h = colorbar;
%    h.Label.String = 'V   [ MV ]';
%    h = colormap('parula');
%    beta = 0.5;
%    brighten(h,beta);
%    set(gca,'fontsize',fs);
%    axis equal
%       
% %%
% figure(2)   % 22222222222222222222222222222222222222222222222222222222222
%    set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
%    
%    contourf(xG,yG,E./1e6,16);
%    
%    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('E  [ V / m ]');
%    title('electric field | E |','fontweight','normal');
%    
%    colorbar
%    shading interp
%    h =  colorbar;
%    h.Label.String = '| E |   [ MV/m ]';
%    
% % charged conductors 
%    col = [0 0 0];
%    pos = [-RS, -RS, 2*RS, 2*RS];
%    h = rectangle('Position',pos,'Curvature',[1,1]);
%    set(h,'FaceColor',col,'EdgeColor',col);
%       
%    col = [1 0 0]; b = 0.2;
%    pos = [-b+xQ, -b+yQ, 2*b, 2*b];
%    h = rectangle('Position',pos,'Curvature',[1,1]);
%    set(h,'FaceColor',col,'EdgeColor',col);
%    
%    axis square
%    set(gca,'fontsize',fs) 
%    box on
% 
%       
% %% 
% figure(3)   % 33333333333333333333333333333333333333333333333333333333333
%    set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
%   
%    xP = y; yP = V(:,701)./1e6;
%    h = plot(xP,yP);
%    set(h,'lineWidth',2);
%    hold on
%    yP = V(:,501)./1e6;
%    h = plot(xP,yP);
%    set(h,'lineWidth',2);
%    
% % labels
%    tm = 'Potentail in Y direction';
%    h = title(tm); set(h,'fontweight','normal')
%    xlabel('y  [m]'); ylabel('MV/m  [ MV ]');   
%    legend('x = 1.5 m','x = 0.5 m');
%    
% % Graphics parameters
%    set(gca,'fontsize',fs);       
%    set(gca,'xLim',[minY maxY]);
%    
%    
% %% 
% figure(4)   % 44444444444444444444444444444444444444444444444444444444444
%    set(gcf,'units','normalized','position',[0.73 0.52 0.23 0.32]);
%   
%    xP = y; yP = E(:,701)./1e6;
%    h = plot(xP,yP);
%    set(h,'lineWidth',2);
%    hold on
%    yP = E(:,501)./1e6;
%    h = plot(xP,yP);
%    set(h,'lineWidth',2);
%    
% % labels
%    tm = 'Electric Field in Y direction';
%    h = title(tm); set(h,'fontweight','normal')
%    xlabel('y  [m]'); ylabel('MV/m  [ MV/m ]');   
%    legend('x = 1.5 m','x = 0.5 m');
%    
% % Graphics parameters
%    set(gca,'fontsize',fs);       
%    set(gca,'xLim',[minY maxY]);
%      
% 
%    
% %%      
% figure(5)   % 55555555555555555555555555555555555555555555555555555555555
%    set(gcf,'units','normalized','position',[0.01 0.1 0.23 0.32]);        
%    grid on
%    hold on
%    plot(A,V01/1e6,'r','lineWidth',2);
%    plot(A,V02/1e6,'m','lineWidth',2);
%    plot(A(1:10:end),V03(1:10:end)/1e6,'ok');
%    plot(A,Vtotal/1e6,'b','lineWidth',2);
%    
%    h = legend('V_Q','V_{QI}','V_C', 'V_{total}');
%    set(h,'fontSize',12);
%    set(h,'Location','southoutside','orientation','horizontal');
%    
% % labels
%    tm = 'potential in circle around sphere';
%    h = title(tm); set(h,'fontweight','normal')
%    xlabel('angle  [deg]'); ylabel('V  [ MV ]');   
% 
% % Graphics parameters
%    set(gca,'xLim',[-180 180]);
%    set(gca,'xTick',-180:45:180);
%    set(gca,'fontsize',fs);
%    grid on
%    
% %%
% figure(6)   % 66666666666666666666666666666666666666666666666666666666666
%    set(gcf,'units','normalized','position',[0.25 0.10 0.23 0.32]);  
%    contourf(xG,yG,V./1e6,16)
%       
%    xlabel('x  [m]'); ylabel('y  [m]');
%    title('potential  /  electric field lines','fontweight','normal');
%   
%    hold on
%      
% % streamline
%    a = 0.2; ang = linspace(0,360,36);
%    sx = xQ + a .* cosd(ang); sy = yQ + a .* sind(ang);
%    p3 = Ex; p4 = Ey;
%    h = streamline(xG,yG,p3,p4,sx,sy);
%    set(h,'linewidth',2,'color',[1 1 1]);
%       
% % charged conductors 
%    col = [0 0 0];
%    pos = [-RS, -RS, 2*RS, 2*RS];
%    h = rectangle('Position',pos,'Curvature',[1,1]);
%    set(h,'FaceColor',col,'EdgeColor',col);
%        
%    col = [1 0 0]; b = 0.2;
%    pos = [-b+xQ, -b+yQ, 2*b, 2*b];
%    h = rectangle('Position',pos,'Curvature',[1,1]);
%    set(h,'FaceColor',col,'EdgeColor',col);
%         
%    shading interp
%    h = colorbar;
%    h.Label.String = 'V   [ MV ]';
%    h = colormap('parula');
%    beta = 0.5;
%    brighten(h,beta);
%       
%    set(gca,'xLim',[minX,maxX]); set(gca,'yLim', [minY, maxY]);
%    set(gca,'fontsize',fs);
%    axis square;
%      box on;
%      
% %%  
% figure(7)    % 7777777777777777777777777777777777777777777777777777777777
%      set(gcf,'units','normalized','position',[0.49 0.1 0.23 0.32]); 
%      hold on
%      index1 = 51 : 50 : 951;
%      index1 = [index1 500 502];
%      index2 = index1;
%           
%      p1 = xG(index1, index2); p2 = yG(index1, index2);
%      
% % scaling of electric field lines: unit length
%      p3 = Ex(index1, index2)./(E(index1,index2));
%      p4 = Ey(index1, index2)./(E(index1,index2));
% 
% % no scaling of electric field lines
% %  p3 = Ex(index1, index2); p4 = Ey(index1, index2); 
%      
%      h = quiver(p1,p2,p3,p4,'autoscalefactor',0.8);
%      set(h,'color',[0 0 1],'linewidth',1.2)
%      
%      hold on
% 
% % charged conductors 
%     col = [0 0 0];
%     pos = [-RS, -RS, 2*RS, 2*RS];
%     h = rectangle('Position',pos,'Curvature',[1,1]);
%     set(h,'FaceColor',col,'EdgeColor',col);
%       
%     col = [1 0 0]; b = 0.2;
%     pos = [-b+xQ, -b+yQ, 2*b, 2*b];
%     h = rectangle('Position',pos,'Curvature',[1,1]);
%     set(h,'FaceColor',col,'EdgeColor',col);
%       
%     xlabel('x  [m]'); ylabel('y  [m]');
%     title('direction of scaled E at grid points','fontweight','normal');
%      
%     set(gca,'xLim',[minX maxX]); set(gca,'yLim', [minY maxY]);
%     set(gca,'fontsize',fs);
%     axis equal
%     box on   
% 
%    
% %%
% figure(8)   % 88888888888888888888888888888888888888888888888888888888888
%    set(gcf,'units','normalized','position',[0.73 0.10 0.23 0.32]);        
%    plot(A, sigmaC*1e6,'b','lineWidth',2);
%    hold on
%    %plot(A(1:10:end),sigma(1:10:end)*1e6,'ro');
%    
%    xlabel('angle  [deg]'); ylabel('sigma  [\muC/m^2]');
%    title('Sphere: surface charge density','fontweight','normal');
%    
%    set(gca,'xLim',[-180 180]);
%    set(gca,'xTick',-180:45:180);
%    set(gca,'fontsize',fs);
%    grid on
%    
% %%
% figure(9)   % 999999999999999999999999999999999999999999999999999999999999
%    %A = 0 : .01 : 2 * pi;
%    %sigma = -(Q/(4*pi)).*(xQ^2 - RS^2)./(RS*(RS^2+xQ^2-2.*RS.*xQ.*cos(A)).^(3/2));
%    AR = deg2rad(A);
%    h = polar(AR, abs(sigmaC)*1e6, 'b');
%    set(h,'LineWidth',2);
%       
%    title('Sphere: surface charge density','fontweight','normal');
%        
% %%   
toc
% 
% 
% 
% 
