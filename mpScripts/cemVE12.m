% cemVE13.m
% 27 Aapril 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% Potential and Electric Field in [2D] region
% METHOD OF IMAGES:  charge and GROUNDED PLATE

clear all
close all
clc
tic

% ========================================================================
% INPUTS  
% ========================================================================

% Number of grid point    [N = 1001]
   N = 1001;
% 1: Charge  Q   charge / position / radius
   Q = 50e-6;
   xQ = 0; yQ = 0;
   radius = 0.1;

% 2: Image charge     charge / position
   QI = -Q;
   xI = -2; yI = 0;

   
% Dimensions of region 
   minX = -1.018;  
   maxX =  1.018;
   minY = -1.018;
   maxY =  1.018;
%    minR = 1e-4;
%    minRx = 1e-4;
%    minRy = 1e-4;
%    Vsat = 5e6;
%    Esat = 20e6;
 
% =======================================================================
% SETUP 
% =======================================================================

% constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);
% Saturation levels
  Vsat = kC * Q / radius;
  Esat = kC * Q / radius^2;

% [2D] region
   x = linspace(minX,maxX,N);
   y = linspace(minY, maxY,N);
   
% Grid positions
   [xG, yG] = meshgrid(x,y);
   dx = xG(2)-xG(1); dy = yG(2)-yG(1);
        
% =======================================================================
% CALCULATION: POTENTIAL & ELECTRIC FIELD 
% =======================================================================

% 1: potential & electric field for Q      
   RQx = xG - xQ;
   RQy = yG - yQ;
   
   RQ = sqrt(RQx.^2 + RQy.^2);
%   RQ(RQ < minR) = minR;
   VQ = (kC * Q)./ RQ;
   VQ(VQ > Vsat) = Vsat;
   
   RQ3 = RQ.^3;
   EQx = kC .* Q .* RQx ./ RQ3;
   EQy = kC .* Q .* RQy ./ RQ3;
   
   EQ = sqrt(EQx.^2 + EQy.^2);
%    if max(max(EQ)) >=  Esat; EQ(EQ >  Esat)     =  Esat; end
%    if max(max(EQx)) >=  Esat; EQx(EQx >  Esat)  =  Esat; end
%    if min(min(EQx)) <= -Esat; EQx(EQx < -Esat)  = -Esat; end
%    if max(max(EQy)) >=  Esat; EQy(EQy >  Esat)  =  Esat; end
%    if min(min(EQy)) <= -Esat; EQy(EQy < -Esat)  = -Esat; end
%    
% 2: potential & electric field for QI
   RIx = xG - xI;
   RIy = yG - yI;
   
   RI = sqrt(RIx.^2 + RIy.^2);
 %  RI(RI < minR) = minR;
   VI = (kC * QI)./ RI;
   VI(abs(VI) > Vsat) = -Vsat;
   
   RI3 = RI.^3;
   EIx = kC .* QI .* RIx ./ RI3;
   EIy = kC .* QI .* RIy ./ RI3;
   
   EI = sqrt(EIx.^2 + EIy.^2);
%    if max(max(EI)) >=  Esat; EI(EI >  Esat)  =  Esat; end
%    if min(min(EI)) <= -Esat; EI(EI < -Esat)  = -Esat; end
%    
%    if max(max(EIx)) >=  Esat; EIx(EIx >  Esat)  =  Esat; end
%    if min(min(EIx)) <= -Esat; EIx(EIx < -Esat)  = -Esat; end
%    
%    if max(max(EIy)) >=  Esat; EIy(EIy >  Esat)  =  Esat; end
%    if min(min(EIy)) <= -Esat; EIy(EIy < -Esat)  = -Esat; end
   
   R = sqrt(xG.^2 + yG.^2);
   
% total potential and electric field
   V = VQ + VI;
   V(V > Vsat) = Vsat;
 
   Ex = EQx + EIx;
   Ey = EQy + EIy;   
   E = sqrt(Ex.^2 + Ey.^2);
   E(E > Esat) = Esat;  



% ======================================================================= 
% GRAPHICS 
% ======================================================================= 
 fs = 12;
%%
figure(1)   % potential --------------------------------------------------
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
% potential
   zP = (V./1e6);
   contourf(xG,yG,zP,16);
   hold on
% charged sphere 
   col = [1 0 0];
   pos = [-radius, -radius, 2*radius, 2*radius];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
% grounded plate
   xP = [minX minX]; yP = [minY maxY];
   plot(xP,yP,'color',[0 0 0],'linewidth',6)
% titles
   xlabel('x  [m]'); ylabel('y  [m]'); 
   title('potential','fontweight','normal');
% graphics parameters
   shading interp
   h = colorbar;
   h.Label.String = 'V   [ MV ]';
   h = colormap('parula');
   beta = 0.5;
   brighten(h,beta);
   set(gca,'fontsize',fs);
   axis equal

%%
figure(2)   % potential --------------------------------------------------
   set(gcf,'units','normalized','position',[0.01 0.10 0.23 0.32]);

   zP = (V./1e6);
   surf(xG,yG,zP);
   hold on
% charged sphere 
   col = [1 0 0];
   pos = [-radius, -radius, 2*radius, 2*radius];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
% grounded plate
   xP = [minX minX]; yP = [minY maxY];
   plot(xP,yP,'color',[0 0 0],'linewidth',6)   
% titles
   xlabel('x  [m]'); ylabel('y  [m]'); 
   title('potential','fontweight','normal');
% graphics parameters
   shading interp
   h = colorbar;
   h.Label.String = 'V   [ MV ]';
   h = colormap('parula');
   beta = 0.5;
   brighten(h,beta);
   set(gca,'fontsize',fs);
   view(-28,47);
   
   
figure(3)   % field ----------------------------------------------------
   set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
   
   contourf(xG,yG,E./1e6,16);
   hold on
% grounded plate
   xP = [minX minX]; yP = [minY maxY];
   plot(xP,yP,'color',[0 0 0],'linewidth',6)
% charged sphere 
   col = [1 0 0];
   pos = [-radius, -radius, 2*radius, 2*radius];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
   
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('E  [ MV / m ]');
   title('electric field | E |','fontweight','normal');
   
   colorbar
   shading interp
   h =  colorbar;
   h.Label.String = '| E |   [ MV/m ]';
      
   axis square
   set(gca,'fontsize',fs) 
   box on

figure(4)   % field ----------------------------------------------------
   set(gcf,'units','normalized','position',[0.25 0.10 0.23 0.32]);
   
   surf(xG,yG,E./1e6);
   hold on
% charged sphere 
   col = [1 0 0];
   pos = [-radius, -radius, 2*radius, 2*radius];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
% grounded plate
   xP = [minX minX]; yP = [minY maxY];
   plot(xP,yP,'color',[0 0 0],'linewidth',6)
   
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('E  [ MV / m ]');
   title('electric field | E |','fontweight','normal');
   
   colorbar
   shading interp
   h =  colorbar;
   h.Label.String = '| E |   [ MV/m ]';
      
   view(-28,47);
   set(gca,'fontsize',fs) 
   box on      

figure(5)   % potential -----------------------------------------------
   set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
   indX1 = find(x>-0.5,1); indX2 = find(x>0.5,1);
   xP = y; yP = V(:,indX1)./1e6;
   h = plot(xP,yP);
   set(h,'color','r','lineWidth',2);
   hold on
   yP = V(:,indX2)./1e6;
   h = plot(xP,yP);
   set(h,'color','b','lineWidth',2);
   
% labels
   tm = 'Potential  V  in Y direction';
   h = title(tm); set(h,'fontweight','normal')
   xlabel('y  [m]'); ylabel('V    [ MV ]');   
   tm1 = 'x = ';
   tm2 = num2str(x(indX1),2);
   tm3 = num2str(x(indX2),2);
   tm4= '  m';
   tmA = [tm1 tm2 tm4]; 
   tmB = [tm1 tm3 tm4]; 
   legend(tmA,tmB);
   
% Graphics parameters
   set(gca,'fontsize',fs);       
   set(gca,'xLim',[minY maxY]);
   grid on
   
figure(6)   % field ----------------------------------------------------
   set(gcf,'units','normalized','position',[0.49 0.10 0.23 0.32]);
   indX1 = find(x>-0.5,1); indX2 = find(x>0.5,1);
   xP = y; yP = E(:,indX1)./1e6;
   h = plot(xP,yP);
   set(h,'color','r','lineWidth',2);
   hold on
   yP = E(:,indX2)./1e6;
   h = plot(xP,yP);
   set(h,'color','b','lineWidth',2);
   
% labels
   tm = 'Electric Field | E | in Y direction';
   h = title(tm); set(h,'fontweight','normal')
   xlabel('y  [m]'); ylabel('| E |    [ MV/m ]');   
   tm1 = 'x = ';
   tm2 = num2str(x(indX1),2);
   tm3 = num2str(x(indX2),2);
   tm4= '  m';
   tmA = [tm1 tm2 tm4]; 
   tmB = [tm1 tm3 tm4]; 
   legend(tmA,tmB);
   
% Graphics parameters
   set(gca,'fontsize',fs);       
   set(gca,'xLim',[minY maxY]);
   grid on 

   
figure(7)   % ----------------------------------------------------------
   set(gcf,'units','normalized','position',[0.75 0.52 0.23 0.32]);  
   contourf(xG,yG,V./1e6,16)
% charged sphere 
   col = [1 0 0];
   pos = [-radius, -radius, 2*radius, 2*radius];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);   
   xlabel('x  [m]'); ylabel('y  [m]');
   title('potential  /  electric field lines','fontweight','normal');
  
   hold on
     
% streamline
   a = 0.2; ang = linspace(0,360,36);
   sx = xQ + a .* cosd(ang); sy = yQ + a .* sind(ang);
   p3 = Ex; p4 = Ey;
   h = streamline(xG,yG,p3,p4,sx,sy);
   set(h,'linewidth',2,'color',[1 1 1]);

% grounded plate
   xP = [minX minX]; yP = [minY maxY];
   plot(xP,yP,'color',[0 0 0],'linewidth',6)
        
   shading interp
   h = colorbar;
   h.Label.String = 'V   [ MV ]';
   h = colormap('parula');
   beta = 0.5;
   brighten(h,beta);
      
   set(gca,'xLim',[minX,maxX]); set(gca,'yLim', [minY, maxY]);
   set(gca,'fontsize',fs);
   axis square;
     box on;
     
 %%
figure(8)    % ---------------------------------------------------------
     set(gcf,'units','normalized','position',[0.75 0.1 0.23 0.32]); 
     hold on
% grounded plate
   xP = [minX minX]; yP = [minY maxY];
   plot(xP,yP,'color',[0.5 0.5 0.5],'linewidth',6)     
     index1 = 51 : 100 : 951;
    % index1 = [index1 500 502];
     index2 = index1;
          
     p1 = xG(index1, index2); p2 = yG(index1, index2);
     
% scaling of electric field lines: unit length
     p3 = Ex(index1, index2)./(E(index1,index2));
     p4 = Ey(index1, index2)./(E(index1,index2));

% no scaling of electric field lines
%  p3 = Ex(index1, index2); p4 = Ey(index1, index2); 
     
     h = quiver(p1,p2,p3,p4,'autoscalefactor',0.5);
     set(h,'color',[0 0 1],'linewidth',1.2)
     
     hold on 
    
% charged sphere 
   col = [1 0 0];
   pos = [-radius, -radius, 2*radius, 2*radius];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);

     
    xlabel('x  [m]'); ylabel('y  [m]');
    title('direction of scaled E at grid points','fontweight','normal');
     
    set(gca,'xLim',[minX maxX]); set(gca,'yLim', [minY maxY]);
    set(gca,'fontsize',fs);
    axis equal
    box on   

   

   toc




