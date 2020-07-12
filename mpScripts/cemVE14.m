% cemVE14.m
% 28 April 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Potential and Electric Field in [2D] region
% METHOD OF IMAGES:  charge and INSULATED conductive sphere

clear all
close all
clc
tic

% ========================================================================
% INPUTS  
% ========================================================================

% Number of grid point    [N = 1001]
   N = 1001;
% 1: Charge  Q   charge / position
   Q = 200e-6;
   xQ = 1.2; yQ = 0;
   
% 2: Conducting Sphere     radius / position
   RS = 0.25;
   xS = 0; yS = 0;
   
% Dimensions of region / saturation levels
   minX = -2;  
   maxX =  3;
   minY = -2.5;
   maxY =  2.5;
   minR = 1e-3;
   minRx = 1e-6;
   minRy = 1e-6;
   Vsat = 10e6;
   Esat = 10e6;
   
% =======================================================================
% SETUP 
% =======================================================================

% constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);

% Image charge     charge / position
   QI = -(RS/xQ) * Q;
   xI = RS^2 / xQ; yI = 0;

% [2D] region
   x  = linspace(minX,maxX,N);
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
   RQ(RQ < minR) = minR;
   VQ = (kC * Q)./ RQ;
   VQ(VQ > Vsat) = Vsat;
   
   RQ3 = RQ.^3;
   EQx = kC .* Q .* RQx ./ RQ3;
   EQy = kC .* Q .* RQy ./ RQ3;
   
   EQ = sqrt(EQx.^2 + EQy.^2);
   if max(max(EQ)) >=  Esat; EQ(EQ >  Esat)     =  Esat; end;
   if max(max(EQx)) >=  Esat; EQx(EQx >  Esat)  =  Esat; end;
   if min(min(EQx)) <= -Esat; EQx(EQx < -Esat)  = -Esat; end;
   if max(max(EQy)) >=  Esat; EQy(EQy >  Esat)  =  Esat; end;
   if min(min(EQy)) <= -Esat; EQy(EQy < -Esat)  = -Esat; end;
   
% 2: potential & electric field for QI
   RIx = xG - xI;
   RIy = yG - yI;
   
   RI = sqrt(RIx.^2 + RIy.^2);
   RI(RI < minR) = minR;
   VI = (kC * QI)./ RI;
   VI(abs(VI) > Vsat) = -Vsat;
   
   RI3 = RI.^3;
   EIx = kC .* QI .* RIx ./ RI3;
   EIy = kC .* QI .* RIy ./ RI3;
   
   EI = sqrt(EIx.^2 + EIy.^2);
   if max(max(EI)) >=  Esat; EI(EI >  Esat)  =  Esat; end;
   if min(min(EI)) <= -Esat; EI(EI < -Esat)  = -Esat; end;
   
   if max(max(EIx)) >=  Esat; EIx(EIx >  Esat)  =  Esat; end;
   if min(min(EIx)) <= -Esat; EIx(EIx < -Esat)  = -Esat; end;
   
   if max(max(EIy)) >=  Esat; EIy(EIy >  Esat)  =  Esat; end;
   if min(min(EIy)) <= -Esat; EIy(EIy < -Esat)  = -Esat; end;
   
   R = sqrt(xG.^2 + yG.^2);

% 3: potential & electric field for -QI at centre of sphere (0,0)   C
   QC = -QI;
   RCx = xG;
   RCy = yG;
   
   RC = sqrt(RCx.^2 + RCy.^2);
   RC(RC < minR) = minR;
   VC = (kC * QC)./ RC;
   VC(abs(VC) > Vsat) = Vsat;
   
   RC3 = RC.^3;
   ECx = kC .* QC .* RCx ./ RC3;
   ECy = kC .* QC .* RCy ./ RC3;
   
   EC = sqrt(ECx.^2 + ECy.^2);
   if max(max(EC)) >=  Esat; EC(EC >  Esat)  =  Esat; end;
   if min(min(EC)) <= -Esat; EC(EC < -Esat)  = -Esat; end;
   
   if max(max(ECx)) >=  Esat; ECx(ECx >  Esat)  =  Esat; end;
   if min(min(ECx)) <= -Esat; ECx(ECx < -Esat)  = -Esat; end;
   
   if max(max(ECy)) >=  Esat; ECy(ECy >  Esat)  =  Esat; end;
   if min(min(ECy)) <= -Esat; ECy(ECy < -Esat)  = -Esat; end;
   
   
   
% total potential and electric field
   V = VQ + VI + VC;
   V(R <= RS) = 0;
   V(V > Vsat) = Vsat;
 
   Ex = EQx + EIx + ECx;
   Ey = EQy + EIy + ECy;   
   E = sqrt(Ex.^2 + Ey.^2);
   Ex(R <= RS) = 0;
   Ey(R <= RS) = 0;
   E(R <= RS)  = 0;
  
%%
% Charge on outer surface of spherical conductor -------------------------
   A = linspace(-180,180,101);    % MUST HAVE AN ODD NUMBER OF ELEMENTS
   %A = 79.5;
   d = RS + 0.00;
   x2 = d .* cosd(A);
   y2 = d .* sind(A);
   R2x = x2 - xI; R2y = y2 - yI;
   R2 = sqrt(R2x.^2 + R2y.^2);
   R23 = R2.^3;
   Ex2 = kC .* QI .* R2x ./ R23;
   Ey2 = kC .* QI .* R2y ./ R23;
   V02 = kC .* QI ./ R2;
   
   R3x = x2; R3y = y2;
   R3 = sqrt(R3x.^2 + R3y.^2);
   R33 = R3.^3;
   Ex3 = kC .* QC .* R3x ./ R33;
   Ey3 = kC .* QC .* R3y ./ R33;
   V03 = kC .* QC ./ R3;

   R1x = x2 - xQ; R1y = y2 - yQ;
   R1 = sqrt(R1x.^2 + R1y.^2);
   R13 = R1.^3;
   Ex1 = kC .* Q .* R1x ./ R13;
   Ey1 = kC .* Q .* R1y ./ R13;
   V01 = kC .* Q ./ R1;

   Vtotal = V01 + V02 + V03;
   
   Ex12 = Ex1 + Ex2 + Ex3; Ey12 = Ey1 + Ey2 + Ey3;
   Etotal = sqrt(Ex12.^2 + Ey12.^2);
   sigmaC = Etotal .* eps0;
   
   angle = atan2(Ey12,Ex12);
   angle = rad2deg(angle);
   
   % calculate angle of radial electric field and set sign for sigma
   for n = 2 : 51
      if angle(n) > 0; sigmaC(n) = -sigmaC(n); end;
   end 
   
   for n = 51 : 101
      if angle(n) < 0; sigmaC(n) = -sigmaC(n); end;
   end 
   
   K = -2*pi*Q/(4*pi);
   T = deg2rad(A);
   fn = sigmaC .* RS^2 .* sin(T);
% Find total surface charge on sphere   
   Qinduced = K * simpson1d(fn,0,pi);
   Qinduced
   
   
 % Force between sphere and charge
    F = kC * Q * QI /(xQ - xI)^2 + kC * Q * QC / xQ^2; 
    
%%
% ======================================================================= 
% GRAPHICS 
% ======================================================================= 
   fs = 12;
%%
figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);

   zP = V./1e6;
   contourf(xG,yG,zP,16);
   
% titles
   xlabel('x  [m]'); ylabel('y  [m]'); 
   title('potential','fontweight','normal');
   
% charged conductors 
   col = [0 0 0];
   pos = [-RS, -RS, 2*RS, 2*RS];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
      
   col = [1 0 0]; b = 0.2;
   pos = [-b+xQ, -b+yQ, 2*b, 2*b];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
      
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
figure(2)   % 22222222222222222222222222222222222222222222222222222222222
   set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
   
   contourf(xG,yG,E./1e6,16);
   
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('E  [ V / m ]');
   title('electric field | E |','fontweight','normal');
   
   colorbar
   shading interp
   h =  colorbar;
   h.Label.String = '| E |   [ MV/m ]';
   
% charged conductors 
   col = [0 0 0];
   pos = [-RS, -RS, 2*RS, 2*RS];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
      
   col = [1 0 0]; b = 0.2;
   pos = [-b+xQ, -b+yQ, 2*b, 2*b];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
   
   axis square
   set(gca,'fontsize',fs) 
   box on

      
%% 
figure(3)   % 33333333333333333333333333333333333333333333333333333333333
   set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
  
   xP = y; yP = V(:,701)./1e6;
   h = plot(xP,yP);
   set(h,'lineWidth',2);
   hold on
   yP = V(:,501)./1e6;
   h = plot(xP,yP);
   set(h,'lineWidth',2);
   
% labels
   tm = 'Potentail in Y direction';
   h = title(tm); set(h,'fontweight','normal')
   xlabel('y  [m]'); ylabel('MV/m  [ MV ]');   
   legend('x = 1.5 m','x = 0.5 m');
   
% Graphics parameters
   set(gca,'fontsize',fs);       
   set(gca,'xLim',[minY maxY]);
   
   
%% 
figure(4)   % 44444444444444444444444444444444444444444444444444444444444
   set(gcf,'units','normalized','position',[0.73 0.52 0.23 0.32]);
  
   xP = y; yP = E(:,701)./1e6;
   h = plot(xP,yP);
   set(h,'lineWidth',2);
   hold on
   yP = E(:,501)./1e6;
   h = plot(xP,yP);
   set(h,'lineWidth',2);
   
% labels
   tm = 'Electric Field in Y direction';
   h = title(tm); set(h,'fontweight','normal')
   xlabel('y  [m]'); ylabel('MV/m  [ MV/m ]');   
   legend('x = 1.5 m','x = 0.5 m');
   
% Graphics parameters
   set(gca,'fontsize',fs);       
   set(gca,'xLim',[minY maxY]);
     

   
%%      
figure(5)   % 55555555555555555555555555555555555555555555555555555555555
   set(gcf,'units','normalized','position',[0.01 0.1 0.23 0.32]);        
   grid on
   hold on
   plot(A,V01/1e6,'r','lineWidth',2);
   plot(A,V02/1e6,'m','lineWidth',2);
   plot(A(1:10:end),V03(1:10:end)/1e6,'ok');
   plot(A,Vtotal/1e6,'b','lineWidth',2);
   
   h = legend('V_Q','V_{QI}','V_C', 'V_{total}');
   set(h,'fontSize',12);
   set(h,'Location','southoutside','orientation','horizontal');
   
% labels
   tm = 'potential in circle around sphere';
   h = title(tm); set(h,'fontweight','normal')
   xlabel('angle  [deg]'); ylabel('V  [ MV ]');   

% Graphics parameters
   set(gca,'xLim',[-180 180]);
   set(gca,'xTick',-180:45:180);
   set(gca,'fontsize',fs);
   grid on
   
%%
figure(6)   % 66666666666666666666666666666666666666666666666666666666666
   set(gcf,'units','normalized','position',[0.25 0.10 0.23 0.32]);  
   contourf(xG,yG,V./1e6,16)
      
   xlabel('x  [m]'); ylabel('y  [m]');
   title('potential  /  electric field lines','fontweight','normal');
  
   hold on
     
% streamline
   a = 0.2; ang = linspace(0,360,36);
   sx = xQ + a .* cosd(ang); sy = yQ + a .* sind(ang);
   p3 = Ex; p4 = Ey;
   h = streamline(xG,yG,p3,p4,sx,sy);
   set(h,'linewidth',2,'color',[1 1 1]);
      
% charged conductors 
   col = [0 0 0];
   pos = [-RS, -RS, 2*RS, 2*RS];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
       
   col = [1 0 0]; b = 0.2;
   pos = [-b+xQ, -b+yQ, 2*b, 2*b];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
        
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
figure(7)    % 7777777777777777777777777777777777777777777777777777777777
     set(gcf,'units','normalized','position',[0.49 0.1 0.23 0.32]); 
     hold on
     index1 = 51 : 50 : 951;
     index1 = [index1 500 502];
     index2 = index1;
          
     p1 = xG(index1, index2); p2 = yG(index1, index2);
     
% scaling of electric field lines: unit length
     p3 = Ex(index1, index2)./(E(index1,index2));
     p4 = Ey(index1, index2)./(E(index1,index2));

% no scaling of electric field lines
%  p3 = Ex(index1, index2); p4 = Ey(index1, index2); 
     
     h = quiver(p1,p2,p3,p4,'autoscalefactor',0.8);
     set(h,'color',[0 0 1],'linewidth',1.2)
     
     hold on

% charged conductors 
    col = [0 0 0];
    pos = [-RS, -RS, 2*RS, 2*RS];
    h = rectangle('Position',pos,'Curvature',[1,1]);
    set(h,'FaceColor',col,'EdgeColor',col);
      
    col = [1 0 0]; b = 0.2;
    pos = [-b+xQ, -b+yQ, 2*b, 2*b];
    h = rectangle('Position',pos,'Curvature',[1,1]);
    set(h,'FaceColor',col,'EdgeColor',col);
      
    xlabel('x  [m]'); ylabel('y  [m]');
    title('direction of scaled E at grid points','fontweight','normal');
     
    set(gca,'xLim',[minX maxX]); set(gca,'yLim', [minY maxY]);
    set(gca,'fontsize',fs);
    axis equal
    box on   

   
%%
figure(8)   % 88888888888888888888888888888888888888888888888888888888888
   set(gcf,'units','normalized','position',[0.73 0.10 0.23 0.32]);        
   plot(A, sigmaC*1e6,'b','lineWidth',2);
   hold on
   %plot(A(1:10:end),sigma(1:10:end)*1e6,'ro');
   
   xlabel('angle  [deg]'); ylabel('sigma  [\muC/m^2]');
   title('Sphere: surface charge density','fontweight','normal');
   
   set(gca,'xLim',[-180 180]);
   set(gca,'xTick',-180:45:180);
   set(gca,'fontsize',fs);
   grid on
   
%%
figure(9)   % 999999999999999999999999999999999999999999999999999999999999
   %A = 0 : .01 : 2 * pi;
   %sigma = -(Q/(4*pi)).*(xQ^2 - RS^2)./(RS*(RS^2+xQ^2-2.*RS.*xQ.*cos(A)).^(3/2));
   AR = deg2rad(A);
   h = polar(AR, abs(sigmaC)*1e6, 'b');
   set(h,'LineWidth',2);
      
   title('Sphere: surface charge density','fontweight','normal');
       
%%   
   toc




