% cemVE07.m
% 18 april 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% Potential and Electric Field in [2D] region
% [2D] quadrupole: 4 charges, equal magnitude, alternating signs at the 
%                  corners of a square                       



clear all
close all
clc
tic



% INPUTS  ================================================================

% Number of grid point    [N = 10001]
   N = 1001;
   
% Charge  Q = [10, 0, 0, 0, 0] 
   Q = [20, -20, 20, -20, 0] .* 1e-6;
   
% Radius of circular charged conductor;   
   a = 0.1;
   
% X & Y components of position of charges  [0, 0, 0, 0, 0]
   xC = [0.5,  0 ,  -0.5,   0 ,   0];
   yC = [ 0 , 0.5,   0  , -0.5, 0];

% 5 random charges   uncomment to run the program for 5 random charges
%   Q = (1 + 9 .* rand(5,1)) .* 1e-6;
%   xC = -2 + 4 .* rand(5,1); 
%   yC = -2 + 4 .* rand(5,1); 

% constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);
   
% Dimensions of region / saturation levels
%   [dimensions of region -2 to 2 / minR = 1e-6 / Esat = 1e6 / Vsat = 1e6]
   minX = -2;  
   maxX =  2;
   minY = -2;
   maxY =  2;
   minR = 1e-6;
   minRx = 1e-6;
   minRy = 1e-6;
   Vsat = kC * max(abs(Q)) / a;
   Esat = kC * max(abs(Q)) / a^2;
   
% SETUP  =================================================================
   
  % fields
    V = zeros(N,N);
    Ex = zeros(N,N); Ey = zeros(N,N);
  % [2D] region
    x  = linspace(minX,maxX,N);
    y = linspace(minY, maxY,N);
    
  % color of charged object  +  red   /   - black
    col1 = [1 0 0];
    col2 = [0 0 0];
    %if Q(1) < 0; col1 = [0 0 0]; end;
        
  % grid positions
    [xG, yG] = meshgrid(x,y);


% CALCULATION: POTENTIAL & ELECTRIC FIELD ================================

for n = 1 : 4
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
   
   if max(max(Ex)) >=  Esat; Ex(Ex >  Esat)  =  Esat; end;
   if min(min(Ex)) <= -Esat; Ex(Ex < -Esat)  = -Esat; end;
   
   if max(max(Ey)) >=  Esat; Ey(Ey >  Esat)  =  Esat; end;
   if min(min(Ey)) <= -Esat; Ey(Ey < -Esat)  = -Esat; end;
   

% GRAPHICS ===============================================================

%%

figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
   surf(xG,yG,V./1e6);
   
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('V  [ V ]');
   title('potential','fontweight','normal');
   
   rotate3d
   view(76,32);
   set(gca,'fontsize',12)
   set(gca,'xLim',[-2, 2]); set(gca,'yLim',[-2, 2]);
   
   shading interp;
   h = colorbar;
   h.Label.String = 'V   [ MV ]';
   colormap(parula);
   axis square
   box on

%%   
figure(2)   %2222222222222222222222222222222222222222222222222222222222222
   set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
   zP = V./1e6;
   contourf(xG,yG,zP,16);
   %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
   %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
   hold on
   
   % charged conductors 
      col = col1;
      if Q(1) < 0; col = col2; end;
      pos1 = [-a+xC(1), -a+yC(1), 2*a, 2*a];
      h = rectangle('Position',pos1,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(2) < 0; col = col2; end;
      pos2 = [-a-xC(2), -a+yC(2), 2*a, 2*a];
      h = rectangle('Position',pos2,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(3) < 0; col = col2; end;
      pos3 = [-a+xC(3), -a+yC(3), 2*a, 2*a];
      h = rectangle('Position',pos3,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
   
   xlabel('x  [m]'); ylabel('y  [m]');
   title('potential','fontweight','normal');
      
   shading interp
   h = colorbar;
   h.Label.String = 'V   [ MV ]';
   colormap(parula);
   
   
   set(gca,'fontsize',12);
   axis square
   box on

%%

figure(3)   % 33333333333333333333333333333333333333333333333333333333333
   set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
       c = 0; yStep = zeros(1,5);
       
      for n = ceil(N/2) : ceil(N/10): N
   %  for n = [1 11 21 31 41 51]
        xP = xG(n,:); yP = V(n,:)./1e6;
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
       xlabel('x  [m]'); ylabel('V  [ MV ]');
       set(gca,'fontsize',12)  
       h_title = title(tm);
       set(h_title,'fontsize',12,'FontWeight','normal');
       %set(gca,'xTick',-5:5);
       tm1 = num2str(yStep(1),2);
       tm2 = num2str(yStep(2),2); tm3 = num2str(yStep(3),2);
       tm4 = num2str(yStep(4),2);  tm5 = num2str(yStep(5),2);
       %tm6 = num2str(yStep(6),3);
       h = legend(tm1,tm2,tm3,tm4,tm5,'Orientation','horizontal');
       set(h,'Location','northOutside');
      % 
%%
figure(4)   % 4444444444444444444444444444444444444444444444444444444444
     set(gcf,'units','normalized','position',[0.73 0.52 0.23 0.32]);
     
     c = 0; yStep = zeros(1,5);
     for n = ceil(N/2) : ceil(N/10): N
   %  for n = [1 11 21 31 41 51]
        xP = xG(n,:); yP = E(n,:)./1e6;
        plot(xP,yP,'linewidth',2);
        hold on
        c = c+1;
        yStep(c) = yG(n,1);
     end
     
    xlabel('x  [m]'); ylabel('| E |  [ MV/m ]');
    tm = '| E | profiles for different y values';
    h = title(tm);
    set(h,'fontsize',12,'FontWeight','normal');
       
    tm1 = num2str(yStep(1),2);
    tm2 = num2str(yStep(2),2); tm3 = num2str(yStep(3),2);
    tm4 = num2str(yStep(4),2);  tm5 = num2str(yStep(5),2);
    h = legend(tm1,tm2,tm3,tm4,tm5);
    set(h,'Orientation','horizontal','Location','northOutside'); 
       
    set(gca,'fontsize',12);
    box on;

%%    
figure(5)   % 5555555555555555555555555555555555555555555555555555555555
set(gcf,'units','normalized','position',[0.01 0.1 0.23 0.32]);
   surf(xG,yG,E./1e6);
   
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('E  [ V / m ]');
   title('electric field | E |','fontweight','normal');
   
  
   colorbar
   shading interp
   h =  colorbar;
   h.Label.String = '| E |   [ MV/m ]';
   rotate3d
   view(30,50);
   
   axis square
   set(gca,'fontsize',12) 
   box on
   
%%
figure(6)   % 66666666666666666666666666666666666666666666666666666666666
   set(gcf,'units','normalized','position',[0.25 0.1 0.23 0.32]);  
   contourf(xG,yG,V./1e6,12)
      
   xlabel('x  [m]'); ylabel('y  [m]');
   title('potential','fontweight','normal');
  
   hold on
     
   % Electric field lines
   %     need to reverse the direction of streamlines for neg charges
   Ns = abs(Q(2)/Q(1));       % ratio of charges
   A = pi/6;   
   
   % Charge 1 
      %ang = pi/3  :pi/3 : 2*pi-pi/3;      % same charge
      %ang = [pi/3 2*pi/6 3*pi/6 4*pi/6 pi/6 0*pi/6 7*pi/6 8*pi/6 9*pi/6  11*pi/6];  
      
      %ang1 = [0*A 1*A 2*A 3*A 4*A 5*A 6*A 7*A 8*A 9*A  ];
      ang1 = linspace(5,355,20); 
      ang1 = deg2rad(ang1);
      if Q(1) > 0;
          p3 = Ex; p4 = Ey;
      else 
          p3 = -Ex; p4 = -Ey;
          
      end
      sx = xC(1) + a .* cos(ang1); sy = yC(1) + a .* sin(ang1);
      %sx(sy==0) = sx(1); sy(sy==0) = sy(1);
      h = streamline(xG,yG,p3,p4,sx,sy);
      set(h,'linewidth',2,'color',[1 1 1]);
          
   % Charge 2 
      %ang = [-A 0*A 1*A 3*A 4*A 5*A 6*A 7*A 8*A 9*A   ];
      ang2 = ang1 + pi/2;
      if Q(2) > 0;
          p3 = Ex; p4 = Ey;
      else 
          p3 = -Ex; p4 = -Ey;
      end
      
      sx = xC(2) + a .* cos(ang2); sy = yC(2) + a .* sin(ang2);
      h = streamline(xG,yG,p3,p4,sx,sy);
      set(h,'linewidth',2,'color',[1 1 1]);
      
   % Charge 3 
      
      if Q(3) > 0;
          p3 = Ex; p4 = Ey;
      else 
          p3 = -Ex; p4 = -Ey;
      end
      
      %ang3 = [1*A 2*A 3*A 4*A 5*A 7*A 8*A 9*A 10*A 11*A  ];
      ang3 = ang1 + pi;
      sx = xC(3) + a .* cos(ang3); sy = yC(3) + a .* sin(ang3);
      h = streamline(xG,yG,p3,p4,sx,sy);
      set(h,'linewidth',2,'color',[1 1 1]);
      
     % Charge 4 
      
      if Q(4) > 0;
          p3 = Ex; p4 = Ey;
      else 
          p3 = -Ex; p4 = -Ey;
      end
      
      %ang3 = [1*A 2*A 3*A 4*A 5*A 7*A 8*A 9*A 10*A 11*A  ];
      ang4 = ang1 - pi/2;
      sx = xC(4) + a .* cos(ang4); sy = yC(4) + a .* sin(ang4);
      h = streamline(xG,yG,p3,p4,sx,sy);
      set(h,'linewidth',2,'color',[1 1 1]);  
      
      % charged conductors 
      col = col1;
      if Q(1) < 0; col = col2; end;
      pos1 = [-a+xC(1), -a+yC(1), 2*a, 2*a];
      h = rectangle('Position',pos1,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(2) < 0; col = col2; end;
      pos2 = [-a-xC(2), -a+yC(2), 2*a, 2*a];
      h = rectangle('Position',pos2,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(3) < 0; col = col2; end;
      pos3 = [-a+xC(3), -a+yC(3), 2*a, 2*a];
      h = rectangle('Position',pos3,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(4) < 0; col = col2; end;
      pos4 = [-a+xC(4), -a+yC(4), 2*a, 2*a];
      h = rectangle('Position',pos4,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
        
   shading interp
   h = colorbar;
   h.Label.String = 'V   [ MV ]';
   h = colormap('parula');
   beta = 0.5;
   brighten(h,beta);
      
   set(gca,'xLim',[minX,maxX]); set(gca,'yLim', [minY, maxY]);
  
   axis square;
   box on;
     
%%
figure(7)   % 7777777777777777777777777777777777777777777777777777777777
   set(gcf,'units','normalized','position',[0.49 0.1 0.23 0.32]);
   pcolor(xG,yG,E./1e6);
   
   xlabel('x  [m]'); ylabel('y  [m]');
   title('electric field | E | ','fontweight','normal');
   
   %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
   %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
   
   shading interp
   h = colorbar;
   h.Label.String = '|E|   [ MV/m ]';
   
   axis square
   set(gca,'fontsize',12)
   box on   
   

 
%%       
figure(8)    % 8888888888888888888888888888888888888888888888888888888888
     set(gcf,'units','normalized','position',[0.73 0.1 0.23 0.32]); 
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
      col = col1;
      if Q(1) < 0; col = col2; end;
      pos1 = [-a+xC(1), -a+yC(1), 2*a, 2*a];
      h = rectangle('Position',pos1,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(2) < 0; col = col2; end;
      pos2 = [-a-xC(2), -a+yC(2), 2*a, 2*a];
      h = rectangle('Position',pos2,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(3) < 0; col = col2; end;
      pos3 = [-a+xC(3), -a+yC(3), 2*a, 2*a];
      h = rectangle('Position',pos3,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
      
      col = col1;
      if Q(4) < 0; col = col2; end;
      pos4 = [-a+xC(4), -a+yC(4), 2*a, 2*a];
      h = rectangle('Position',pos4,'Curvature',[1,1]);
      set(h,'FaceColor',col,'EdgeColor',col);
     
     xlabel('x  [m]'); ylabel('y  [m]');
     title('direction of scaled E at grid points','fontweight','normal');
     
     set(gca,'xLim',[-2,2]); set(gca,'yLim', [-2, 2]);
     axis([-2 2 -2 2]);
     axis equal
     box on
 
     %%

     toc
     
   