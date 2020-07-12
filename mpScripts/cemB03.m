% cemB03.m
% 02 June 2016
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Magnetic field from a set of parallel wires
%    Biot-Savart Law

% SI units used unless stated otherwise

clear all
close all
clc

tic


% ========================================================================
% INPUTS  
% ========================================================================

% Wires: current I / xy (x,y) coordinates / Z length / Z partitions
   Nwire = 1;
   I = [10,-10,10,10];
   xW = [0.000, -0.005, -0.005,  0.005];
   yW = [0.000, -0.005,  0.005, -0.005];
   minZ = -5; maxZ = 5;
   nZ = 4101;
   radiusW = 0.001;
   
% Grid dimensions: Grid points / X range / y range
   N = 401;
   L = 0.005*2;
   minX = -L; maxX = L; 
   minY = -L; maxY = L;
   minR = 1e-3;
   
% Calculation of line integral and magentic flux  L < maxX or maxY
   L = 5e-3;
   minLx = -L; maxLx = L;
   minLy = -L; maxLy = L;
   
   
% Calculation of the divergence at the point a grid closest to the values xD, yD
   xD = 0.05e-3;
   yD = 0.05e-3;
   
   
% =======================================================================
% SETUP 
% =======================================================================

% constants   permeability of free space
   mu0 = 4*pi*1e-7;
   K = I .* (mu0/(4*pi));
   
   %Bsat = mu0 * max(I) / minR;
   Bsat = 2e-3;
   
   zW = linspace(minZ, maxZ, nZ);
   dLz = zW(2)-zW(1);
   dLx = 0; dLy = 0;
   dL = [dLx dLy dLz];
   
 % [2D] region
   x  = linspace(minX,maxX,N);
   y = linspace(minY, maxY,N);
   dx = x(2)-x(1); dy = y(2)-y(1);
   
% Grid positions
   [xG, yG] = meshgrid(x,y);
   dx = x(2)-x(1); dy = y(2)-y(1);
   
   RG = sqrt(xG.^2 + yG.^2);

  for n1 = 1 : Nwire
     for n2 = 1 : N
       if abs(x(n2) - xW(n1)) < dx/2 ; xW(n1) = x(n2) + dx/2; end
       if abs(y(n2) - yW(n1)) < dy/2 ; yW(n1) = y(n2) + dy/2; end
     end
  end
   
   
% =======================================================================
% CALCULATIONS: MAGNETIC FIELD 
% =======================================================================

Bx = zeros(N,N); By = zeros(N,N); Bz = zeros(N,N);

for n = 1 : Nwire
    for m = 1 : nZ
      Rx = xG - xW(n);
      Ry = yG - yW(n);
      Rz = 0  - zW(m);
    %  Rx(Rx<minR) = minR;
    %  Ry(Ry<minR) = minR;
      R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
     % R(R<minR) = minR;
      R3 = R.^3;
      
      Bx = Bx + K(n) .* (dLy .* Rz - dLz .* Ry) ./ R3; 
      By = By + K(n) .* (dLz .* Rx - dLx .* Rz) ./ R3;
      Bz = Bz + K(n) .* (dLx .* Ry - dLy .* Rx) ./ R3;
    end
end
  B = sqrt(Bx.^2 + By.^2 + Bz.^2);
  
%  B(B > Bsat) = Bsat;

% Line integral  B.dL ----------------------------------------------------
%       fn must have an odd number of elements
  X1 = find(x>minLx,1);X2 = find(x>maxLx,1);
  Y1 = find(y>minLy,1);Y2 = find(y>maxLy,1);
  if mod(X2-X1,2) ~= 0; X2 = X2+1; end;
  if mod(Y2-Y1,2) ~= 0; Y2 = Y2+1; end;
          
  a = x(X1); b = x(X2); fn = Bx(Y1,X1:X2); 
  BdL(1) = simpson1d(fn,a,b);
  a = y(Y1); b = y(Y2); fn = By(Y1:Y2,X2)'; 
  BdL(2) = simpson1d(fn,a,b);
  a = x(X1); b = x(X2); fn = -Bx(Y2,X1:X2); 
  BdL(3) = simpson1d(fn,a,b);
  a = y(Y1); b = y(Y2); fn = By(Y1:Y2,X1)'; 
  BdL(4) = -simpson1d(fn,a,b);
  
  Ienclosed = sum(BdL)/mu0;
  
% Gauss's Law for Magentisim:  Magnetic flux through a closed surface
    a = y(Y1); b = y(Y2); fn = Bx(Y1:Y2,X2)';
      phi(1) = simpson1d(fn,a,b);
    a = y(Y1); b = y(Y2); fn = -Bx(Y1:Y2,X1)';
      phi(2) = simpson1d(fn,a,b);
    a = x(X1); b = x(X2); fn = By(Y2,X1:X2);
      phi(3) = simpson1d(fn,a,b);
    a = x(X1); b = x(X2); fn = -Bx(Y1,X1:X2);
      phi(4) = simpson1d(fn,a,b);
  
    phi_total = sum(phi);
   
% Divergence  divB and Curl  curlB  -------------------------------------
  Nx = find(x > xD,1);
  Ny = find(y > yD,1);
  
  divB = divergence(xG,yG,Bx,By);
  min_divB = min(min(divB));
  max_divB = max(max(divB));

  curlB = curl(xG,yG,Bx,By);
  min_curlB = min(min(curlB));
  max_curlB = max(max(curlB));

  dBdx = (Bx(Ny,Nx+1) - Bx(Ny,Nx-1))/(x(Nx+1)-x(Nx-1));
  dBdy = (By(Ny+1,Nx) - By(Ny-1,Nx))/(y(Ny+1)-y(Ny-1));
  divBp = dBdx+dBdy;
  
% curl of B at the position of wire 1
   curlB_1 = mu0*I(1)/(dx/2)^2;
  
% ANALTYICAL CALC. Magentic field surrounding a single infinite
%                  straight conductor with current I(1) 
   N_r = 501;
   min_r = 0.001;
   max_r = maxX;
   r = linspace(min_r,max_r, N_r);
   B_A = mu0*I(1)./(2.*pi.*r);  

   
% ========================================================================
% COMMAND WINDOW OUTPUT
% ========================================================================
disp('   ');
fprintf('Number of wires =  %2.0f   \n',Nwire);
disp('   ');
disp('x location of wires [mm] ');
fprintf('   %2.2f   ',1e3*xW(1:Nwire));
disp('   ');
disp('y location of wires [mm] ');
fprintf('   %2.2f   ',1e3*yW(1:Nwire));
disp('    ');
disp('Wire currents  [A] ');
fprintf('   %2.3f   \n',I(1:Nwire));
disp('    ');
fprintf('   I_total =  %2.3f A  \n',sum(I(1:Nwire)));
disp('   ');
disp('Circulation path and closed surface: line integral & flux');
disp('  ');
fprintf('   I_enclosed  =  %2.2f  A \n',Ienclosed);
disp('   ');
disp('Gauss Law for magnetism: total flux through the closed surface    ');
fprintf('   phi_B =  %2.3e   a.u.\n',phi_total);
disp('   ');
disp('Divergence and curl of B at the grid point (xG,yG)')
fprintf('   xG =  %2.3f   mm\n',1e3*x(Nx));
fprintf('   yG =  %2.3f   mm\n',1e3*y(Ny));
disp('  ');
disp('Divergence of B at (xG,yG) - an element of the array divB'); 
fprintf('   div B =  %2.2e   T/m\n',divB(Ny,Nx));
disp('Divergence of B at (xG,yG) - calculated from definition of divergence'); 
fprintf('   div B =  %2.2e   T/m\n',divBp);
disp('  ');
disp('Curl of B at (xG,yG) - an element of the array curlB'); 
fprintf('   curl B =  %2.2e   T/m\n',curlB(Ny,Nx));
disp('  ');
disp('Ampere-Maxwell Law: curlB at location of wire 1');
fprintf(   'Wire 1: curlB  =  %2.2e   T/m\n',curlB_1);
disp('   ');
disp('To display in the Command Window the elements of the matrices');
disp('for the divergence or curl of the magnetic field B type')
disp('   divB   or  curlB   ');
disp('  ');


% % ======================================================================= 
% % GRAPHICS 
% % ======================================================================= 
   fs = 12;
   B(B > Bsat) = Bsat;
figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
   xP = 1e3.*xG; yP = 1e3.*yG; zP = 1e3 .* B;
   surf(xP,yP,zP);
   shading interp
   
   h = colorbar;
   h.Label.String = 'B   [ mT ]';
   colormap(parula);
   xlabel('x','fontsize',12); ylabel('y','fontsize',12);
   
   grid on
   box on
   
   set(gca,'fontsize',12)
     
figure(2)   %2222222222222222222222222222222222222222222222222222222222222
   set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
   xP = 1e3.*xG; yP = 1e3.*yG; zP = 1e3.*B;
  % pcolor(xP,yP,zP);
  % shading interp
  %  contourf(xP,yP,zP,0: 0.25: 5);
   contourf(xP,yP,zP,10);
   h = colorbar;
   h.Label.String = 'B   [ mT ]';
   colormap(parula);
   %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
   %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
   hold on
   % conductors
   for n = 1 : Nwire   
   col = [1 0 0];
   if I(n) < 0; col = [0 0 0]; end;
   p1 = xW(n)-radiusW; p2 = yW(n)-radiusW; p3 = 2*radiusW; p4 = 2*radiusW;
   pos = 1e3.*[p1, p2, p3, p4];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
   end
   
   xlabel('x  [mm]'); ylabel('y  [mm]');
   title('Magnetic field','fontweight','normal');
   
   axis square
   
figure(3)   % 33333333333333333333333333333333333333333333333333333333333
   set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
   
    hold on
     index1 = 1 : 20 : N;
    % index1 = [index1 500 502];
     index2 = index1;
          
     p1 = 1e3.*xG(index1, index2); p2 = 1e3.*yG(index1, index2);
    % Bscale = (B(index1,index2));
     Bscale = sqrt(Bx(index1,index2).^2 + By(index1,index2).^2);
     
     scale = 0.00001;
     Bscale(Bscale < scale)= scale; 
     % scaling of electric field lines: unit length
        p3 = 1e3.*Bx(index1, index2)./Bscale;
        p4 = 1e3.*By(index1, index2)./Bscale;
         % no scaling of electric field lines
       %p3 = Bx(index1, index2); p4 = By(index1, index2); 
     
     h = quiver(p1,p2,p3,p4,'autoscalefactor',0.6);
     set(h,'color',[0 0 1],'linewidth',1.2)
   
% conductors
   for n = 1 : Nwire   
   col = [1 0 0];
   if I(n) < 0; col = [0 0 0]; end;
   p11 = xW(n)-radiusW; p22 = yW(n)-radiusW; p33 = 2*radiusW; p44 = 2*radiusW;
   pos = 1e3.*[p11, p22, p33, p44];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
   end
     
     set(gca,'xLim', 1e3.*[minX maxX])
     set(gca,'yLim', 1e3.*[minY maxY])
     box on
     grid on
     hold on
     
     xlabel('x  [mm]'); ylabel('y  [mm]');
     title('Magnetic field vectors (unit length)','fontweight','normal');
     
     axis square
     
figure(4)   % 44444444444444444444444444444444444444444444444444444444
   set(gcf,'units','normalized','position',[0.73 0.52 0.23 0.32]);
   xP = 1e3.*xG; yP = 1e3.*yG; zP = 1e3.*B;
  % pcolor(xP,yP,zP);
  % shading interp
    contourf(xP,yP,zP,0: 0.25: 5);
   %contourf(xP,yP,zP,10);
   h = colorbar;
   h.Label.String = 'B   [ mT ]';
   colormap(parula);
   %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
   %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
   hold on
   % conductors
   for n = 1 : Nwire   
   col = [1 0 0];
   if I(n) < 0; col = [0 0 0]; end;
   p1 = xW(n)-radiusW; p2 = yW(n)-radiusW; p3 = 2*radiusW; p4 = 2*radiusW;
   pos = 1e3.*[p1, p2, p3, p4];
   h = rectangle('Position',pos,'Curvature',[1,1]);
   set(h,'FaceColor',col,'EdgeColor',col);
   end
   
   col = [1 1 1];
   pos = 1e3 .* [minLx minLy maxLx-minLx maxLy-minLy];
   h = rectangle('Position',pos);
   set(h,'EdgeColor',col);
  
   t1 = '\Sigma I = ';
   t2 = num2str(sum(I),'%2.3f A'); t3 = '    ';
   t4 = num2str(Ienclosed, 'I_{enclosed} =  %2.3f A');
   tm = [t1 t2 t3 t4];
   xlabel('x  [mm]'); ylabel('y  [mm]');
   title(tm,'fontweight','normal','fontsize',12);
   
   axis square
     
figure(5)   % 55555555555555555555555555555555555555555555555555555555555
   set(gcf,'units','normalized','position',[0.01 0.1 0.23 0.32]);        
   xP = r.*1e3; yP = B_A.*1e3;
   plot(xP,yP,'b','linewidth',2);
   hold on
   xP = -r.*1e3; yP = B_A.*1e3;
   plot(xP,yP,'b','linewidth',2);
   xP = x.*1e3; yP = B(ceil(N/2),:).*1e3;
   plot(xP,yP,'r','linewidth',2);
   
   xlabel('x  [mm]'); ylabel('B  [mT]');
   title('B-field blue(A: single wire)  red(N: multiple wires)','fontweight','normal');
   
   grid on

figure(6)   % 66666666666666666666666666666666666666666666666666666666666
   set(gcf,'units','normalized','position',[0.25 0.1 0.23 0.32]);
   xP = 1e3.*xG; yP = 1e3.*yG; zP = divB;
   surf(xP,yP,zP);
   shading interp
   
   h = colorbar;
   h.Label.String = 'divB   [ T/m ]';
   colormap(parula);
   xlabel('x','fontsize',12); ylabel('y','fontsize',12);
   title('Divergence of B');
   grid on
   box on
   
   set(gca,'fontsize',12)

   figure(7)   % 77777777777777777777777777777777777777777777777777777777
   set(gcf,'units','normalized','position',[0.49 0.1 0.28 0.32]);
   xP = 1e3.*xG; yP = 1e3.*yG; zP = curlB;
   surf(xP,yP,zP);
   shading interp
   
   h = colorbar;
   h.Label.String = 'curlB   [ T/m ]';
   colormap(parula);
   xlabel('x','fontsize',12); ylabel('y','fontsize',12);
   title('Curl of B');
   grid on
   box on
   
   set(gca,'fontsize',12)

   % =====================================================================
toc
