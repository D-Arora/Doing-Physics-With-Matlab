% cemEDR01.m

% ELECTRIC DIPOLE RADIATION
%   GENERATION OF ELECTROMAGNETIC WAVES
%   Dipole in Z direction

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Notes
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemEDR.htm
% SCRIPTS
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb


close all
clc
clear

% Speed of light
 c = 2.99792458e8;

% Time t0
  t0 = 1e-9;

% wavelength wL / propagation constant k /
% freq f / period / ang freq  w / w*t
  wL =  1;
  k = 2*pi/wL;
  f = c/wL;
  T = 1/f;
  w = 2*pi*f;
  wt = w*t0;

% Spatial Grid  YZ plane
% polar angle TH  measured from Z axis 
% radial displacement r R from Origin
  N = 195;            % Grid points
  phi = pi/2;         % Azimuthal angle
  r1 = 0.2; r2 = 4;   % Radial displacement limits
  th = linspace(-pi,pi,N);    % polar angle theta
  r  = linspace(r1,r2,N);
  [R, TH] = meshgrid(r,th);
  Y = R.*sin(TH);
  Z = R.*cos(TH);
  wt = w*t0;
  kR = k.*R;

% Electric field Polar component & Cartesian compoents
  Ep = -(sin(TH)./R) .* sin(kR - wt);
  Ex = cos(phi).*cos(TH).*Ep;
  Ey = sin(phi).*cos(TH).*Ep;
  Ez = -sin(TH).*Ep;
  E = sqrt(Ey.^2 + Ez.^2);

% Magnetic field
  theta = pi/2;
  phi = linspace(0,2*pi,N);
  [R, PHI] = meshgrid(r,phi);

  Bx = sin(PHI) .* sin(theta)./R.^2 .* (sin(kR - wt) - kR.*cos(kR - wt));
  By = cos(PHI) .* sin(theta)./R.^2 .* (sin(kR - wt) - kR.*cos(kR - wt));
  B = sin(theta)./R.^2 .* (sin(kR - wt) - kR.*cos(kR - wt));
  xx = R.* sin(PHI);
  yy = R.*cos(PHI);

% GRAPHICS  ========================================================
figure(1)   % YZ plot of electric field
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.25,0.35])
    hold on
    contour(Y,Z,E.^0.2,10,'linewidth',1.2)  
   % contourf(Y,Z,E.^0.2,12,'linewidth',1.2)
   
    shading interp
    colormap("parula")
    hbar = colorbar;
    set(hbar,'Ticks', [ ])
% Comment / Uncommnent to plot E-field vectors
        a = 1:12:N; b = 1:12:N;
        [aa,bb] = meshgrid(a,b);
    for ny = 1:length(a)
    for nz = 1: length(a)    
        yD = Y(a(ny),b(nz)); zD = Z(a(ny),b(nz)) ;
        if abs(yD) > 0.5 && abs(zD) > 0.5
        ED = E(a(ny),b(nz));
        EDy = Ey(a(ny),b(nz));
        EDz = Ez(a(ny),b(nz));
        Dangle = atan2(EDz,EDy);
        zT = yD+1i*zD; 
        DrawArrow(zT,0.3,Dangle,0.15,0.1,2,[0 0 1])
        end
    end
    end
       DrawArrow(-0.5i,1,pi/2,0.20,0.2,4,[1 0 0])

    axis square
    xticks(-4:1:4); yticks(-4:1:4)
    xlim([-4 4])
    box on
    xlabel('y  [ m ]');  ylabel('z  [ m ]')
    title('Electric Field E')
    set(gca,'FontSize',12)

% 22222222222222222222222222222222222222222222222222222222222222222222
figure(2)   % YZ plot of magnetic field
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.36 0.1 0.25,0.30])
    pcolor(xx,yy,B)
    %contourf(xx,yy,abs(B.^0.2),16);
    shading interp
    colormap("parula")
    hbar = colorbar;
    set(hbar,'Ticks', [ ])
    xticks(-4:1:4); yticks(-4:1:4)
    xlim([-4 4])
    box on
    xlabel('x  [ m ]');  ylabel('y  [ m ]')
    title('Magnetic Field B')
    set(gca,'FontSize',12)
    axis square

%%  3333333333333333333333333333333333333333333333333333333333333333333
figure(3)   % VECTORS E-field, B-field, Poynting vector S
            % YZ & XY planes
     set(gcf,'Units','normalized');
     set(gcf,'Position',[0.63 0.1 0.20,0.50])

     subplot(2,1,1)
     xP = [0 2]; yP = [0 2.5];
     plot(xP,yP,'k','linewidth',1)
     hold on
     zT = 2+2.5i; zTmag = 1.3; zTA = atan2(2.5,2);
     HL = 0.2; HW = 0.2; LW = 2; col = [1 0 0];
     DrawArrow(zT, zTmag, zTA, HL, HW, LW, col)  % S

     zT = 2+2.5i; zTmag = 1.3; zTA = pi/2+atan2(2.5,2);
     HL = 0.2; HW = 0.2; LW = 2; col = [0 0 1];
     DrawArrow(zT, zTmag, zTA, HL, HW, LW, col)  % B

     plot(2,2.5,'kx','markersize',20,'LineWidth',4);
    
     text(2.4,2.4,'E','FontWeight','bold','FontSize',14)
     text(0.4,3.4,'B','FontWeight','bold','FontSize',14,'Color','b')
     text(3.0,3.52,'S','FontWeight','bold','FontSize',14,'Color','r')
     xlim([-1 4]); ylim([-1 4])
     xticks(-1:1:4); yticks(-1:1:4)
     xlabel('y');  ylabel('z')
     set(gca,'FontSize',12)
     grid on, box on
     axis square

 subplot(2,1,2)
     A = linspace(0,2*pi,199);
     xP = cos(A); yP = sin(A);
     plot(xP,yP,'k','linewidth',1)
     hold on
     xP = [0 1]; yP = [0 0];
     plot(xP,yP,'k','linewidth',1)
     zT = 1+0i; zTmag = 1; zTA = 0; HL = 0.2; HW = 0.2; LW = 2; col = [1 0 0];
     DrawArrow(zT, zTmag, zTA, HL, HW, LW, col)  % S
     zT = 1+0i; zTmag = 1; zTA = pi/2; HL = 0.2; HW = 0.2; LW = 2; col = [0 0 1];
     DrawArrow(zT, zTmag, zTA, HL, HW, LW, col)  % S
     plot(1, 0,'kx','markersize',20,'LineWidth',4);
      plot(0,0,'mo','markersize',7,'LineWidth',4);
     text(1,-0.6,'E','FontWeight','bold','FontSize',14)
     text(0.9,1.3,'B','FontWeight','bold','FontSize',14,'Color','b')
     text(1.5,0.5,'S','FontWeight','bold','FontSize',14,'Color','r')
       xlim([-2 2]); ylim([-2 2])
       xticks(-2:1:2); yticks(-2:1:2)
     xlabel('x');  ylabel('y')
     set(gca,'FontSize',12)
     grid on, box on
     axis square
