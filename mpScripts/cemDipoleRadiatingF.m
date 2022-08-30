% cemDipoleRadiatingE.m


close all
clc
clear

% Speed of light
 c = 2.99792458e8;

% Time t0
  t0 = 0;

% wavelength wL / propagation constant k /
% freq f / period / ang freq  w / w*t
  wL =  1;
  k = 2*pi/wL;
  f = c/wL;
  T = 1/f;
  w = 2*pi*f;
  wt = w*t0;

% Grid
% polar angle p P  (elevation)  cos(theta) = sin(p) sin(theta) = cos(p) 
% radial displacement r R
  N = 95;
  phi = pi/2;         % Azimuthal angle
  r1 = 0.2; r2 = 4;   % Radial displacement limits
  p = linspace(-pi,pi,N);    % elevation  p
  r  = linspace(r1,r2,N);
    [R, P] = meshgrid(r,p);
    Y = R.*sin(P);
    Z = R.*cos(P);
    
% Setup

  wt = w*t0;
  kR = k.*R;

% Electric field Polar component & Cartesian compoents
  Ep = (sin(P)./R) .* sin(kR - wt);
  Ex = cos(phi).*sin(P).*Ep;
  Ey = sin(phi).*sin(P).*Ep;
  Ez = sin(P).*Ep;
  E = sqrt(Ey.^2 + Ez.^2);

  Dy = 71; Dz = 30;
  yD = Y(Dy,1); zD = Z(Dz,1) ;
  ED = E(Dy,Dz);
  EDy = E(Dy,Dz);
  EDz = E(Dy,Dz);
  Dangle = atan2(EDz,EDz);

  zT = yD+1i*zD

  

% GRAPHICS  ========================================================
figure(1)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.30,0.30])

    contour(Y,Z,E.^0.2,12)    
    shading interp
    colormap(hot)
    colorbar

    hold on
    DrawArrow(zT,1,Dangle,0.2,0.2,2,[0 0 1])
    axis square
    xticks(-4:1:4); yticks(-4:1:4)
    xlim([-4 4])
    set(gca,'FontSize',12)

%figure(2)
%     set(gcf,'Units','normalized');
%     set(gcf,'Position',[0.1 0.35 0.30,0.30])
% 
%    [startY,startZ] = meshgrid(-2:2);
%     Eyn = Ey./sqrt(Ey.^2 + Ez.^2);
%     Ezn = Ez./sqrt(Ey.^2 + Ez.^2);
%     D = 1:10:N;
%     h2 = quiver(Y(D,D),Z(D,D),Eyn(D,D),Ezn(D,D));
%     set(h2,'AutoScale','on', 'AutoScaleFactor',0.5)
 %   stream2(Y,Z,Ey,Ez,startY,startZ)
  %  stream2(Y,Z,Ey,Ez,[1,1], [-1 1])
  %    starty= -2;% -2:1:2;
   %   startx = -2:1:2;
  %    startx=ones(size(starty))*min(Y,[],"all");
 %    hold on;
 %     stream2(Y,Y,Ey,Ey,1,1);
     
%      axis square
%      xlim([-2 2]); ylim([-2 2])
%    

  