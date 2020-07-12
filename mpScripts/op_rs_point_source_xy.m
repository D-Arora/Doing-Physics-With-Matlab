close all;
clear all;
clc;

% op_rs_point_source_xy.m
% Circular apertures - irradiance in XY observation plane
% POINT SOURCE illumination of aperture
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]
% Calculation of energy enclosed in circles
% Uses functions
%     simpson1d.m  fn_distancePQ.m 
% 20 nov 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
    n1 = 60;              % Aperture: grid points for inner ring
    n2 = 120;             % Aperture: grid points for outer ring
    nR = 101;             % Aperture: number of rings    must be ODD
    nP = 121;             % Screen (Observation plane XY): must be odd
    wL = 632.8e-9;        % wavelength [m]
    a = 10*wL;            % radius of circular aperture  [m]

% Source
    xS = 0*wL; yS = 0*wL;
    zS = -50*wL;
    %zS = -1;
    ES = 1;

% Observation Space  [m]    
    yPmin = -55*wL;
    yPmax = 55*wL; 
    xPmin = -55*wL;
    xPmax = 55*wL; 
    zP = 100 * wL;               

% % Default values
%     n1 = 60;              % Aperture: grid points for inner ring
%     n2 = 120;             % Aperture: grid points for outer ring
%     nR = 101;             % Aperture: number of rings    must be ODD
%     nP = 121;             % Screen (Observation plane XY): must be odd
%     wL = 632.8e-9;        % wavelength [m]
%     a = 10*wL;            % radius of circular aperture  [m]
% 
% % Source
%     xS = 0*wL; yS = 0*wL;
%     zS = -50*wL;
%     %zS = -1;
%     ES = 1;
% 
% % Observation Space  [m]    
%     yPmin = -55*wL;
%     yPmax = 55*wL; 
%     xPmin = -55*wL;
%     xPmax = 55*wL; 
%     zP = 100 * wL;               


% ========================================================================
% SETUP 
% ========================================================================
    cL = 2.99792458e8;      % speed of light
    eps0 = 8.854187e-12;    % permittivity of free space
    nRI = 1;                % refractive index

    k = 2*pi/wL;            % propagation constant  [rad/s]
    ik = 1i*k;              % j k

    d_RL = 4*a^2/wL;        % Rayleigh distance  

% Aperture Space -------------------------------------------------------
    zQ = 0;
    A = zeros(nR,1);     % intgeral for each ring in aperture space
    n = zeros(nR,1);     % number of points for ring in aperture space

% Ring structure
%   radius of ring r  [m]     no. of data points around a ring n
%   Greater the circumference of a ring -->  more grid points
%   Width of each ring  dr
%   Total no. grid points for Aperture  nQ
rMax = a;
rMin = eps;
r = linspace(rMin, rMax, nR);
dr = r(2)-r(1);
m = (n2-n1) / (nR-1);
b = n2 - m * nR;

for c = 1 : nR
   n(c) = 2*round(0.5*(m * c + b))+1;
end
nQ = sum(n);

  
% Observation Space -----------------------------------------------------
    yP = linspace(yPmin,yPmax,nP);
    dyP = yP(2)-yP(1);
    xP = linspace(xPmin,xPmax,nP);
    dxP = xP(2)-xP(1);
    [xx,yy] = meshgrid(xP./wL,yP./wL);
    EP = zeros(nP,nP);

% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral

for cx = 1 : nP
for cy = 1 : nP    
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t);
      rSQ = fn_distancePQ(xS,yS,zS,xQ,yQ,zQ);
      EQ = ES * exp(ik .* rSQ) ./ rSQ;
      
      rPQ = fn_distancePQ(xP(cx),yP(cy),zP,xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      unit = ones(1,n(cQ));
      MP2 = zP .* (ik .* rPQ - unit);
      
      f = EQ .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cx,cy) = dr * sum(A)/ (2*pi);
end
end

uP = (cL*nRI*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));


% =======================================================================
% OUTPUT to Command Window
% =======================================================================
disp('Parameter summary  [SI units]');
fprintf('wavelength [m]  =  %3.5g \n',wL);
fprintf('no. of rings  nRings =  %3.3d \n',nR);
fprintf('inner ring n1  =  %3.3d \n',n1);
fprintf('outer ring n2  =  %3.3d \n',n2);
fprintf('nQ  =  %3.3d \n',nQ);
fprintf('nP  =  %3.3d \n',nP);

disp('  ')
fprintf('source  xS/wL  =  %3.3d \n',xS/wL);
fprintf('source  yS/wL  =  %3.3d \n',yS/wL);
fprintf('source  zS/wL  =  %3.3d \n',zS/wL);
disp('  ')
disp('Aperture Space');
fprintf('radius of aperture [m]  =  %3.3e \n',a);

disp('  ')
disp('Observation Space');
fprintf('observation space zP/wL  =  %3.3d \n',zP/wL);
disp('  ');


% =======================================================================
% GRAPHICS
% =======================================================================
 figure(1)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 14;
    pcolor(xx,yy,(uP./max(max(uP))).^0.3);
    shading interp
    tx = 'x_P / \lambda';
    ty = 'y_P / \lambda';
    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
    set(gca,'fontsize',fs);
    axis square
    colorbar 
    set(gca,'Xtick',-40:20:40); set(gca,'Ytick',-40:20:40);
    
        
figure(2)
     plot(xP./wL, uP(1,:));

figure(3)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 14;
    surf(xx,yy,(uP./max(max(uP))).^0.3);
    shading interp
    tx = 'x_P / \lambda';
    ty = 'y_P / \lambda';
    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
    set(gca,'fontsize',fs);
    axis square
    set(gca,'Xtick',-90:20:10); set(gca,'Ytick',-90:20:10);    
    axis off 
toc



