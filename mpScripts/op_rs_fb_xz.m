close all;
clear all;
clc;

% op_rs_point_source_xz.m
% Circular apertures - irradiance in ZY MERIDONAL observation plane
% POINT SOURCE illumination of aperture
%    point source must be on optical  xS = 0  ys = 0
% a
% Irradiance pattern has to be to be symmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]
% Calculation of energy enclosed in circles
% Uses functions
%     simpson1d.m  fn_distancePQ.m   turningPoints.m

% 13 nov oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
n1 = 131;             % Aperture: grid points for inner ring
n2 = 221;            % Aperture: grid points for outer ring
nR = 121;             % Aperture: number of rings    must be ODD
nP = 181;             % Screen (Observation plane XY): must be odd
wL = 632.8e-9;        % wavelength [m]
a = 0.01;             % radius of circular aperture  [m]

% Source
    xS = 0; yS = 0;
    zS = 10;
    %zS = 0.1;
    ES = 1e-3;

% Observation Space  [m]    
delta_x = 2e-3;
delta_z = 3.8;
% delta_x = 0.015e-3;
% delta_z = 0.4e-3;
yP = 0;               


% Default values --------------------------------------------------------
% n1 = 100;             % Aperture: grid points for inner ring
% n2 = 1000;            % Aperture: grid points for outer ring
% nR = 1000;             % Aperture: number of rings
% nP = 509;             % Screen (Observation plane XY): must be odd
% wL = 632.8e-9;        % wavelength [m]
% a = 1e-4;             % radius of circular aperture  [m]
% uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]
% zP = 1;               % Aperture to Screen distance  [m] 
% yP = 0;               % Observation point P
% xPmax = 2e-2;         % Max radial distance from Z axis


% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

%EQmax = sqrt(2*uQmax/(cL*nRI*eps0));   % Aperture: Electric Field {V/m]

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

d_RL = 4*a^2/wL;          % Rayleigh distance  

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
    zPmin = zS - delta_z;   zPmax = zS + delta_z;
    zP = linspace(zPmin,zPmax,nP);
    dzP = zP(2)-zP(1);
    xPmin = xS - delta_x;   xPmax = xS + delta_x;
    xP = linspace(xPmin,xPmax,nP);
    dxP = xP(2)-xP(1);
    [zz,xx] = meshgrid(zP,xP);

     % optical coordinate
    uP = (k*a^2 / zS^2) .* (zz - zS);
    
    vP =  (k*a / zS) .* xx;
    
    EP = zeros(nP,nP);


% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    uP energy density (irradiance) along the X axis  [W/m^2]
% aperture electric field EQ / energy density / energy transmitted
     UQring = zeros(1,nR);  EQ = zeros(1,nR);
     for cQ = 1 : nR
        rSQ = sqrt(zS^2 + r(cQ)^2);
        EQ(cQ) = ES * exp(-ik .* rSQ) ./ rSQ;
     end


for cz = 1 : nP
for cx = 1 : nP    
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t); 
      rPQ = fn_distancePQ(xP(cx),yP,zP(cz),xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      unit = ones(1,n(cQ));
      MP2 = zP(cz) .* (ik .* rPQ - unit);
      
      fn = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(fn,0,2*pi);
   end
     EP(cx,cz) = dr * sum(A)/ (2*pi);
end
end

we = (cL*nRI*eps0/2) .* abs(EP).^2;
wemax = max(max(we));
wedB = 10.* log10(we./wemax);

NF = a^2 / (wL * zS);       % FRESNEL NUMBER
NA = a / sqrt(a^2 + zS^2);  % NUMERICAL APERTURE

% =======================================================================
% OUTPUT to Command Window
% =======================================================================
disp('Parameter summary  [SI units]');
fprintf('wavelength [m]  =  %3.5g \n',wL);
fprintf('nQ  =  %3.3d \n',nQ);
fprintf('nP  =  %3.3d \n',nP);
disp('  ')
disp('Aperture Space');
fprintf('radius of aperture [m]  =  %3.3e \n',a);
fprintf('   Fresnel Number  N_F  =  %3.3f   \n',NF);
fprintf('   Numerical Aperture  N.A.  =  %3.3f   \n',NA);

disp('  ')
disp('SOURCE  ')
fprintf('   xS  =  %3.3d  m \n',xS);
fprintf('   yS  =  %3.3d  m \n',yS);
fprintf('   zS  =  %3.3d  m \n',zS);
fprintf('   focal length  f  =  %3.3d  m \n',zS);
disp('   ');




% =======================================================================
% GRAPHICS
% =======================================================================
 figure(1)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 14;
    x = zz-zS; y = xx; z = wedB;
    pcolor(x,y,z);
    shading interp
    tx = 'axial position  z_P - z_S  [m]';
    ty = 'radial position r_P   [m]';
    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
    set(gca,'fontsize',fs);
    colormap(hot)
    colorbar
    
figure(3)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 14;
    x = uP./pi; y = vP./ pi; z = wedB;   %we.^0.1;
    pcolor(x,y,z);
    shading interp
    tx = 'u_P / \pi';
    ty = 'v_P / \pi';
    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
    set(gca,'fontsize',fs);
    set(gca,'Xtick',-12:4:12);
    colormap(hot)
    set(gca,'fontsize',fs);   
    
figure(4)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 14;
    x = uP./pi; y = vP./ pi; z = wedB;
    contourf(x,y,z,32);
    %shading interp
    tx = 'u_P / \pi';
    ty = 'v_P / \pi';
    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
    set(gca,'fontsize',fs);
    set(gca,'Xtick',-12:4:12);
    colormap(hot)
    hold on
    x = -4:1:4;
    y = x;
    plot(x,y,'k','linewidth',2);
    plot(x,-y,'k','linewidth',2);
    colorbar
    set(gca,'fontsize',fs);
    
toc



