close all;
clear all;
clc;

% op_rs_source.m
% Circular apertures - irradiance in XY observation plane
% POINT SOURCE illumination of aperture
%    point source must be on optical  xS = 0  ys = 0
% Irradiance pattern has to be to be symmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  irradiance = intensity = energy density w_e [W.m^-2]
%           radiant flux (flux)  aperture --> observation screen  W_E [W or J/s]
% Calculation of radiant flux enclosed in circles
% Uses functions
%     simpson1d.m  fn_distancePQ.m   turningPoints.m
% 20 nov oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
n1 = 200;             % Aperture: grid points for inner ring
n2 = 800;            % Aperture: grid points for outer ring
nR = 801;             % Aperture: number of rings    must be ODD
nP = 689;             % Screen (Observation plane XY): must be odd
wL = 632.8e-9;        % wavelength [m]
a = 10*wL;             % radius of circular aperture  [m]

% Source
    xS = 0; yS = 0;
    %zS = -1;
    zS = - 50 * wL;
    ES = 1;


zP = 600*wL;               % Aperture to Screen distance  [m]
yP = 0;               % Observation point P
xPmax = 200*wL;         % Max radial distance from Z axis


% Default values --------------------------------------------------------


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
xP = linspace(0,xPmax,nP);
dxP = xP(2)-xP(1);
% optical coordinates
vP = (2*pi*a/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
EP = zeros(1,nP);


% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    uP energy density (irradiance) along the X axis  [W/m^2]
%    UP energy from aperture to screen  [J/s]

% aperture electric field EQ / energy density / energy transmitted
UQring = zeros(1,nR);  EQ = zeros(1,nR);
     for cQ = 1 : nR
        rSQ = sqrt(zS^2 + r(cQ)^2);
        EQ(cQ) = ES * exp(ik .* rSQ) ./ rSQ;
     end
  
       we_Qring = zeros(1,nR);
    for c = 1 : nR
       we_Q = (cL*nRI*eps0/2) .* abs(EQ(c)).^2;
       we_Qring(c) = 2 * pi * r(c) * dr * we_Q;
    end
       we_Qmax = max(we_Q);
       WE_Q = sum(we_Qring);
            
for cP = 1 : nP
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t); 
      rPQ = fn_distancePQ(xP(cP),yP,zP,xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      unit = ones(1,n(cQ));
      MP2 = zP .* (ik .* rPQ - unit);
      
      f = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cP) = dr * sum(A)/ (2*pi);
end

we = (cL*nRI*eps0/2) .* abs(EP).^2;
we_max = max(we);
we_dB = 10 .* log10(we./we_max);

% Energy received - XY Observation plane  -------------------------------
%   AP      energy within a ring  [J/s]
%   UPsum   energy within XY Observation plane  [J/s]
%   UPsum   energy enclosed within rings of increasing radius  
AP    = zeros(1,nP); 
WE_Psum = zeros(1,nP);
for cP = 1 : nP
    AP(cP) = 2 * pi * xP(cP) * we(cP) * dxP ;
end
WE_P = sum(AP);

for cP = 1 : nP-1
 WE_Psum(cP+1) = WE_Psum(cP) + AP(cP);
end
 WE_Pmax = max(WE_Psum);

% ======================================================================
% Turning Point Indices
   xData = xP; yData = we;
   [indexMin, indexMax] = turningPoints(xData, yData) ;
% Zeros
   xP_zeros = xP(indexMin);
% Maxima
   xP_maxs = xP(indexMax);
   we_maxs = we(indexMax);
   

% =======================================================================
% OUTPUT to Command Window
% =======================================================================
disp('Parameter summary  [SI units]');
fprintf('wavelength  wL  =  %3.5g  m  \n',wL);
disp('PARTITIONS');
fprintf('   inner ring n1  =  %3.3d \n',n1);
fprintf('   outer inner ring n2  =  %3.3d \n',n2);
fprintf('   rings nR  =  %3.3d \n',nR);
fprintf('   nQ  =  %3.3d \n',nQ);
fprintf('   nP  =  %3.3d \n',nP);

disp('  ')
disp('SOURCE  ')
fprintf('   xS  =  %3.3d  m \n',xS);
fprintf('   yS  =  %3.3d  m \n',yS);
fprintf('   zS  =  %3.3d  m \n',zS);
fprintf('   xS / wL   =  %3.3d \n',xS/wL);
fprintf('   yS / wL   =  %3.3d \n',yS/wL);
fprintf('   zS / wL   =  %3.3d \n',zS/wL);

disp('   ');
disp('APERTURE SPACE');
fprintf('   radius of aperture a / wL  =  %3.3f \n',a/wL);
fprintf('   max irradiance  we_Q  =  %3.3e  W/m^2 \n',we_Q);
fprintf('   radiant flux max irradiance  WE_Q  =  %3.3e  W \n',WE_Q);
disp('  ')

disp('OBSERVATION SPACE');
fprintf('   max radius rP  =  %3.3e  m \n',xPmax);
fprintf('   distance aperture to observation plane  zP = %3.3e  m \n',zP);
fprintf('   Rayleigh distance  d_RL = %3.3e  m \n',d_RL);
fprintf('   distance aperture to observation plane  zP / wL = %3.3e \n',zP/wL);
fprintf('   Rayleigh distance  d_RL / wL = %3.3e \n',d_RL/wL);
fprintf('   max irradiance  w_e    = %3.3e  W/m^2 \n',we_max);
fprintf('   radiant flux  W_E = %3.3e  W \n',WE_P);
disp('   ');
disp('   ');

disp('Radial coordinates - zero positions in energy density');
for c = indexMin
   fprintf('     %3.3f \n',vP(c));
end

disp('   ');
disp('Radial coordinates - max positions in energy density');
disp('Relative intensities of peaks      ')
for c = indexMax
   fprintf('     %3.3f      %3.4f \n',vP(c), we(c)./we_max);
end

% =======================================================================
% GRAPHICS
% =======================================================================
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 14;

% energy density   uP vs Xp & yP
   subplot(2,2,1);
   tx = 'radial coordinate r_P / \lambda ';
   ty = 'irradiance  w_e  [ W.m^{-2} ]';
   x = xP./wL; y = we;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Fontsize',fs);
   
% energy density   uPdB  vs xP & yP
   subplot(2,2,3);   
   tx = 'radial coordinate  r_P / \lambda ';
   ty = 'irradiance  w_e  [ dB ]';
   x = xP./wL;   y = we_dB;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Fontsize',fs);   

% energy density   uP vs vPx & vPy
   subplot(2,2,2);
   tx = 'optical coordinate  v_P';
   ty = 'irradiance  w_e  [a.u.]';
   x = vP;   y = we./we_max;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   grid on
   set(gca,'Xlim',[0,20]);
   set(gca,'Fontsize',fs);
   
% energy density   uPdB  vs xP & yP
   subplot(2,2,4);   
   tx = 'optical coordinate   v_P';
   ty = 'irradaince  w_e  [dB]';
   x = vP;   y = we_dB;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Xlim',[0,20]);
   grid on
   set(gca,'Fontsize',fs);

% Energy enclosed --------------------------------------------------------
figure(2)
   tx = 'optical coordinate  v_P ';
   ty = '% flux within circle of radius v_P ';
   x = vP; y = 100 * WE_Psum./WE_Q;
   plot(x,y,'lineWidth',2);
   xlabel(tx);   ylabel(ty);
   grid on

figure(22)
   fs = 12;
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.3 0.3]);
   tx = 'radial coordinate  r_P / \lambda';
   ty = '% flux within circle of radius r_P ';
   x = xP./wL; y = 100 * WE_Psum./WE_Q;
   plot(x,y,'lineWidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   hold on
   x = [a/wL, a/wL]; y = [0,100];
   plot(x,y,'r','lineWidth',1);
   grid on
   set(gca,'Fontsize',fs);

% [3D] plots --------------------------------------------------------
  r = xP;
  t = linspace(0,2*pi,nP);

  xx = zeros(nP,nP);
  yy = zeros(nP,nP);

  for c = 1: nP
    xx(c,:) = r .* cos(t(c));
    yy(c,:) = r .* sin(t(c));
  end

  we_xy = meshgrid(we,we);

  
figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.3]);
pcolor(xx,yy,10.*(we_xy).^0.3);
shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')

figure(4)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.4 0.2 0.2 0.3]);
surf(xx,yy,10.*(we_xy).^0.3);
shading interp
colormap(jet)
axis off
set(gcf,'color','b')


figure(5)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 16;

% irradiance uP vs Xp & yP
   tx = 'radial coordinate r_P / \lambda';
   ty = 'irradiance  w_e  [ W.m^{-2} ]';
   x = xP ./ wL; y = we;
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);

% =======================================================================




% figure(1)
% for c = 1 : nR
%    t = linspace(0,2*pi,n(c));  
%    x = r(c) .* cos(t);
%    y = r(c) .* sin(t); 
%    plot(x,y,'o');
%    axis square
%    hold on
% end

toc



