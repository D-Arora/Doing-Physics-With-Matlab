close all;
clear all;
clc;

% op_rs_point_source_z.m
% Circular apertures - irradiance along optical axis   Z axis
% POINT SOURCE illumination of aperture
%    point source must be on optical  xS = 0  ys = 0
% Irradiance pattern has to be to be symmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  we irradiance [W.m^-2]
% CALLS     simpson1d.m  fn_distancePQ.m   
% 19 nov oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
    n1 = 500;             % Aperture: grid points for inner ring
    n2 = 1200;            % Aperture: grid points for outer ring
    nR = 801;             % Aperture: number of rings    must be ODD
    nP = 809;             % Screen (Observation plane XY): must be odd
    wL = 632.8e-9;        % wavelength [m]
    a = 10*wL;             % radius of circular aperture  [m]

% Source
    xS = 0; yS = 0;
    zS = -1;
    %zS = -50*wL;
    ES = 1;

    zPmin = 0.1*wL;
    %zPmax = 120*wL;               % Aperture to Screen distance  [m]
    zPmax = 600*wL;
    yP = 0;                        % Observation point P
    xP = 0;                        % Max radial distance from Z axis

% % Default values
%  n1 = 500;             % Aperture: grid points for inner ring
%     n2 = 1200;            % Aperture: grid points for outer ring
%     nR = 801;             % Aperture: number of rings    must be ODD
%     nP = 809;             % Screen (Observation plane XY): must be odd
%     wL = 632.8e-9;        % wavelength [m]
%     a = 10*wL;             % radius of circular aperture  [m]
% 
% % Source
%     xS = 0; yS = 0;
%     zS = -1;
%     %zS = -50*wL;
%     ES = 1;
% 
%     zPmin = 0.1*wL;
%     %zPmax = 120*wL;               % Aperture to Screen distance  [m]
%     zPmax = 600*wL;
%     yP = 0;                        % Observation point P
%     xP = 0;                        % Max radial distance from Z axis


% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

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
%zP = linspace(zPmin,zPmax,nP);
%dzP = zP(2)-zP(1);
    nc = 1 : nP;
    zP = (zPmax - zPmin) .* (1-(cos(2*pi*nc / (4*nP))).^2) + zPmin;
    EP = zeros(1,nP);

% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    w_e energy density (irradiance) along the Z axis  [W/m^2]

     EQ = zeros(1,nR);
     for cQ = 1 : nR
        rSQ = sqrt(zS^2 + r(cQ)^2);
        EQ(cQ) = ES * exp(ik .* rSQ) ./ rSQ;
     end    

for cP = 1 : nP
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t); 
      rPQ = fn_distancePQ(xP,yP,zP(cP),xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      unit = ones(1,n(cQ));
      MP2 = zP(cP) .* (ik .* rPQ - unit);
      
      f = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cP) = dr * sum(A)/ (2*pi);
end

we = (cL*nRI*eps0/2) .* abs(EP).^2;
we_max = max(max(we));
we_dB = 10 .* log10(we./we_max);

% =======================================================================
% Theortical w_e
we_envelope_p = max(we) .* (1 + zP ./ sqrt(a^2 + zP.^2)).^2 ./ 4; 
we_envelope_m = max(we) .* (1 - zP ./ sqrt(a^2 + zP.^2)).^2 ./ 4;
we_theory =  max(we) .* (1 - 2.*zP ./ sqrt(a^2 + zP.^2) .* cos(k .* (sqrt(a^2 + zP.^2)-zP)) + zP.^2 ./(a^2 + zP.^2)) ./ 4;

% =======================================================================
% OUTPUT to Command Window
% =======================================================================
disp('Parameter summary  [SI units]');
fprintf('wavelength [m]  =  %3.5g \n',wL);
disp('PARTITIONS');
fprintf('   inner ring n1  =  %3.3d \n',n1);
fprintf('   outer inner ring n2  =  %3.3d \n',n2);
fprintf('   rings nR  =  %3.3d \n',nR);
fprintf('   nQ  =  %3.3d \n',nQ);
fprintf('   nP  =  %3.3d \n',nP);
disp('  ')
disp('SOURCE  ')
fprintf('   xS [m]    =  %3.3d \n',xS);
fprintf('   yS [m]    =  %3.3d \n',yS);
fprintf('   zS [m]    =  %3.3d \n',zS);
fprintf('   xS / wL   =  %3.3d \n',xS/wL);
fprintf('   yS / wL   =  %3.3d \n',yS/wL);
fprintf('   zS / wL   =  %3.3d \n',zS/wL);
disp('   ');
disp('APERTURE SPACE');
fprintf('   radius of aperture a / wL  =  %3.3f \n',a/wL);
disp('  ')
disp('Observation Space');
fprintf('   Rayleigh distance  d_RL / wL = %3.3f \n',d_RL/wL);
disp('  ');


% =======================================================================
% GRAPHICS
% =======================================================================

figure(1)  
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.6 0.6]);
   fs = 16;
   tx = 'axial position  z_P / \lambda';
   ty = 'irradiance   w_e  [W.m^{-2}]';
   x = zP ./ wL; y = we;
  
   plot(x(1:15:end),y(1:15:end) ,'bo','linewidth',1);
   hold on
   plot(x,y ,'b','linewidth',2);
   y = we_envelope_p;
   plot(x,y ,'r','linewidth',1);
   y = we_envelope_m;
   plot(x,y ,'r','linewidth',1);
   y = we_theory;
   plot(x,y ,'b','linewidth',1);
   set(gca,'Fontsize',fs);  
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   %set(gca,'Xlim',[0 120]);

figure(2)   
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.6 0.6]);
   fs = 16;
   tx = 'axial position  z_P / \lambda';
   ty = 'irradiance   w_e  [W.m^{-2}]';
   x = zP ./ wL; y = we;
  
   plot(x(1:15:end),y(1:15:end) ,'bo','linewidth',1);
   hold on
   plot(x,y ,'b','linewidth',2);
   y = we_envelope_p;
   plot(x,y ,'r','linewidth',1);
   y = we_envelope_m;
   plot(x,y ,'r','linewidth',1);
   y = we_theory;
   plot(x,y ,'b','linewidth',1);
   set(gca,'Fontsize',fs);  
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Xlim',[0 120]);
   
   
figure(3)   
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
    fs = 16;
    tx = 'axial position  z_P / \lambda';
    ty = 'irradiance   w_e  [W.m^{-2}]';
    x = zP ./ wL; y = we;
    plot(x,y,'b','linewidth',1);
   
   set(gca,'Fontsize',fs);  
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Xlim',[0 120]);
   

toc



