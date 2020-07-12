close all;
clear all;
clc;

% op_rs_fb_z.m
% Circular apertures - irradiance along optical axis   Z axis
% FOCUSSED BEAM - POINT SOURCE illumination of aperture
%    point source must be on optical  xS = 0  ys = 0
% Irradiance pattern has to be to be symmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
%     with increasing data points as radius increases 
% S.I. units used unless otherwise stated
% SYMBOLS:  we irradiance [W.m^-2]
% CALLS     simpson1d.m  fn_distancePQ.m   
% 21 nov oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
    n1 = 201;             % Aperture: grid points for inner ring
    n2 = 601;            % Aperture: grid points for outer ring
    nR = 221;             % Aperture: number of rings    must be ODD
    nP = 1809;             % Screen (Observation plane XY): must be odd
    wL = 1000*6e-9;        % wavelength [m]
    a = 300e-6;             % radius of circular aperture  [m]
    %a = 0.004472135955;
% Source S
    xS = 0;
    yS = 0;
   % zS = 0.5;
    zS = 0.15;
    %zS = 10;
    ES = 1e-3;
    
% Observation space P   range zP = (zS +/- delta_z)
    %delta_z = 0.003;
    %delta_z = 0.018;
    delta_z = 40e-3;
    zPmin = 0.001;
    zPmax = 0.20;
   % delta_z = 5;
    yP = 0;                        
    xP = 0;                        

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

NF = a^2 / (wL * zS);       % FRESNEL NUMBER
NA = a / sqrt(a^2 + zS^2);  % NUMERICAL APERTURE

% Observation Space -----------------------------------------------------
    %zPmin = zS - delta_z;
    %zPmax = zS + delta_z;            
    zP = linspace(zPmin,zPmax,nP);
    dzP = zP(2)-zP(1);
    %nc = 1 : nP;
    %zP = (zPmax - zPmin) .* (1-(cos(2*pi*nc / (4*nP))).^2) + zPmin;
    EP = zeros(1,nP);

    % optical coordinate
    uP = (k*a^2 / zS^2) .* (zP - zS);
    
    % theoretical irradiance IRR
    
    IRR = (sin(eps+uP/4)./(eps+uP/4)).^2;
% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    w_e energy density (irradiance) along the Z axis  [W/m^2]

     EQ = zeros(1,nR);
     for cQ = 1 : nR
        rQS = sqrt(zS^2 + r(cQ)^2);
        EQ(cQ) = ES * exp(-ik .* rQS) ./ rQS;
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
      
      fn = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(fn,0,2*pi);
   end
     EP(cP) = dr * sum(A)/ (2*pi);
end

we = (cL*nRI*eps0/2) .* abs(EP).^2;
we_max = max(we);
we_dB = 10 .* log10(we./we_max);


% ======================================================================
% % Turning Point Indices
%    % Numerical
%    xData = uP; yData = we;
%    [indexMin, indexMax] = turningPoints(xData, yData) ;
% % Zeros
%    uP_zeros = uP(indexMin);
% % Maxima
%    uP_maxs = uP(indexMax);
%    we_maxs = we(indexMax);
%     
%    % Analytical
%    xData = uP;
%    yData = IRR;
%    [indexMin2, indexMax2] = turningPoints(xData, yData) ;
% % Zeros
%    uP_IRR_zeros = uP(indexMin2);
% % Maxima
%    uP__IRR_maxs = uP(indexMax2);
%    we_IRR_maxs = we(indexMax2);


% =======================================================================
% OUTPUT to Command Window
% =======================================================================
disp('Parameter summary  [SI units]');
fprintf('wavelength  wL =  %3.5g  m \n',wL);
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
fprintf('   focal length  f  =  %3.3d  m \n',zS);
disp('   ');
disp('APERTURE SPACE');
fprintf('   radius of aperture a  =  %3.3f  m \n',a);
fprintf('   Fresnel Number  N_F  =  %3.3f   \n',NF);
fprintf('   Numerical Aperture  N.A.  =  %3.3f   \n',NA);
disp('  ')
disp('Observation Space - optical axis');
fprintf('max distance from focal plane (|zP-zS| =  %3.5g  m \n',delta_z);
%disp('  ');

% % TURNING POINTS
% disp('N: Radial coordinates - zero positions in energy density');
% for c = indexMin
%    fprintf('     %3.3f \n',uP(c)./pi);
% end
% disp('   ');
% disp('N Radial coordinates - max positions in energy density');
% disp('N: Relative intensities of peaks      ')
% for c = indexMax
%    fprintf('     %3.3f      %3.4f \n',uP(c)./pi, we(c)./we_max);
% end
% 
% disp('A: Radial coordinates - zero positions in energy density');
% for c = indexMin2
%    fprintf('     %3.3f \n',uP(c)./pi);
% end
% 
% disp('   ');
% disp('A: Radial coordinates - max positions in energy density');
% disp('Relative intensities of peaks      ')
% 
% for c = indexMax2
%    fprintf('     %3.3f      %3.4f \n',uP(c)./pi, IRR(c));
% end

% =======================================================================
% GRAPHICS
% =======================================================================

figure(1)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.6 0.6]);
   fs = 14;
   
   subplot(2,2,1);
   X1 = 0.198; X2 = 0.202;
   tx = 'z_P  [m] ';
   ty = 'irradiance  w_e  [ W.m^{-2} ]';
   x = zP; y = we;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Fontsize',fs);
   %set(gca,'Xlim',[X1, X2]);
   
   subplot(2,2,3); 
   ty = 'irradiance  w_e  [ dB ]';
   x = zP;   y = we_dB;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Fontsize',fs);   
   %set(gca,'Xlim',[X1, X2]);
   
   subplot(2,2,2);
   X3 = -20:4:20;
   tx = 'optical coordinate  u_P / \pi';
   ty = 'irradiance  w_e  [a.u.]';
   x = uP / pi;   y = we./we_max;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   grid on
   set(gca,'Xlim',[-20,20]);
   set(gca,'Fontsize',fs);
   set(gca,'Xtick',X3);
   
   subplot(2,2,4);   
   ty = 'irradaince  w_e  [dB]';
   x = uP /pi;   y = we_dB;
   plot(x,y ,'linewidth',2);
   xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
   set(gca,'Xlim',[-20,20]);
   grid on
   set(gca,'Fontsize',fs)
   set(gca,'Xtick',X3);
   
% figure(2)   
%     set(gcf,'Units','normalized');
%     set(gcf,'Position',[0.1 0.1 0.6 0.6]);
%     fs = 14;
%     tx = 'axial position  u_P / \pi';
%     ty = 'irradiance   w_e  [W.m^{-2}]';
%     
%     subplot(2,1,1)
%     x = uP./pi; y = we ./ we_max;
%     plot(x,y,'b','linewidth',2);
%     hold on
%      x = uP./pi; y = IRR;
%     plot(x,y,'r','linewidth',1);
%     xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%     set(gca,'Fontsize',fs);  
%     set(gca,'Xlim',[-20,20]);
%     set(gca,'Xtick',X3);
%     grid on
%     legend('N','A')
%     subplot(2,1,2)
%      ty = 'irradiance   w_e  [dB]';
%     x = uP./pi;  y = we_dB;
%     plot(x,y,'b','linewidth',2);
%     hold on
%     x = uP./pi; y = 10 .* log10(IRR);
%     plot(x,y,'r','linewidth',1);
%     xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%     set(gca,'Fontsize',fs);  
%     set(gca,'Xlim',[-20,20]);
%     set(gca,'Ylim',[-80,0]);
%     set(gca,'Xtick',X3);
%     grid on
%     legend('N','A')
%     
%   
%     
% figure(3)
% x = uP/pi;
% y = real(EP);
% plot(x,y)
% hold on
% y = imag(EP);
% plot(x,y);
% 
% figure(4)
% x = uP/pi;
% y = angle(EP);
% plot(x,y)
toc



