close all;
clear all;
clc;

% op_rs_fb_xy1.m
% Circular aperture
% focused beam
% [1D] calculation - irradiance in XY observation plane
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
% 28 nov oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================
    n1 = 101;             % Aperture: grid points for inner ring
    n2 = 101;            % Aperture: grid points for outer ring
    nR = 185;             % Aperture: number of rings    must be ODD
    nP = 285;             % Screen (Observation plane XY): must be odd
    wL = 632.8e-9;        % wavelength [m]
    a = 0.01;             % radius of circular aperture  [m]

% Source  S
    zS = 0.1;
    %zS = 10;
    xS = 0; yS = 0;
    ES = 1e-3;

% Observation space  P
    %xPmax = 1.98e-5;         
    %xPmax = 200e-5;
    %zP = zS - 5e-5;   
    %zP = 10;
    zP = 0.1;
    yP = 0;               
    xP = 0.5e-5;
% Default values --------------------------------------------------------


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
rMin = eps;               % full aperture
%rMin = a/2;                % annulur aperture
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

% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    uP energy density (irradiance) along the X axis  [W/m^2]
%    UP energy from aperture to screen  [J/s]

% aperture electric field EQ / energy density / energy transmitted
UQring = zeros(1,nR);  EQ = zeros(1,nR);
     for cQ = 1 : nR
        rSQ = sqrt(zS^2 + r(cQ)^2);
        EQ(cQ) = ES * exp(-ik .* rSQ) ./ rSQ;
     end
  
       we_Qring = zeros(1,nR);
    for c = 1 : nR
       we_Q = (cL*nRI*eps0/2) .* abs(EQ(c)).^2;
       we_Qring(c) = 2 * pi * r(c) * dr * we_Q;
    end
       we_Qmax = max(we_Q);
       WE_Q = sum(we_Qring);

       
Asum = 0;       
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t); 
      rPQ = fn_distancePQ(xP,yP,zP,xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      unit = ones(1,n(cQ));
      MP2 = zP .* (ik .* rPQ - unit);
      
      fn = EQ(cQ) .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(fn,0,2*pi);
      
%      if (angle(A(cQ)/pi)) > 0.5;
%          A(cQ) = 0;
%      end
   end  
     EP = dr * sum(A)/ (2*pi);
     %EP(cP) = dr * Asum / (2*pi);
     %Asum = 0;
  
we = (cL*nRI*eps0/2) .* abs(EP).^2;
we_max = max(we);
we_dB = 10 .* log10(we./we_max);


% % ======================================================================
% % Turning Point Indices
%    xData = xP; yData = we;
%    [indexMin, indexMax] = turningPoints(xData, yData) ;
% % Zeros
%    xP_zeros = xP(indexMin);
% % Maxima
%    xP_maxs = xP(indexMax);
%    we_maxs = we(indexMax);
% 
%     xData = xP; yData = IRRdB;
%    [indexMinA, indexMaxA] = turningPoints(xData, yData) ;
% % Zeros
%    xP_zerosA = xP(indexMinA);
% % Maxima
%    xP_maxsA = xP(indexMaxA);
%    we_maxsA = we(indexMaxA);

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
fprintf('   focal length  f  =  %3.3d  m \n',zS);
fprintf('   source strength  E_S =  %3.3e  V/m \n',ES);


disp('   ');
disp('APERTURE SPACE');
fprintf('   radius of aperture a  =  %3.3f m \n',a);
fprintf('   max irradiance  we_Q  =  %3.3e  W/m^2 \n',we_Q);
fprintf('   radiant flux max irradiance  WE_Q  =  %3.3e  W \n',WE_Q);
fprintf('   Fresnel Number  N_F  =  %3.3f   \n',NF);
fprintf('   Numerical Aperture  N.A.  =  %3.3f   \n',NA);
disp('  ')

disp('OBSERVATION SPACE');
fprintf('   distance aperture to observation plane  zP = %3.3e  m \n',zP);
fprintf('   max irradiance  w_e    = %3.3e  W/m^2 \n',we_max);
%fprintf('   radiant flux  W_E = %3.3e  W \n',WE_P);

disp('   ');

% disp('N: radial coordinates  vP / pi  zero positions in irradiance');
% for c = indexMin
%    fprintf('     %3.3f \n',vP(c)/pi);
% end
% 
% disp('   ');
% disp('N: radial coordinates vP / pi - max positions in irradaince');
% disp('N: relative intensities of peaks      ')
% for c = indexMax
%    fprintf('     %3.3f      %3.4f \n',vP(c)/pi, we(c)./we_max);
% end
% 
% disp('A: radial coordinates  vP / pi  zero positions in irradaince');
% for c = indexMinA
%    fprintf('     %3.3f \n',vP(c)/pi);
% end
% 
% disp('   ');
% disp('A: radial coordinates vP / pi  max positions in irradaince');
% disp('A: relative intensities of peaks      ')
% for c = indexMaxA
%    fprintf('     %3.3f      %3.4f \n',vP(c)/pi, IRR(c));
% end


% =======================================================================
% GRAPHICS
% =======================================================================
figure(1)
 set(gcf,'Units','normalized');
 set(gcf,'Position',[0.1 0.1 0.6 0.6]);
 fs = 14;

 subplot(2,1,1)
 tx = 'ring number ';
 ty = 'Electric field ';
 x = 1:nR; y = real(A);
 plot(x,y ,'linewidth',2);
 hold on
 y = imag(A);
 plot(x,y ,'linewidth',2);
 y = conj(A) .* A;
 %plot(x,y ,'linewidth',2);

 subplot(2,1,2)
 tx = 'ring number ';
 ty = 'phase angle ';
 x = 1:nR; y = angle(A)/pi;
 plot(x,y ,'linewidth',2);
 

 
 %    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%    set(gca,'Fontsize',fs);
%    
% % energy density   uPdB  vs xP & yP
%    subplot(2,2,3);   
%    tx = 'radial coordinate  r_P  [m] ';
%    ty = 'irradiance  w_e  [ dB ]';
%    x = xP;   y = we_dB;
%    plot(x,y ,'linewidth',2);
%    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%    set(gca,'Fontsize',fs);   
% 
% % energy density   uP vs vPx & vPy
%    subplot(2,2,2);
%    tx = 'optical coordinate  v_P / \pi';
%    ty = 'irradiance  w_e  [a.u.]';
%    x = vP./pi;   y = we./we_max;
%    plot(x,y ,'linewidth',2);
%    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%    grid on
%    set(gca,'Xlim',[0,6]);
%    set(gca,'Fontsize',fs);
%    
% % energy density   uPdB  vs xP & yP
%    subplot(2,2,4);   
%    tx = 'optical coordinate   v_P / \pi';
%    ty = 'irradaince  w_e  [dB]';
%    x = vP./pi;   y = we_dB;
%    plot(x,y ,'linewidth',2);
%    xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%    set(gca,'Xlim',[0,6]);
%    grid on
%    set(gca,'Fontsize',fs);
% 
% % Energy enclosed --------------------------------------------------------
% figure(2)
%    set(gcf,'Units','normalized');
%     set(gcf,'Position',[0.2 0.2 0.3 0.3]);
%    fs = 12;
%    tx = 'optical coordinate  v_P / \pi ';
%    ty = '% flux within circle of radius v_P ';
%    x = vP./pi; y = 100 * WE_Psum./WE_Q;
%    plot(x,y,'b','lineWidth',2);
%    xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
%    grid on
%    set(gca,'fontsize',fs);
%    hold on
%    %x = [vP(min(indexMin))/pi,vP(min(indexMin))/pi];
%    x = [1.22,1.22];
%    y = [0, 100];
%    plot(x,y,'r');
%    set(gca,'Xlim',[0,6]);
%    
% % [3D] plots --------------------------------------------------------
%   r = xP;
%   t = linspace(0,2*pi,nP);
% 
%   xx = zeros(nP,nP);
%   yy = zeros(nP,nP);
% 
%   for c = 1: nP
%     xx(c,:) = r .* cos(t(c));
%     yy(c,:) = r .* sin(t(c));
%   end
% 
%   we_xy = meshgrid(we,we);
% 
%   
% figure(3)
% set(gcf,'Units','normalized');
% set(gcf,'Position',[0.2 0.2 0.2 0.3]);
% pcolor(xx,yy,10.*(we_xy).^0.3);
% shading interp
% axis equal
% colormap(gray)
% axis off
% set(gcf,'color','k')
% 
% figure(4)
% set(gcf,'Units','normalized');
% set(gcf,'Position',[0.4 0.2 0.2 0.3]);
% surf(xx,yy,10.*(we_xy).^0.3);
% shading interp
% colormap(jet)
% axis off
% set(gcf,'color','b')
% 
% 
% figure(5)
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
%    fs = 14;
%    tx = 'radial coordinate r_P  [m]';
%    ty = 'irradiance  w_e  [ W.m^{-2} ]';
%    x = xP; y = we;
%    plot(x,y ,'linewidth',2);
%    xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
%    set(gca,'fontsize',fs);
%    
% figure(6)   
%     set(gcf,'Units','normalized');
%     set(gcf,'Position',[0.1 0.1 0.6 0.6]);
%     fs = 14;
%     tx = 'radial coordinate  v_P / \pi';
%     ty = 'irradiance   w_e  [W.m^{-2}]';
%     
%     subplot(2,1,1)
%     x = vP./pi; y = we ./ we_max;
%     plot(x,y,'b','linewidth',2);
%     hold on
%      x = vP./pi; y = IRR;
%     plot(x,y,'r','linewidth',1);
%     xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%     set(gca,'Fontsize',fs);  
%     set(gca,'Xlim',[0,6]);
%     %set(gca,'Xtick',X3);
%     grid on
%     legend('N','A')
%     
%     subplot(2,1,2)
%      ty = 'irradiance   w_e  [dB]';
%     x = vP./pi;  y = we_dB;
%     plot(x,y,'b','linewidth',2);
%     hold on
%     x = vP./pi; y = 10 .* log10(IRR);
%     plot(x,y,'r','linewidth',1);
%     xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
%     set(gca,'Fontsize',fs);  
%     set(gca,'Xlim',[0,6]);
%     %set(gca,'Ylim',[-80,0]);
%     %set(gca,'Xtick',X3);
%     grid on
%     legend('N','A')
%        
%    
% % =======================================================================




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



