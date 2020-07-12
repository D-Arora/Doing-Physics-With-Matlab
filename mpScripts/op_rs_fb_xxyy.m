close all;
clear all;
clc;

% op_rs_fb_xxyy.m
% Circular aperture - % focused beam
% [1D] calculation - irradiance in XY observation plane
% POINT SOURCE illumination of aperture - convergent spherical waves
% Irradiance pattern can be asymmetrical about the optical axis (Z) 
% Numerical integration of the Rayleigh-Sommerfeld diffraction integral of
%    the first kind - Simpson's 1/3 rule
% Integration performed by dividing the aperture into rings 
% S.I. units used unless otherwise stated
% SYMBOLS:  irradiance = intensity = energy density w_e [W.m^-2]
%           radiant flux (flux)  aperture --> observation screen  W_E [W or J/s]
% Calculation of radiant flux enclosed in circles
% Uses functions
%     simpson1d.m  fn_distancePQ.m   turningPoints.m
% 11 dec oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic


% =======================================================================
% INPUTS 
% =======================================================================

% Aperture space
    wL = 632.8e-9;        % wavelength [m]
    a = 0.01;             % radius of circular aperture  [m]
    nQR = 15;               % number of rings - must be odd
    nQA = 13 * 4 + 1;               % additional partitons for each ring
    zQ = 0;
    
% Source  S
    zS = 0.1;
    xS = 0.2;
    yS = 0;
    ES = 1e-3;

% Observation space  P
    %nP = 14 * 4 + 1;             % Screen (Observation plane XY): must be odd  N * 4 + 1
    nPR = 55;
    nPA = 13 * 4 + 1;
    rPmax = 1.98e-5;         
    %xPmax = 200e-5;
    %zP = zS - 5e-5;   
    %zP = 10;
    zP = 0.1;
    % Default values --------------------------------------------------------


% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

NF = a^2 / (wL * zS);       % FRESNEL NUMBER
NA = a / sqrt(a^2 + zS^2);  % NUMERICAL APERTURE

% ========================================================================
% Aperture Space / Ring Structure 
% ========================================================================
    %   radius of ring r  [m]     no. of data points around a ring n
    %   Greater the circumference of a ring -->  more grid points
    %   Width of each ring  dr / grid points for a ring nQ
    %   Total no. grid points for Aperture  nQ_sum

    rQ = linspace(eps,a,nQR);  % radii of each ring
    drQ = rQ(2) - rQ(1);
    tQ = linspace(0,2*pi,nQA);  % angles for grid points
    xQ = zeros(nQR,nQA);  % X values for grid points
    yQ = zeros(nQR,nQA);  % Y values for grid points
    rSQ = zeros(nQR,nQA); % distance Q to S
    EQ = zeros(nQR,nQA);  % electric field at Q
    WE_Qring = zeros(1,nQR);  % flux for a ring


for c = 1 : nQR
    xQ(c,:) = rQ(c) .* cos(tQ);
    yQ(c,:) = rQ(c) .* sin(tQ);
   % tQ = tQ + 2*pi*c/nQR;    
end




% distance between aperture grid points and source point
% electric field at each aperture grid point [V/m]
for c = 1 : nQR
    rSQ(c,:) = fn_distancePQ(xS,yS,zS,xQ(c,:),yQ(c,:),zQ);
    EQ(c,:) = ES * exp(-ik .* rSQ(c,:)) ./ rSQ(c,:);
end
% irradiance distribution
    we_Q = (cL*nRI*eps0/2) .* abs(EQ).^2;
% radiant flux
    for c = 1 : nQR
        WE_Qring(c) = rQ(c) * drQ * simpson1d(we_Q(c,:),0,2*pi);
    end
        WE_Q = sum(WE_Qring);

        
% ========================================================================
% Observation / Ring Structure 
% ========================================================================    
    
    rP = linspace(eps,rPmax,nPR);  % radii of each ring
    drP = rP(2) - rP(1);

    tP = linspace(0,2*pi,nPA);  % angles for grid points
    xP = zeros(nPR,nPA);  % X values for grid points
    yP = zeros(nPR,nPA);  % Y values for grid points
       
for c = 1 : nPR
    xP(c,:) = rP(c) .* cos(tP)+xS;
    yP(c,:) = rP(c) .* sin(tP);
end

% X & Y axes / Optical coordinates - radial
%  nXY:   +X (1)   -X (3)   +Y (2)   -Y (4)

    nXY = [1 (nPA+1)/2 ; 1+(nPA-1)/4 nPA - (nPA-1)/4];
    vPxp = (2*pi*a/wL) .* xP(:,nXY(1)) ./ sqrt(xP(:,nXY(1)).^2 + zP^2);  
    vPxm = (2*pi*a/wL) .* xP(:,nXY(3)) ./ sqrt(xP(:,nXY(3)).^2 + zP^2);  
    vPyp = (2*pi*a/wL) .* yP(:,nXY(2)) ./ sqrt(yP(:,nXY(2)).^2 + zP^2);  
    vPym = (2*pi*a/wL) .* yP(:,nXY(4)) ./ sqrt(yP(:,nXY(4)).^2 + zP^2);  
%   
A = zeros(1,nQR);
EP = zeros(nPR,nPA);
for c1 = 1 : nPR            % P rings
for c2 = 1 : nPA            % P angles
for c3 = 1 : nQR            % Q rings    
    rPQ = fn_distancePQ(xP(c1,c2),yP(c1,c2),zP,xQ(c3,:),yQ(c3,:),zQ);
    rPQ3 = rPQ .* rPQ .* rPQ;
    kk = ik .* rPQ;
    MP1 = exp(kk);
    MP1 = MP1 ./ rPQ3;
    %unit = ones(1,nQ(c3));
    MP2 = zP .* (ik .* rPQ - 1);
    f = EQ(c3,:) .* MP1 .* MP2;
    A(c3) = rQ(c3) * simpson1d(f,0,2*pi); 

end
    EP(c1,c2) = drQ * sum(A)/ (2*pi);
end
end
   %EP(cP1,cP2) = drQ * sum(sum(A))/ (2*pi);
% end
% end
we_P = (cL*nRI*eps0/2) .* abs(EP).^2;
we_P_max = max(max(we_P));
we_P_dB = 10 .* log10(we_P./we_P_max);
 
% Energy received radiant flux - XY Observation plane ---------------------
    Psum = 0;
    WE_Psum = zeros(1,nPR);  WE_Pring = zeros(1,nPR);
    for c = 1 : nPR
        WE_Pring(c) = rP(c) * drP * simpson1d(we_P(c,:),0,2*pi);
        Psum = Psum + WE_Pring(c);
        WE_Psum(c) = Psum;
    end
        WE_P = sum(WE_Pring);
   
       nWE = 1; c = 1;     
   while  vPxp(c) < 1.22*pi
      c = c+1; 
      nWE = nWE + 1; 
   end

% % ======================================================================
% % % Turning Point Indices
% %    xData = xP; yData = we;
% %    [indexMin, indexMax] = turningPoints(xData, yData) ;
% % % Zeros
% %    xP_zeros = xP(indexMin);
% % % Maxima
% %    xP_maxs = xP(indexMax);
% %    we_maxs = we(indexMax);
% % 
% %     xData = xP; yData = IRRdB;
% %    [indexMinA, indexMaxA] = turningPoints(xData, yData) ;
% % % Zeros
% %    xP_zerosA = xP(indexMinA);
% % % Maxima
% %    xP_maxsA = xP(indexMaxA);
% %    we_maxsA = we(indexMaxA);
% 
% % =======================================================================
% % OUTPUT to Command Window
% % =======================================================================
disp('Parameter summary  [SI units]');
disp('APERTURE SPACE');
   fprintf('   Wavelength  wL         =  %3.5g  m  \n',wL);
   fprintf('   Radius of aperture a   =  %3.3f m \n',a);
   fprintf('   Fresnel Number  N_F    =  %3.3f   \n',NF);
   fprintf('   Numerical Aperture  N.A.  =  %3.3f   \n',NA);
   fprintf('   Number of rings  nQR   =  %3.3d \n',nQR);
   fprintf('   Number of angles nQA   =  %3.3d \n',nQA);
   fprintf('   Number grid points nQ  =  %3.3d \n',nQR * nQA);
   fprintf('   Max irradiance  we_Q   =  %3.3e  W/m^2 \n',max(max(we_Q)));
   fprintf('   Radiant flux WE_Q      =  %3.3e  W \n',WE_Q);

disp('  ')
   disp('SOURCE  ')
   fprintf('   xS  =  %3.3d  m \n',xS);
   fprintf('   yS  =  %3.3d  m \n',yS);
   fprintf('   zS  =  %3.3d  m \n',zS);
   fprintf('   Focal length  f  =  %3.3d  m \n',zS);
   fprintf('   Source strength  ES =  %3.3e  V/m \n',ES);

disp('  ')
   disp('OBSERVATION SPACE');
   fprintf('   Max radius rP  =  %3.3e  m \n',rPmax);
   fprintf('   Distance aperture to observation plane  zP = %3.3e  m \n',zP);
   fprintf('   Max irradiance  we_P    = %3.3e  W/m^2 \n',we_P_max);
   fprintf('   Radiant flux  WE_P = %3.3e  W \n',WE_P);
   fprintf('   Radiant flux enclosed within (v_P = 1.22 pi) = %3.2f percent \n',100 * WE_Psum(nWE) / WE_Q);
% 
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
% 
% 
% % =======================================================================
% % GRAPHICS
% % =======================================================================
% figure(99)
%  set(gcf,'Units','normalized');
%  set(gcf,'Position',[0.1 0.1 0.6 0.6]);
%  fs = 14;
%   x = xP(:,1); y = we_P(:,1);
%  plot(x,y);

% figure(1)
%  set(gcf,'Units','normalized');
%  set(gcf,'Position',[0.1 0.1 0.6 0.6]);
%  fs = 14;
%   x = xP; y = yP; z = we_P_dB;
%  surf(x,y,z);
%  shading interp
%  axis square
%  
%  figure(2)
%  set(gcf,'Units','normalized');
%  set(gcf,'Position',[0.1 0.1 0.6 0.6]);
%  fs = 14;
% 
%  contourf(xxP,yyP,we);
%  axis equal
%  
figure(1) % --------------------------------------------------------------
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.8 0.8]);
   fs = 15;

subplot(2,2,1)
   tx = 'x_P (blue) &  y_P (red)  [m]';
   ty = 'irradiance w_e  [W.m^{-2}]';
   x = xP(:,nXY(1))-xS;   % +X axis                
   y = we_P(:, nXY(1));
   plot(x,y,'b','linewidth',2)
   hold on
   x = xP(:,nXY(3))-xS;   % -X axis
   y = we_P(:, nXY(3));
   plot(x,y,'b','linewidth',2)

   x = yP(:,nXY(2));   % +Y axis
   y = we_P(:, nXY(2));
   plot(x,y,'r','linewidth',0.5)

   x = yP(:,nXY(4));   % -Y axis
   y = we_P(:, nXY(4));
   plot(x,y,'r','linewidth',0.5)

  xlabel(tx);
  ylabel(ty);
  grid on
  set(gca,'fontsize',fs);
  
subplot(2,2,2)
   tx = 'v_{Px} (blue) &  v_{Py} (red)  [m]';
   ty = 'irradiance w_e  [W.m^{-2}]';
   x = vPxp;   % +X axis                
   y = we_P(:, nXY(1));
   plot(x,y,'b','linewidth',2)
   hold on
   x = vPxm;   % -X axis
   y = we_P(:, nXY(3));
   plot(x,y,'b','linewidth',2)

   x = vPyp;   % +Y axis
   y = we_P(:, nXY(2));
   plot(x,y,'r','linewidth',0.5)

   x = vPym;   % -Y axis
   y = we_P(:, nXY(4));
   plot(x,y,'r','linewidth',0.5)

  xlabel(tx);
  ylabel(ty);
  grid on
  set(gca,'fontsize',fs);
  
  subplot(2,2,3)
   tx = 'x_P (blue) &  y_P (red)  [m]';
   ty = 'irradiance w_e  [dB]';
   x = xP(:,nXY(1));   % +X axis                
   y = we_P_dB(:, nXY(1));
   plot(x,y,'b','linewidth',2)
   hold on
   x = xP(:,nXY(3));   % -X axis
   y = we_P_dB(:, nXY(3));
   plot(x,y,'b','linewidth',2)

   x = yP(:,nXY(2));   % +Y axis
   y = we_P_dB(:, nXY(2));
   plot(x,y,'r','linewidth',0.5)

   x = yP(:,nXY(4));   % -Y axis
   y = we_P_dB(:, nXY(4));
   plot(x,y,'r','linewidth',0.5)

  xlabel(tx);
  ylabel(ty);
  grid on
  set(gca,'fontsize',fs);
  
subplot(2,2,4)
   tx = 'v_{Px} (blue) &  v_{Py} (red)  [m]';
   ty = 'irradiance w_e  [dB]';
   x = vPxp;   % +X axis                
   y = we_P_dB(:, nXY(1));
   plot(x,y,'b','linewidth',2)
   hold on
   x = vPxm;   % -X axis
   y = we_P_dB(:, nXY(3));
   plot(x,y,'b','linewidth',2)

   x = vPyp;   % +Y axis
   y = we_P_dB(:, nXY(2));
   plot(x,y,'r','linewidth',0.5)

   x = vPym;   % -Y axis
   y = we_P_dB(:, nXY(4));
   plot(x,y,'r','linewidth',0.5)

  xlabel(tx);
  ylabel(ty);
  grid on
  set(gca,'fontsize',fs);
  
figure(2)   % aperture EQ ------------------------------------------------
   x = xQ; y = yQ; z = we_Q ./ max(max(we_Q));
   pcolor(x,y,z)
   shading interp
   axis off
   colormap(hot)
   colorbar
   axis equal
   
figure(3)   % EP --------------------------------------------------------
   x = xP; y = yP; z = we_P_dB;
   pcolor(x,y,z)
   shading interp
   colormap(hot)
   colorbar 
   xlabel('xP  [m])');
   ylabel('yP  [m]')
   axis equal
   
figure(4)
   x = xP; y = yP; z = we_P_dB;
   surf(x,y,z)
   shading interp
   colormap(hot)
   %colorbar 
   axis off

figure(5)   % energy enclosed ------------------------------------------
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.2 0.2 0.3 0.3]);
   fs = 12;
   tx = 'optical coordinate  v_P / \pi ';
   ty = '% flux within circle of radius v_P ';
   x = vPxp./pi; y = 100 * WE_Psum./WE_Q;
   plot(x,y,'b','lineWidth',2);
   hold on
   x = [vPxp(nWE), vPxp(nWE)]/pi; y = [0,100];
   plot(x,y,'r');
   xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
   grid on
% % % Energy enclosed --------------------------------------------------------
% % figure(2)
% %    set(gcf,'Units','normalized');
% %     set(gcf,'Position',[0.2 0.2 0.3 0.3]);
% %    fs = 12;
% %    tx = 'optical coordinate  v_P / \pi ';
% %    ty = '% flux within circle of radius v_P ';
% %    x = vP./pi; y = 100 * WE_Psum./WE_Q;
% %    plot(x,y,'b','lineWidth',2);
% %    xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
% %    grid on
% %    set(gca,'fontsize',fs);
% %    hold on
% %    %x = [vP(min(indexMin))/pi,vP(min(indexMin))/pi];
% %    x = [1.22,1.22];
% %    y = [0, 100];
% %    plot(x,y,'r');
% %    set(gca,'Xlim',[0,6]);
% %    
% % % [3D] plots --------------------------------------------------------
% %   r = xP;
% %   t = linspace(0,2*pi,nP);
% % 
% %   xx = zeros(nP,nP);
% %   yy = zeros(nP,nP);
% % 
% %   for c = 1: nP
% %     xx(c,:) = r .* cos(t(c));
% %     yy(c,:) = r .* sin(t(c));
% %   end
% % 
% %   we_xy = meshgrid(we,we);
% % 
% %   
% % figure(3)
% % set(gcf,'Units','normalized');
% % set(gcf,'Position',[0.2 0.2 0.2 0.3]);
% % pcolor(xx,yy,10.*(we_xy).^0.3);
% % shading interp
% % axis equal
% % colormap(gray)
% % axis off
% % set(gcf,'color','k')
% % 
% % figure(4)
% % set(gcf,'Units','normalized');
% % set(gcf,'Position',[0.4 0.2 0.2 0.3]);
% % surf(xx,yy,10.*(we_xy).^0.3);
% % shading interp
% % colormap(jet)
% % axis off
% % set(gcf,'color','b')
% % 
% % 
% % figure(5)
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',[0.1 0.1 0.6 0.6]);
% %    fs = 14;
% %    tx = 'radial coordinate r_P  [m]';
% %    ty = 'irradiance  w_e  [ W.m^{-2} ]';
% %    x = xP; y = we;
% %    plot(x,y ,'linewidth',2);
% %    xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
% %    set(gca,'fontsize',fs);
% %    
% % figure(6)   
% %     set(gcf,'Units','normalized');
% %     set(gcf,'Position',[0.1 0.1 0.6 0.6]);
% %     fs = 14;
% %     tx = 'radial coordinate  v_P / \pi';
% %     ty = 'irradiance   w_e  [W.m^{-2}]';
% %     
% %     subplot(2,1,1)
% %     x = vP./pi; y = we ./ we_max;
% %     plot(x,y,'b','linewidth',2);
% %     hold on
% %      x = vP./pi; y = IRR;
% %     plot(x,y,'r','linewidth',1);
% %     xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
% %     set(gca,'Fontsize',fs);  
% %     set(gca,'Xlim',[0,6]);
% %     %set(gca,'Xtick',X3);
% %     grid on
% %     legend('N','A')
% %     
% %     subplot(2,1,2)
% %      ty = 'irradiance   w_e  [dB]';
% %     x = vP./pi;  y = we_dB;
% %     plot(x,y,'b','linewidth',2);
% %     hold on
% %     x = vP./pi; y = 10 .* log10(IRR);
% %     plot(x,y,'r','linewidth',1);
% %     xlabel(tx,'Fontsize',fs);   ylabel(ty,'Fontsize',fs);
% %     set(gca,'Fontsize',fs);  
% %     set(gca,'Xlim',[0,6]);
% %     %set(gca,'Ylim',[-80,0]);
% %     %set(gca,'Xtick',X3);
% %     grid on
% %     legend('N','A')
% %        
%    
% % =======================================================================
% 
   

% figure(3)
%    x = xQ; y = yQ;
%    plot(x,y,'o');
%    axis square
% 
% figure(4)
%    x = xP; y = yP;
%    plot(x,y,'o');
%    axis equal
%    
toc



