% op_rs_xy_circle.m

% Circular aperture aperture
% Numerical integration of the Rayleigh-Sommerfield diffraction integral of
% the first kind
% Uses cartesain coordinates / S.I. units
% Calculate of energy enclosed in circles
% 
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]


% Uses functions
%     simpson2d.m  fn_distancePQ.m

% 25 oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************
clear all
close all
%clc

tic

% ========================================================================
% INPUT PARAMETERS 
num = 50;             % number for observation space
nP = num*4+1;         % observation points for P  format  integer * 4 + 1
nQ = 29;              % aperture points for Q  must be ODD

wL = 632.8e-9;        % wavelength [m]
 
a = 10*wL;             % radius of circular aperture  [m]

uQmax = 1e-3;         % incident energy density  [W.m^-2]

xPmax = 100*wL;         % half width - observation space   [m] 
yPmax = 100*wL;
zP = 600*wL;               % distance: aperture plan to observation plane [m]

% default values
% num = 80;             
% nP = num*4+1;         % observation points for P  format  integer * 4 + 1
% nQ = 99;              % aperture points for Q  must be ODD
% wL = 632.8e-9;        
% a = 1e-4; 
% xPmax = 2e-2;          
% yPmax = 2e-2;
% zP = 1;               

% ========================================================================
% SETUP 
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nR = 1;                 % refractive index

k = 2*pi/wL;            % propagation constant
ik = 1i*k;              % j k

d_RL = a^2/wL;          % Rayleigh distance  

% Initialize matrices ----------------------------------------------------
unit = ones(nQ,nQ); 
rPQ = zeros(nQ,nQ); rPQ3 = zeros(nQ,nQ);
MP1 = zeros(nQ,nQ); MP2 = zeros(nQ,nQ); kk = zeros(nQ,nQ); 
MP = zeros(nQ,nQ);
EP = zeros(nP,nP);  

% Aperture space --------------------------------------------------------- 
xQmin = -a;  xQmax = a;
yQmin = -a;  yQmax =  a;
zQ = 0;
xQ1 = linspace(xQmin,xQmax,nQ);    
yQ1 = linspace(yQmin,yQmax,nQ);
[xQ, yQ] = meshgrid(xQ1,yQ1);
rQ = sqrt(xQ.^2 + yQ.^2);

% Electric Eield EQ / Energy density uQ / Energy from aperture UQ  [J/s]
EQmax = sqrt(2*uQmax/(cL*nR*eps0)); 
EQ = EQmax .* ones(nQ,nQ); 
EQ(rQ > a) = 0;
  
uQ = (cL*nR*eps0/2) .* EQ .*EQ;       
UQ = simpson2d(uQ,xQmin,xQmax,yQmin,yQmax);  
UQtheory = uQmax * pi* a^2;

% Observation space -----------------------------------------------------
xPmin = -xPmax;
yPmin =  -yPmax;

xP1 = linspace(xPmin,xPmax,nP);     
yP1 = linspace(yPmin,yPmax,nP);
[xP, yP] = meshgrid(xP1,yP1);
indexXY = num*2+1;          % defines X and Y axes  y = 0 and x = 0;
rP = sqrt(xP.^2 + yP.^2);

% optical coordinates
vP = (2*pi*a/wL) .* xP ./ sqrt(xP.^2 + zP^2);  

% first min xP and yP values
theta = asin(wL*3.831/(2*pi*a));
x1 = zP * tan(theta);


% =======================================================================
% COMPUATION OF DIFFRACTION INTEGRAL 

% Electric Field EP - observation space   P
for c1 = 1 : nP
    for c2 = 1 :nP
       rPQ = fn_distancePQ(xP(c1,c2),yP(c1,c2),zP,xQ,yQ,zQ);
       rPQ3 = rPQ .* rPQ .* rPQ;
       kk = ik .* rPQ;
       MP1 = exp(kk);
       MP1 = MP1 ./ rPQ3;
       MP2 = zP .* (ik .* rPQ - unit);
             
       MP = EQ .* MP1 .* MP2;
       EP(c1,c2) = simpson2d(MP,xQmin,xQmax,yQmin,yQmax) ./(2*pi) ;
   end 
end

% Irradiance (energy density) u /  energy U  -   observation space   P 
uP = (cL*nR*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));
uPdB = 10 .* log10(uP./uPmax);

UP = simpson2d(uP,xPmin,xPmax,yPmin,yPmax);

% % percentage of energy enclosed within a circle of radius of 1st min in x direction
% Rc = 3.831;
% uPenclosed = uP;
% %uPenclosed(rP > 3.875e-3) = 0;
% uPenclosed(abs(vP) > Rc) = 0;
% UPenclosed = simpson2d(uPenclosed,xPmin,xPmax,yPmin,yPmax);
% UPenclosed_perc = 100*UPenclosed/UQ;

%theta_c = asin(Rc/(k*ax));     % ax effective diameter of circular aperture
%xPc = zP * tan(theta_c);

% % Draw circle
% nc = 500;
% theta = linspace(0,2*pi,nc);
% xc = x1 .* cos(theta); yc = x1 .* sin(theta);


% GRAPHICS --------------------------------------------------------------
% figure(1)
% set(gcf,'Units','normalized');
% set(gcf,'Position',[0.1 0.1 0.6 0.6]);
% fs = 12;
% %tx = 'x_P (blue)  or y_P (red)  [m]';
% %ty = 'irradiance  (a.u.)';
% 
% % energy density   uP vs Xp & yP
% subplot(2,2,1);
% tx = 'x_P (blue)  or y_P (red)  [m]';
% ty = 'energy density  u  [W.m^{-2}]';
% x = xP(indexXY,:); y = uP(indexXY,:);
% plot(x,y ,'linewidth',2);
% hold on
% x = yP(:,indexXY); y = uP(:,indexXY);
% plot(x,y,'linewidth',2);
% xlabel(tx);   ylabel(ty);
% 
% % energy density   uPdB  vs xP & yP
% subplot(2,2,3);   
% tx = 'x_P (blue)  or y_P (red)  [m]';
% ty = 'energy density  u  [dB]';
% x = xP(indexXY,:);
% y = uPdB(indexXY,:);
% plot(x,y ,'linewidth',2);
% % hold on
% % x = yP(:,indexXY);
% % y = uPdB(:,indexXY);
% % plot(x,y,'linewidth',2)
% xlabel(tx);   ylabel(ty);
% 
% % energy density   uP vs vPx & vPy
% subplot(2,2,2);
% tx = 'optical coordinate  v_P';
% ty = 'energy density  u  [a.u.]';
% x = vP(indexXY,:);
% y = uP(indexXY,:);
% plot(x,y ,'linewidth',2);
% xlabel(tx);   ylabel(ty);
% grid on
% 
% % energy density   uPdB  vs xP & yP
% subplot(2,2,4);   
% tx = 'optical coordinate   v_P';
% ty = 'energy density  u  [dB]';
% x = vP(indexXY,:);
% y = uPdB(indexXY,:)./uPmax;
% plot(x,y ,'linewidth',2);
% xlabel(tx);   ylabel(ty);
% grid on
% 
% figure(2)      % photograph plot  /  [3D] plot
% set(gcf,'Units','normalized');
% set(gcf,'Position',[0.2 0.2 0.3 0.3]);
% pcolor(xP,yP,10.*(uP).^0.3);
% shading interp
% axis equal
% colormap(gray)
% axis off
% set(gcf,'color','k')
% 
% hold on
% x = xc; y = yc;
% plot(x,y,'y','lineWidth',2);
% 
% figure(3)
% set(gcf,'Units','normalized');
% set(gcf,'Position',[0.3 0.3 0.3 0.3]);
% surf(xP,yP,10.*(uP).^0.3);
% shading interp
% axis equal
% colormap(jet)
% axis off
% %set(gcf,'color','k')
% xlabel('x');
% ylabel('y')
% 
% 
% % figure(4)
% % set(gcf,'Units','normalized');
% % set(gcf,'Position',[0.3 0.3 0.15 0.2]);
% % %tx = 'x_P (blue)  or y_P (red)  [m]';
% % tx = 'x_P  [m]';
% % %ty = 'energy density  u  [W.m^{-2}]';
% % ty = 'u  [W.m^{-2}]';
% % tm1 ='    z_P  =   ';
% % tm2 = num2str(zP,3);
% % tm3 = '  m';
% % tm = [tm1 tm2 tm3];
% % x = xP(indexXY,:); y = uP(indexXY,:);
% % plot(x,y ,'b','linewidth',2);
% % %hold on
% % %x = yP(:,indexXY); y = uP(:,indexXY);
% % %plot(x,y,'r','linewidth',2);
% % xlabel(tx);   ylabel(ty);
% % title(tm);
% % grid on
% % %set(gca,'Ylim',[0 1e-6]);
% % %set(gca,'Xlim',[-4e-4 4e-4]);
% 
% 
figure(6)
fs = 14;
x = xQ; y = yQ;
pcolor(x,y,real(EQ));
shading flat
tx = 'x_Q /  \lambda';
ty = 'y_Q /  \lambda';
set(gca,'fontsize',fs);
xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
axis square
% 
% % COMMAND WINDOW OUTPUTS ------------------------------------------------
% disp('Parameter summary  [SI units]');
% fprintf('wavelength [m]  =  %3.3g \n',wL);
% fprintf('nQ  =  %3.3d \n',nQ);
% fprintf('nP  =  %3.3d \n',nP);
% disp('  ')
% disp('Aperature Space');
% fprintf('radius of aperture [m]  =  %3.3e \n',a);
% fprintf('energy density [W/m2] uQmax  =  %3.3e \n',uQmax);
% fprintf('energy from aperture [J/s]   UQ(theory) = %3.3e \n',UQtheory);
% fprintf('energy from aperture [J/s]   UQ(calculated) = %3.3e \n',UQ);
% disp('  ')
% disp('Observation Space');
% fprintf('X width [m] =  %3.3e \n',2*xPmax);
% fprintf('Y width [m]  =  %3.3e \n',2*yPmax);
% fprintf('distance aperture to observation plane [m]   zP = %3.3e \n',zP);
% fprintf('Rayleigh distance  [m]   d_RL = %3.3e \n',d_RL);
% %fprintf('Fraunhofer: position of 1st min in X direction [m]  =  %3.3e \n',x1);
% %fprintf('Fraunhofer: position of 1st min in Y direction [m]  =  %3.3e \n',y1);
% disp('  ');
% fprintf('energy to aperture [J/s]   UP = %3.3e \n',UP);
% fprintf('max energy density  [W./m2]   uPmax = %3.3e \n',uPmax);
% disp('   ');
% fprintf('percentage energy enclosed within circle of radius (vP = 3.831)  = %3.1f \n',UPenclosed_perc);
% disp('   ');

figure(5)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 16;

% energy density   uP vs Xp & yP
   tx = 'radial coordinate r_P / \lambda';
   ty = 'energy density  u  [W.m^{-2}]';
  % x = xP ./ wL; y = uP;
  x = xP(indexXY,:)./wL; y = uP(indexXY,:);
   plot(x,y ,'linewidth',2);
   xlabel(tx);   ylabel(ty);

toc


