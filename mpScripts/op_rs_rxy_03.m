% op_rs_xy_03.m

% Rectangular aperture space: miscellaneous shaped aperture openings
% Numerical integration of the Rayleigh-Sommerfield diffraction integral of
% the first kind
% Uses cartesain coordinates / S.I. units

% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]


% Uses functions
%     simpson2d.m  fn_distancePQ.m

% 17 oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************
clear all
close all
clc

tic


% INPUT PARAMETERS -------------------------------------------------------
num = 100;             % number for observation space
nP = num*4+1;         % observation points for P  format  integer * 4 + 1
nQ = 59;              % aperture points for Q  must be ODD

wL = 650e-9;        % wavelength [m]
 
%ax = 4e-4;             % width - aperture space  [m]
ay = 1e-4;  

% double slit inputs
ax1 = 0.015e-3;
ax2 = 0.06e-3;


uQmax = 1e-3;                % incident energy density  [W.m^-2]

xPmax = 10e-2;         % half width - observation space   [m] 
yPmax = 10e-2;
%zP = 1e-3;                % distance: aperture plan to observation plane [m]
zP = 1;

% default values
% num = 30;             
% nP = num*4+1;         % observation points for P  format  integer * 4 + 1
% nQ = 99;              % aperture points for Q  must be ODD
% wL = 632.8e-9;        
% xQmax = 1e-4;        
% yQmax = 2e-4;  
% uQmax = 1e-3;               
% xPmax = 2e-2;          
% yPmax = 2e-2;
% zP = 1;               


% SETUP -----------------------------------------------------------------
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nR = 1;                 % refractive index

k = 2*pi/wL;            % propagation constant
ik = 1i*k;              % j k

% Aperture space

ax = ax1 + ax2;
bx = 0.5*(ax2-ax1);
ay = ax1;
xQmin = -ax/2;  xQmax = ax/2;
yQmin =  -ay/2; yQmax =  ay/2;

unit = ones(nQ,nQ);    % unit matrix
rPQ = zeros(nQ,nQ); rPQ3 = zeros(nQ,nQ);
MP1 = zeros(nQ,nQ); MP2 = zeros(nQ,nQ); kk = zeros(nQ,nQ); 
MP = zeros(nQ,nQ);

EQmax = sqrt(2*uQmax/(cL*nR*eps0));    % Electric field within aperture
EQ = EQmax .* ones(nQ,nQ);  

zQ = 0;
xQ1 = linspace(xQmin,xQmax,nQ);    
yQ1 = linspace(yQmin,yQmax,nQ);
[xQ, yQ] = meshgrid(xQ1,yQ1);
a = max([ax ay]);
d_RL = a^2/wL;                        % Rayleigh distance                  


% Observation space
xPmin = -xPmax;
yPmin =  -yPmax;

EP = zeros(nP,nP);         % electric field

xP1 = linspace(xPmin,xPmax,nP);     
yP1 = linspace(yPmin,yPmax,nP);
[xP, yP] = meshgrid(xP1,yP1);
indexXY = num*2+1;               % defines X and Y axes  y = 0 and x = 0;

%          optical coordinates
vPx = (pi*(xQmax - xQmin)/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
vPy = (pi*(yQmax - yQmin)/wL) .* yP ./ sqrt(yP.^2 + zP^2);

% =======================================================================
% Specifying the shape of the aperture openings 


EQ(abs(xQ) < bx) = 0;

% =======================================================================

% COMPUATION OF DIFFRACTION INTEGRAL ------------------------------------

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


% GRAPHICS --------------------------------------------------------------
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 14;
set(gca,'fontSize',fs);

%tx = 'x_P (blue)  or y_P (red)  [m]';
%ty = 'irradiance  (a.u.)';

% energy density   uP vs Xp & yP
%subplot(2,2,1);
tx = 'x_P [m] blue - double slit / red - single slit envelope';
ty = 'energy density  u  [W.m^{-2}]';
x = xP(indexXY,:); y = uP(indexXY,:)./uPmax;
plot(x,y ,'b','linewidth',2);
hold on
x = yP(:,indexXY); y = uP(:,indexXY)./uPmax;
plot(x,y,'r','linewidth',2);
xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);


figure(2)      % photograph plot  /  [3D] plot
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.3 0.3]);
pcolor(xP,yP,10.*(uP).^0.3);
shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')

% hold on
% x = xc; y = yc;
% plot(x,y,'y','lineWidth',2);

figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.3 0.3 0.3 0.3]);
surf(xP,yP,10.*(uP).^0.2);
shading interp
axis equal
colormap(jet)
axis off
%set(gcf,'color','k')
xlabel('x');
ylabel('y')


figure(4)
surf(xQ,yQ,real(EQ));


% COMMAND WINDOW OUTPUTS ------------------------------------------------
disp('Parameter summary  [SI units]');
fprintf('wavelength [m]  =  %3.3g \n',wL);
fprintf('nQ  =  %3.3d \n',nQ);
fprintf('nP  =  %3.3d \n',nP);
disp('  ')
disp('Aperture Space');
fprintf('X width [m] =  %3.3e \n',ax);
fprintf('Y width [m]  =  %3.3e \n',ay);

disp('  ')
disp('Observation Space');
fprintf('X width [m] =  %3.3e \n',2*xPmax);
fprintf('Y width [m]  =  %3.3e \n',2*yPmax);
fprintf('distance aperture to observation plane [m]   zP = %3.3e \n',zP);


toc


