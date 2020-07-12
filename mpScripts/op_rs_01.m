% op_rs_01.m

% Rectangular aperture
% Numerical integration of the Rayleigh-Sommerfield diffraction integral of
% the first kind
% Uses cartesain coordinates / S.I. units
% Calculate of energy enclosed in circles
% Calls   fn_distancePQ.m
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]


% Uses functions
%     simpson2d.m  distancePQ.m

% 12 oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************
clear all
close all
clc

tic


% INPUT PARAMETERS -------------------------------------------------------
num = 30;             % number for observation space
nP = num*4+1;         % observation points for P  format  integer * 4 + 1
nQ = 49;              % aperture points for Q  must be ODD

wL = 632.8e-9;        % wavelength [m]
 
ax = 2e-4;             % width - aperture space  [m]
ay = 2e-4;  

uQmax = 1e-3;                % incident energy density  [W.m^-2]

xPmax = 2e-2;         % half width - observation space   [m] 
yPmax = 2e-2;
zP = 1;                % distance: aperture plan to observation plane [m]


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

uQ = (cL*nR*eps0/2) .* EQ .*EQ;       % energy density within aperture
UQ = simpson2d(uQ,xQmin,xQmax,yQmin,yQmax);    % energy transmitted from aperture per s [J/s]
UQtheory = uQmax * (xQmax-xQmin) * (yQmax - yQmin);

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


%          first min xP and yP values
thetax1 = asin(wL/ax); thetay1 = asin(wL/ay);
x1 = zP * tan(thetax1); y1 = zP * tan(thetay1);
 


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

UP = simpson2d(uP,xPmin,xPmax,yPmin,yPmax);


% percentage of energy enclosed within a circle of radius of 1st min in x direction
Rc = 3.831*2;
uPenclosed = uP;
uPenclosed(vPx > Rc) = 0;
UPenclosed = simpson2d(uPenclosed,xPmin,xPmax,yPmin,yPmax);
UPenclosed_perc = 100*UPenclosed/UP;

%theta_c = asin(Rc/(k*ax));     % ax effective diameter of circular aperture
%xPc = zP * tan(theta_c);

% Draw circle
nc = 500;
theta = linspace(0,2*pi,nc);
xc = x1 .* cos(theta); yc = x1 .* sin(theta);


% GRAPHICS --------------------------------------------------------------
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 12;
%tx = 'x_P (blue)  or y_P (red)  [m]';
%ty = 'irradiance  (a.u.)';

% energy density   uP vs Xp & yP
subplot(2,2,1);
tx = 'x_P (blue)  or y_P (red)  [m]';
ty = 'energy density  u  [W.m^{-2}]';
x = xP(indexXY,:); y = uP(indexXY,:);
plot(x,y ,'linewidth',2);
hold on
x = yP(:,indexXY); y = uP(:,indexXY);
plot(x,y,'linewidth',2);
xlabel(tx);   ylabel(ty);

% energy density   uPdB  vs xP & yP
subplot(2,2,3);   
tx = 'x_P (blue)  or y_P (red)  [m]';
ty = 'energy density  u  [dB]';
x = xP(indexXY,:);
y = uPdB(indexXY,:);
plot(x,y ,'linewidth',2);
hold on
x = yP(:,indexXY);
y = uPdB(:,indexXY);
plot(x,y,'linewidth',2)
xlabel(tx);   ylabel(ty);

% energy density   uP vs vPx & vPy
subplot(2,2,2);
tx = 'optical coordinate  v_P / \pi';
ty = 'energy density  u  [a.u.]';
x = vPx(indexXY,indexXY:end)./pi;
y = uP(indexXY,indexXY:end)./uPmax;
plot(x,y ,'linewidth',2);
xlabel(tx);   ylabel(ty);
grid on

% energy density   uPdB  vs xP & yP
subplot(2,2,4);   
tx = 'optical coordinate   v_P / \pi';
ty = 'energy density  u  [dB]';
x = vPx(indexXY,indexXY:end)./pi;
y = uPdB(indexXY,indexXY:end)./uPmax;
plot(x,y ,'linewidth',2);
xlabel(tx);   ylabel(ty);
grid on

figure(2)      % photograph plot  /  [3D] plot
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.3 0.3]);
pcolor(xP,yP,10.*(uP).^0.3);
shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')

hold on
x = xc; y = yc;
plot(x,y,'y','lineWidth',2);

figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.3 0.3 0.3 0.3]);
surf(xP,yP,10.*(uP).^0.3);
shading interp
axis equal
colormap(jet)
axis off
%set(gcf,'color','k')
xlabel('x');
ylabel('y')


% COMMAND WINDOW OUTPUTS ------------------------------------------------
disp('Parameter summary  [SI units]');
fprintf('wavelength  =  %3.3g \n',wL);
fprintf('nQ  =  %3.3d \n',nQ);
fprintf('nP  =  %3.3d \n',nP);
disp('  ')
disp('Aperature Space');
fprintf('X width  =  %3.3e \n',ax);
fprintf('Y width  =  %3.3e \n',ay);
fprintf('energy density uQmax  =  %3.3e \n',uQmax);
fprintf('energy from aperture / s   UQ(theory) = %3.3e \n',UQtheory);
fprintf('energy from aperture / s   UQ(calculated) = %3.3e \n',UQ);
disp('  ')
disp('Observation Space');
fprintf('X width  =  %3.3e \n',2*xPmax);
fprintf('Y width  =  %3.3e \n',2*yPmax);
fprintf('distance aperture to observation plane   zP = %3.3e \n',zP);
fprintf('Fraunhofer: position of 1st min in X direction  =  %3.3e \n',x1);
fprintf('Fraunhofer: position of 1st min in Y direction  =  %3.3e \n',y1);
disp('   ');
fprintf('energy to aperture / s   UP = %3.3e \n',UP);
disp('   ');
fprintf('pecentage energy enclosed within circle of radius (vPx = 3.831)  = %3.1f \n',UPenclosed_perc);
disp('   ');


toc


