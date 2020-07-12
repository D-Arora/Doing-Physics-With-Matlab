% op_rs_xy_04.m

% Rectangular aperture space: miscellaneous shaped aperture openings
% Cross shaped aperture
% Numerical integration of the Rayleigh-Sommerfield diffraction integral of
% the first kind
% Uses cartesain coordinates / S.I. units

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
clc

tic


% INPUT PARAMETERS -------------------------------------------------------
num = 80;            % number for observation space
nP = num*4+1;         % observation points for P  format  integer * 4 + 1
nQ = 59;              % aperture points for Q  must be ODD

wL = 650e-9;        % wavelength [m]
 

% aperture space
ax = 20*wL;             % half widths  [m]
ay = 20*wL;



uQmax = 1e-3;                % incident energy density  [W.m^-2]

xPmax = 2000*wL;         % half width - observation space   [m] 
yPmax = 2000*wL;
%zP = 1e-3;                % distance: aperture plan to observation plane [m]
zP = 6000*wL;

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
xQmin = -ax;  xQmax = ax;
yQmin =  -ay; yQmax =  ay;

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
f = 0.4;               % fraction
numQ = round(f*nQ);
EQ(1:numQ,nQ-numQ:nQ) = 0;
EQ(1:numQ,1:numQ) = 0;
EQ(nQ-numQ:nQ,1:numQ) = 0;
EQ(nQ-numQ:nQ,nQ-numQ:nQ) = 0;

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
set(gcf,'Position',[0.1 0.1 0.4 0.4]);
fs = 14;

tx = 'x_P /  \lambda   or   y_P /  \lambda';
ty = 'irradiance  (a.u.)';

x = xP(indexXY,:); y = uP(indexXY,:)./uPmax;
x = x ./wL;
plot(x,y ,'b','linewidth',2);
%hold on
%x = yP(:,indexXY); y = uP(:,indexXY)./uPmax;
%plot(x,y,'r','linewidth',2);
xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
axis([xPmin./wL xPmax./wL 0 1]);
set(gca,'fontSize',fs);


figure(2)      % photograph plot  /  [3D] plot
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.5 0.5]);
x = xP ./ wL; y = yP ./wL;
pcolor(x,y,10.*(uP).^0.2);
shading interp
axis equal
colormap(gray)
colormap(jet)
axis on
%set(gcf,'color','k')
axis([xPmin xPmax yPmin yPmax]./wL)
tx = 'x_Q /  \lambda';
ty = 'y_Q /  \lambda';
set(gca,'fontsize',fs);
xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
%axis square
set(gca,'Xtick',[-2000:1000:2000]);
set(gca,'Ytick',[-2000: 1000: 2000]);


figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.3 0.3 0.5 0.5]);
x = xP ./ wL; y = yP ./wL;
%surf(x,y,10.*(uP).^0.05);
surf(x,y,(uP).^0.01);
%shading interp
shading flat
%axis equal
colormap(jet)
axis off
%set(gcf,'color','k')
xlabel('x');
ylabel('y')


figure(4)
x = xQ./wL; y = yQ ./ wL;
pcolor(x,y,real(EQ));
shading flat
tx = 'x_Q /  \lambda';
ty = 'y_Q /  \lambda';
set(gca,'fontsize',fs);
xlabel(tx,'fontsize',fs);   ylabel(ty,'fontsize',fs);
axis square

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


