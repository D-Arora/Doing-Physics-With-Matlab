% op_rs1_circular_03.m

% Circular aperture - EQ constant in aperature space
% Numerical integration of the Rayleigh-Sommerfield diffraction integral

% Users functiond
%     simpsonxy_coeff.m

% 23 aug 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************
clear all
close all
clc

tic

% INPUTS IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
% Default values for Fraunhoffer diffraction
%    nP = 24*4+1   nQ = 97  wL = 632.8e-9
%    aQ = 1e-4   zP = 1   xPmin = 0   xPmax = 2e-2

nP = 21*4+1;                  % observation points for P  format  integer * 4 + 1
nQ = 97;                      % aperture points for Q  must be ODD

wL = 632.8e-9;               % wavelength [m]
aQ = 10 * wL;                   % radius of circular aperture [m]
zP = 1;                      % z distance between aperture plan and observation plane [m]

xPmin = 0;                   % defintes observation screen size
xPmax = 0.2;

PQ = 1e-3;                   % incident power

% Setup   SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
ik = 1i*2*pi/wL;             % j k

EP = zeros(nP,nP);           % intialize matrices
WP = zeros(nP,nP);
WPs = zeros(nP,nP);
area = zeros(nP,1);

cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
n = 1;                  % refractive index
SQ = pi * aQ^2;         % area of aperture
WQ = PQ / SQ;           % uniform irradiance of aperture
EQ = sqrt(2*WQ / (cL * n * eps0));    % uniform electric field within aperture

% Q aperture space   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
tQmin = 0;                    % angle thetaQ  [rad]
tQmax = 2*pi;

rQmin = 0;                    % radius rhoQ  [m]
rQmax = aQ;           


tQ1 = linspace(tQmin,tQmax,nQ);      % polar coordinates for Q
rQ1 = linspace(rQmin,rQmax,nQ);

ht = tQ1(2) - tQ1(1);        % dx for integraltion
hr = rQ1(2) - rQ1(1);        % dy for inegration

KQ = (ht/3)*(hr/3)* EQ;       % constant for integration

[tQ, rQ] = meshgrid(tQ1,rQ1);

costQ = cos(tQ);
sintQ = sin(tQ);

xQ = rQ .* costQ;                    % cartesian coordinates for Q         
yQ = rQ .* sintQ;
zQ = 0;

unit = ones(nQ,nQ);                  % unit matrix 
K = KQ .* ones(nQ,nQ);

S = simpsonxy_coeff(nQ);             % calculates Simpson 2D coefficients
MQ = K .* rQ .* S;                        % forms matrix

% P observation space   PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
tPmin = 0;                    % angle thetaP  [rad]
tPmax = 2*pi;

rPmin = 0;                    % radius rhoP  [m]
rPmax = xPmax;           

tP1 = linspace(tPmin,tPmax,nP);      % polar coordinates for P
rP1 = linspace(rPmin,rPmax,nP);
[tP, rP] = meshgrid(tP1,rP1);

costP = cos(tP);
sintP = sin(tP);

xP = rP .* costP;                    % cartesian coordinates for P         
yP = rP .* sintP;

% radial optical coordinate v along x and y axes 
indexyP = (nP+3)/4;                 % defines y=axis  theta = 90 deg
vPx = (2*pi*aQ/wL) .* xP(:,1) ./ sqrt(xP(:,1).^2 + zP^2);
vPy = (2*pi*aQ/wL) .* yP(:,indexyP) ./ sqrt(yP(:,indexyP).^2 + zP^2);

vPxx = (2*pi*aQ/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
vPyy = (2*pi*aQ/wL) .* yP ./ sqrt(yP.^2 + zP^2);

% Calculates integral for irradiance CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

% calculates sum to approximate integral
for c1 = 1 : nP
    for c2 = 1 :nP
       rPQ = fn_distancePQ(xP(c1),yP(c2),zP,xQ,yQ,zQ);
       rPQ3 = rPQ .* rPQ .* rPQ;
       MP1 = zP .* (ik .* rPQ - unit);
       MP2 = exp(ik .* rPQ) ./ rPQ3;
       MP = MP1 .* MP2;
       EP(c1,c2) = sum(sum(MQ .* MP));
    end 
end

% calculation of irradiance (intensity or energy density)
WP = (cL * n * eps0 / 2) .* real(conj(EP) .* EP);
WPmax = max(max(WP));
WP = WP ./ WPmax;
WPs = 10 .* log10(WP);


% Energy enclosed within a circle of prescribed radius -------------------
rho = sqrt(vPxx.^2 + vPyy.^2);   % magnitude of radial coordinate 

tMin = 0;                 % angle theta
tMax = 2*pi;

rhoMin = min(min(rho));   % radius rho
rhoMax = max(max(rho));

fn = rho .* WP;            % integration - polar coordinates  rho*drho* dphi 


for c = 1 : nP
radius = vPx(c);                % prescribed radius
fn(rho > radius) = 0;           % function = 0 outside prescribed radius 
area(c) = simpson2d(fn,rhoMin,rhoMax,tMin,tMax);   % performs integration
fn = rho .* WP;    
end
%area = 100 .* area / max(area);   % normalize area to 100 %


% GRAPHICS GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.4 0.4]);
fs = 12;
tx = 'x_P (blue)  or y_P (red)  [m]';
ty = 'irradiance  (a.u.)';

subplot(2,2,1);              % IRR vs xP or yP
plot(xP(:,1), WP(:,1),'linewidth',2);
hold on
plot(-xP(:,1), WP(:,1),'linewidth',2);
plot(yP(:,indexyP), WP(:,indexyP),'r','linewidth',2);
plot(-yP(:,indexyP), WP(:,indexyP),'r','linewidth',2);
grid on
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
%axis([0 xPmax -0.1 1.1]);
%set(gca,'xlim',[0 xPmax]);
set(gca,'Ylim',[-0.1 1.1]);


subplot(2,2,3);               %log scale IRR vs xP or yP
plot(xP(:,1), WPs(:,1),'linewidth',2);
hold on
plot(-xP(:,1), WPs(:,1),'linewidth',2);
plot(yP(:,indexyP), WPs(:,indexyP),'r','linewidth',2);
plot(-yP(:,indexyP), WPs(:,indexyP),'r','linewidth',2);
grid on
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
%axis([0 xPmax -0.1 1.1]);
%set(gca,'xlim',[0 xPmax]);
set(gca,'Ylim',[-75 5]);


subplot(2,2,2);                % IRR vs vPx   vPy
tx = 'v_{Px} (blue)  or y_{Py} (red)';
plot(vPx, WP(:,1),'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
%set(gca,'Ylim',[-0.1 1.1]);
hold on
plot(-vPx, WP(:,1),'linewidth',2);
plot(vPy, WP(:,indexyP),'r','linewidth',2);
plot(-vPy, WP(:,indexyP),'r','linewidth',2);
%set(gca,'Ylim',[-75 5]);
%axis([0 xPmax -90 10]);
%set(gca,'xlim',[0 xPmax]);
grid on

subplot(2,2,4);              % log(IRR) vs vPx   vPy
plot(vPx, WPs(:,1),'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
%set(gca,'Ylim',[-0.1 1.1]);
hold on
plot(-vPx, WPs(:,1),'linewidth',2);
plot(vPy, WPs(:,indexyP),'r','linewidth',2);
plot(-vPy, WPs(:,indexyP),'r','linewidth',2);
%axis([0 xPmax -90 10]);
%set(gca,'xlim',[0 xPmax]);
set(gca,'Ylim',[-75 5]);
grid on

figure(2)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.3]);
pcolor(xP,yP,10.*(WP).^0.3);

shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')


figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.4 0.4]);
fs = 12;
tx = 'vP';
ty = 'Energy enclosed (a.u.)';
plot(vPx, area,'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);

toc



