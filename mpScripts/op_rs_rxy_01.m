% op_rs_rxy_01.m

% Rectangular aperture
% Numerical integration of the Rayleigh-Sommerfield diffraction integral
% Uses Cartesian coordinates
% 
% Users functiond
%     simpsonxy_coeff.m     distancePQ.m

% 07 sep 2014
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

nP = 99;                  % observation points for P  must be ODD
nQ = 99;                  % aperture points for Q  must be ODD

wL = 632.8e-9;            % wavelength [m]

aQx = 2e-4;               % defines aperture size  [m]
aQy = 2e-4; 

aP = 2e-2;               % defines size of observation space  [m]
%aP = 10*wL;

zQ = 0; zP = 1;           % z distance between aperture plan and observation plane [m]
%zP = 100 * wL;

xQmin = -aQx/2; xQmax = +aQx/2;     % aperture screen size
yQmin = -aQy/2; yQmax = +aQy/2;
 
xPmin = -aP; xPmax = aP;             % observation screen size
yPmin = -aP; yPmax = aP;

% Setup   SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
ik = 1i * (2*pi/wL);

EP = zeros(nP,nP);           % intialize matrices
WP = zeros(nP,nP);
WPs = zeros(nP,nP);
unit = ones(nQ,nQ);          % unit matrix 

% aperture
xQ1 = linspace(xQmin,xQmax,nQ);      
yQ1 = linspace(yQmin,yQmax,nQ);
[xQ, yQ] = meshgrid(xQ1,yQ1);

% observation screen
xP1 = linspace(xPmin,xPmax,nP);      
yP1 = linspace(yPmin,yPmax,nP);
[xP, yP] = meshgrid(xP1,yP1);

% calculates Simpson 2D coefficients
S = simpsonxy_coeff(nQ);             

% radial optical coordinate v along x and y axes 
indexyP = floor(nP/2) + 1;                 % defines y-axis 
vPx = (pi*aQx/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
vPy = (pi*aQy/wL) .* yP ./ sqrt(yP.^2 + zP^2);

% Calculation of irradiance (evlaution of integral)CCCCCCCCCCCCCCCCCCCCCCC
for c1 = 1 : nP
    for c2 = 1 :nP
       rPQ = fn_distancePQ(xP(c1,c2),yP(c1,c2),zP,xQ,yQ,zQ);
       rPQ3 = rPQ .* rPQ .* rPQ;
       MP1 = exp(ik .* rPQ) ./ rPQ3;
       MP2 = (ik .* rPQ - unit);
       MP = S .* MP1 .* MP2;
       EP(c1,c2) = sum(sum(MP));
    end 
end

% Calculation of irradiance (intensity or energy density)
WP = conj(EP) .* EP;
WPmax = max(max(WP));
WP = WP ./ WPmax;
WPs = 10 .* log10(WP);

% GRAPHICS GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.4 0.6]);
fs = 12;
tx = 'x-axis: x_P (blue)    y-axis: y_P (red)  [m]';
ty = 'irradiance  (a.u.)';

subplot(2,2,1);              % IRR vs xP or yP
plot(xP(indexyP,:), WP(indexyP,:),'b','linewidth',2);
hold on
plot(yP(:,indexyP), WP(:,indexyP),'r','linewidth',1);
grid on
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
set(gca,'xlim',[-aP aP]);
set(gca,'Ylim',[-0.1 1.1]);
set(gca,'fontsize',fs);

subplot(2,2,3);               %log scale IRR vs xP or yP
plot(xP(indexyP,:), WPs(indexyP,:),'linewidth',2);
hold on
plot(yP(:,indexyP), WPs(:,indexyP),'r','linewidth',1);
grid on
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
set(gca,'xlim',[-aP aP]);
set(gca,'Ylim',[-75 5]);
set(gca,'fontsize',fs);

subplot(2,2,2);                % IRR vs vPx   vPy
tx = 'x-axis  v_{P} / \pi ';
plot(vPx(indexyP,:)./pi, WP(indexyP,:),'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
hold on
plot(vPy(:,indexyP)./pi, WP(:,indexyP),'r','linewidth',1);
set(gca,'Ylim',[-0.1 1.1]);
grid on
set(gca,'fontsize',fs);

subplot(2,2,4);              % log(IRR) vs vPx   vPy
plot(vPx(indexyP,:)./pi, WPs(indexyP,:),'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
hold on
plot(vPy(:,indexyP)./pi, WPs(:,indexyP),'r','linewidth',1);
set(gca,'Ylim',[-75 5]);
grid on
set(gca,'fontsize',fs);

figure(2)
sf = 0.2;
set(gcf,'Units','normalized');
set(gcf,'Position',[0.5 0.5 0.3 0.3]);
pcolor(xP,yP,100.*(WP).^sf);

shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')


toc



