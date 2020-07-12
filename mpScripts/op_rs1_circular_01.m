% rs1_circular_01.m

% Circular aperture - EQ constant in aperature space
% Numerical integration of the Rayleigh-Sommerfield diffraction integral 1

% Users functiond
%     simpsonxy_coeff.m   fn_distancePQ.m

% 29 aug 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************
clear all
close all
clc

tic

% INPUTS IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
nP = 99;                      % observation points for P
nQ = 55;                      % aperture points for Q  must be ODD

wL = 632.8e-9;               % wavelength [m]
aQ = 10*wL;                   % radius of circular aperture [m]
zP = 2;                    % z distance between aperture plan and observation plane [m]

xPmin = 0;
xPmax = 2e-1;

% Setup   SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
ik = 1i*2*pi/wL;             % j k

% Q aperture space   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
tQmin = 0;                    % angle thetaQ  [rad]
tQmax = 2*pi;

rQmin = 0;                    % radius rhoQ  [m]
rQmax = aQ;           

tQ1 = linspace(tQmin,tQmax,nQ);      % polar coordinates for Q
rQ1 = linspace(rQmin,rQmax,nQ);
[tQ, rQ] = meshgrid(tQ1,rQ1);

costQ = cos(tQ);
sintQ = sin(tQ);

xQ = rQ .* costQ;                    % cartesian coordinates for Q         
yQ = rQ .* sintQ;
zQ = 0;

S = simpsonxy_coeff(nQ);
MQ = rQ .* S;

% P observation space   PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
%tPmin = 0;                    % angle thetaP  [rad]
%tPmax = 2*pi;

%rPmin = 0;                    % radius rhoP  [m]
%rPmax = 0.1;  

xP = linspace(xPmin,xPmax,nP);   % xP coordinates
yP = 0;
unit = ones(nQ,nQ);
EP = zeros(nP,1);

for cP = 1 : nP
    rPQ = fn_distancePQ(xP(cP),yP,zP,xQ,yQ,zQ);
    rPQ3 = rPQ .* rPQ .* rPQ;
    MP1 = zP .* (ik .* rPQ - unit);
    MP2 = exp(ik .* rPQ) ./ rPQ3;
    MP = MP1 .* MP2;
    EP(cP) = sum(sum(MQ .* MP));
end


WP = conj(EP) .* EP;
WPmax = max(WP);
WP = WP ./ WPmax;
WPs = 10 .* log10(WP);


r = xP;
t = linspace(0,2*pi,nP);

xx = zeros(nP,nP);
yy = zeros(nP,nP);

for c = 1: nP
    xx(c,:) = r .* cos(t(c));
    yy(c,:) = r .* sin(t(c));
end

WPP = meshgrid(WP,WP);

vP = (2*pi*aQ/wL) .* xP ./ sqrt(xP.^2 + zP^2);
J1 = besselj(1,vP);           % Bessel function of the first kind
IRR = (J1 ./ vP).^2;          % irradiance      
IRR = IRR ./ max(IRR);
IRRs = 10 .* log10(IRR);

figure(1)
fs = 12;
tx = 'xP   (m)';
ty = 'irradiance  (a.u.)';

subplot(2,1,1);
plot(xP, WP,'linewidth',2);
grid on
hold on
plot(xP, IRR,'r','linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
axis([0 xPmax -0.1 1.1]);
set(gca,'xlim',[0 xPmax]);
set(gca,'Ylim',[-0.1 1.1]);
legend('FD', 'RS1');

subplot(2,1,2);
ty = 'irradiance  (dB)';
plot(xP, WPs,'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
hold on
plot(xP, IRRs,'r','linewidth',2);
axis([0 xPmax -90 10]);
set(gca,'xlim',[0 xPmax]);
set(gca,'Ylim',[-90 10]);
grid on


figure(2)
fs = 12;
tx = 'vP';
ty = 'irradiance  (a.u.)';
subplot(2,1,1);
plot(vP, WP,'linewidth',2);
set(gca,'Ylim',[-0.1 1.1]);
grid on
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
%axis([0 xPmax -0.1 1.1]);
%set(gca,'xlim',[0 xPmax]);
%set(gca,'Ylim',[-0.1 1.1]);
legend('RS1')

subplot(2,1,2);
ty = 'irradiance  (dB)';
plot(vP, WPs,'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
%axis([0 xPmax -90 10]);
%set(gca,'xlim',[0 xPmax]);
%set(gca,'Ylim',[-90 10]);
grid on

figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.2 0.3]);
pcolor(xx,yy,10.*(WPP).^0.3);
shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')

toc



