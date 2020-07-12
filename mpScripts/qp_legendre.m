% qp_legendre2.m
% associated legendre polynomials
%polar plots

close all
clear all
clc

L = input('Enter the orbital L (L = 0, 1, 2, 3, ...)  ')
m_L = input('Enter the magnetic quantum number mL (mL = 0, 1, 2, 3, ... L)  ')

tmin = -pi;
tmax = pi;
numt = 511;
dt = (tmax-tmin)/(numt-1);
t = tmin : dt : tmax;

%define cos(t)
z = cos(t);

%definite associated legendre polynomials

Plm = legendre(L,z);

A = Plm(m_L+1,:)/max(abs(Plm(m_L+1,:)));
B = -A;

figure(1)
fs = 10;
xA = A .* sin(t);
yA = A .* cos(t);
xB = B .* sin(t);
yB = B .* cos(t);

tm1 = 'P({\itl},{\itm_l }) =  P(';
tm2 = num2str(L);
tm3 = ', ';
tm4 = num2str(m_L);
tm5 = ')';
tm = [tm1 tm2 tm3 tm4 tm5];

set(gcf,'Units','normalized','position',[0.2,0.2,0.4,0.2]);
subplot(1,2,1);
plot(xA,yA,'k', 'linewidth',3)
hold on
plot(xB,yB,'k', 'linewidth',3)
plot([0 0],[0,1.2],'k','linewidth',2);
h_plot = plot([0 0],[1.18,1.18],'^');
set(h_plot,'MarkerSize',8,'MarkerFaceColor','k', 'MarkerEdgeColor','k');
text(0.15,1.2,'{\itZ} axis','Fontsize',fs);
text(-1.8,1,tm,'Fontsize',fs);
axis equal tight
axis off

xA = A.^2 .* sin(t);
yA = A.^2 .* cos(t);
xB = B.^2 .* sin(t);
yB = B.^2 .* cos(t);

tm1 = 'P({\itl},{\itm_l }) ^2 =  P(';
tm2 = num2str(L);
tm3 = ', ';
tm4 = num2str(m_L);
tm5 = ') ^2';
tm = [tm1 tm2 tm3 tm4 tm5];

%set(gcf,'Units','centimeters','position',[16,6,12,11]);
subplot(1,2,2);
plot(xA,yA,'k', 'linewidth',3)
hold on
plot(xB,yB,'k', 'linewidth',3)
plot([0 0],[0,1.2],'k','linewidth',2);
h_plot = plot([0 0],[1.18,1.18],'^');
set(h_plot,'MarkerSize',8,'MarkerFaceColor','k', 'MarkerEdgeColor','k');
text(0.15,1.2,'{\itZ} axis','Fontsize',fs);
text(-1.8,1,tm,'Fontsize',fs);

axis equal tight
axis off


figure(2)
set(gca,'fontsize',4);
s = sprintf('Associated Legendre Polynomials    L = %.0g     m_L = %.0g ',L,m_L);
subplot(2,1,1)
h_p = polar(t,A);
set(h_p, 'linewidth',2);
hold on
h_p = polar(t,-A);
set(h_p, 'linewidth',2);
grid off
title(s);
subplot(2,1,2)
h_p = polar(t,A.^2,'r');
set(h_p, 'linewidth',2);
hold on
h_p = polar(t,-A.^2,'r');
set(h_p, 'linewidth',2);
grid off
title(s);
axis off


figure(3)    % not complete
Phi=0:pi/10:2*pi;
Theta=-pi/2:pi/10:pi;
[PHI,THETA]=meshgrid(Phi,Theta);
r0=0.5;
R=sin(PHI).*cos(THETA)+r0; % for example
[X,Y,Z]=sph2cart(THETA,PHI,R); %get cartesian values
mesh(X,Y,Z); %display

%surf(X,Y,Z) % colored faces

