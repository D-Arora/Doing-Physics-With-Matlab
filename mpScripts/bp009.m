% bp009.m

clear
close all
clc

Gbasal = 85;
GdotBasal = 2;
G1 = linspace(Gbasal,200,501) - Gbasal;
G2 = -linspace(50,Gbasal,501)  + Gbasal;
G = linspace(25,225,501);
Gmax = 125;
Gmin = 70;

k = 1/(Gmax - Gmin);
k = 0.1;
%Gdot1 = GdotBasal .*(1 - (1 - exp(-G1./tau)));
%Gdot2 = GdotBasal .*(1+(exp(-G2./tau)-1));

Gdot1 = G1.^2./(sqrt(Gbasal^2) + G1.^2);
Gdot2 = G2.^2./(sqrt(Gbasal^2) + G2.^2);


Gdot = 4./(1+exp(k*(G-Gbasal)));

figure(1)
plot(G,Gdot)
grid on

%hold on
%plot(G2,Gdot2)
