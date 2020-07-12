% tau_nmh_inf_plot

clear all
close all
clc

global Vr

Vr = -70;  T = 6.3;
V = linspace(-110,50,1000);

[Tn Tm Th] = tau(V,T);

[ n_inf m_inf h_inf ] = N_inf(V,T);

[ n_0 m_0 h_0 ] = N_0(T);

figure(2)
set(gcf,'unit','normalized')
set(gcf,'position',[0.1 0.1 0.58 0.85])
fs = 14;
% n m h plots
subplot(3,2,1)
plot(V,n_inf,'lineWidth',2)
set(gca,'fontsize',fs);
hold on
plot(Vr,n_0,'o','markersize',6)
xlabel('membrane potential   V   (mV)')
ylabel('K   n _{inf}')
axis([-100 50 0 1]);

subplot(3,2,3)
hold on
plot(Vr,m_0,'o')
set(gca,'fontsize',fs);
plot(V,m_inf,'lineWidth',2)
xlabel('membrane potential   V   (mV)')
ylabel('Na   m _{inf}')
axis([-100 50 0 1]);

subplot(3,2,5)
plot(V,h_inf,'lineWidth',2)
set(gca,'fontsize',fs);
hold on
plot(Vr,h_0,'o')
xlabel('membrane potential   V   (mV)')
ylabel('Na   h _{inf}')
axis([-100 50 0 1]);

% tau plots
subplot(3,2,2)
plot(V,Tn,'lineWidth',2)
set(gca,'fontsize',fs);
xlabel('membrane potential   V   (mV)')
ylabel('K+   tau _n   (ms)')
axis([-100 50 0 10]);

subplot(3,2,4)
plot(V,Tm,'lineWidth',2)
set(gca,'fontsize',fs);
xlabel('membrane potential   V   (mV)')
ylabel('Na+   tau _m   (ms)')
axis([-100 50 0 1]);

subplot(3,2,6)
plot(V,Th,'lineWidth',2)
set(gca,'fontsize',fs);
xlabel('membrane potential   V   (mV)')
ylabel('Na   tau _h   (ms)')
axis([-100 50 0 10]);