% h_lithium
% data saved from qn_hydrogen.m
% Ian Cooper  Scholl of Physics, University of Sydney

clc
clear all 
close all

load('qp_hL')

figure(1)
set(gcf,'units','normalized');
set(gcf,'position',[0.1 0.1 0.6 0.8]);
set(gca,'fontsize',14);

xmax = 2;
subplot(4,1,1)

x = 1e9 .* r; y = (0.5 .*(p1s1 - p1s2)).^2;
plot(x,y,'linewidth',2);
ylabel('inner shell','fontsize',14);
h_title = title('Li atom:   Probability density functions ');
h_t = text(0.3,2e10,'blue 1s / red s-states / black p-states / magenta d-states / green f-states');  
set(h_t,'fontsize',14);
set(h_title,'fontsize',14);
set(gca,'Xlim',[0 xmax]);
set(gca,'fontsize',14);

subplot(4,1,2)
y = p2s.^2;
plot(x,y,'r','linewidth',2);
hold on
y = p2p.^2;
plot(x,y,'k','linewidth',2);
ylabel('shell n = 2','fontsize',14);
set(gca,'Xlim',[0 xmax]);
set(gca,'fontsize',14);


subplot(4,1,3)
y = p3s.^2;
plot(x,y,'r','linewidth',2);
hold on
y = p3p.^2;
plot(x,y,'k','linewidth',2);
y = p3d.^2;
plot(x,y,'m','linewidth',2);
ylabel('shell n = 3','fontsize',14);
set(gca,'Xlim',[0 xmax]);
set(gca,'fontsize',14);

subplot(4,1,4)
y = p4s.^2;
plot(x,y,'r','linewidth',2);
hold on
y = p4p.^2;
plot(x,y,'k','linewidth',2);
y = p4d.^2;
plot(x,y,'m','linewidth',2);
y = p4f.^2;
plot(x,y,'g','linewidth',2);
set(gca,'Xlim',[0 xmax]);

ylabel('shell n = 4','fontsize',14);
xlabel('radial position  r  (nm)','fontsize',14)
set(gca,'fontsize',14);

