% qp_balmer.m
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm
% Plot of the line spectrum for the Balmer series


%%
clear all
close all
clc

c = 3e8;
G = 6.67e-11;
NA = 6.02e23;
k = 1.38e-23;
e = 1.602e-19;
sigma = 5.67e-8;
u0 = 4*pi*1e-7;
h = 6.625e-34;
hP = 6.63e-34;
me = 9.109e-31;
mp = 1.673e-27;
mn = 1.675e-31;
g = 9.81;

RH = 10967757;




%% hydrogen spectrum - Balmer series
%close all
%clear all
%clc



RH = 10967759;

ni = [3:15];

nf = 2 .* ones(length(ni),1)';


RwL = RH .* (1./nf.^2 - 1./ni.^2);

wL = (1./ RwL)*1e9;

num = length(wL);

y1 = zeros(1,num);
y2 = ones(1,num);


figure(1)
fs = 14;
set(gcf,'color',[1 1 1]);
set(gcf,'units','normalized'); 
set(gcf,'position',[0.1 0.1 0.6 0.1]);
set(gca,'fontsize',fs);
for cn = 1 : num
    h_plot = plot([wL(cn) wL(cn)],[y1(cn) y2(cn)],'k','linewidth',2);
    hold on
end
xlabel('wavelength   (nm)');
set(gca,'YColor',[1 1 1]);

%%   series limit  Rydberg  nf = 1
wL_min = 1e9*1/RH;
dE = h * c * RH;
dE = dE / e;


%%  hydrogen spectral lines

n = 1:6;
E = -13.6 ./ n.^2;

num = length(E);
E(6) = 0;
x1 = zeros(1,num); x2 = ones(1,num);

figure(2)
fs = 14;
set(gcf,'color',[1 1 1]);
set(gcf,'units','normalized'); 
set(gcf,'position',[0.1 0.1 0.6 0.6]);
set(gca,'fontsize',fs);
for cn = 1 : num-1
    h_plot = plot([x1(cn) x2(cn)],[E(cn) E(cn)],'k','linewidth',2);
    hold on
end
ylabel('energy   (eV)');
set(gca,'XColor',[1 1 1]);
h_plot = plot([x1(6) x2(6)],[E(6) E(6)],'b','linewidth',3);

%% spiral
%clear all
%close all
%clc

nt = 1200;
tmax = 6;
T = 1;
Ao = 1;
b = 1.5*1/tmax;
w = 2*pi/T;

t = linspace(0,tmax,nt);

A = Ao .* exp(-b * t);

x = A .* sin(w*t);
y = A .* cos(w*t);

figure(3)

set(gcf,'color',[1 1 1]);
set(gcf,'units','normalized'); 
set(gcf,'position',[0.1 0.1 0.6 0.6]);

plot(x,y,'linewidth',3);
axis equal
axis off







