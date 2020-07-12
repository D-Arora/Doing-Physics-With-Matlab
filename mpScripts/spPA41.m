% spPA41.m

% Visual Physics Online
% Circuit Experiments
% PA41A.docx

clear all
close all
clc


R = [1500 4700];
emf = [2.2 2.9 4.0 5.1 5.8 7.0 8.1 9.1 10.1 11.8]';

Ic = zeros(length(emf), 2);
Ic = 1000.*emf./R;


%Theoertical values
V = linspace(0,12,200);
I1 = 1e3.*V./R(1);
I2 = 1e3.*V./R(2);

figure(1)
plot(Ic(:,1),emf,'bo','linewidth',2)
hold on
plot(I1,V,'b','linewidth',2)
plot(Ic(:,2),emf,'ro','linewidth',2)
plot(I2,V,'r','linewidth',2)
set(gca,'yLim',[0 12.1]);
grid on
xlabel('current I  [mA]');
ylabel('potential difference V  [V]')

set(gca,'fontsize',14)

%%
emf = 9.8
R1 = 1500
R2 = 4700
Ic = emf / (R1+R2)
V1 = Ic*R1
V2 = Ic*R2

