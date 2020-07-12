% sp6017.m

clear
close all
clc


N =100;
B = 1.25;
R = 135;
L1 = 12e-3; L2 = 6e-3;
w = 20;

A = L1*L2;
w = w*(2*pi/60);
f = w/(2*pi);
T = 1/f;

t = linspace(0,3*T,500);

phi1 = (B*A) .* cos(w*t);
phi2 = (B*A) .* cos(2*w*t);
emf1 = -N*w*B*A*sin(w*t);
emf2 = -N*2*w*B*A*sin(2*w*t);
I1 = emf1/R;
I2 = emf2/R;

fprintf('A = %2.2e \n',A)
fprintf('w = %2.2e \n',w)
fprintf('T = %2.2e \n',T)
fprintf('phi1_max = %2.2e \n',max(phi1))
fprintf('emf1_max = %2.2e \n',max(emf1))
fprintf('I1_max = %2.2e \n',max(I1))

figure(1)
pos = [0.07 0.05 0.32 0.50];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 

xP = t;
tx = 'time  t  [ s ]';

subplot(3,1,1)
yP = phi1;
plot(xP,yP,'b','linewidth',2)
hold on
yP = emf2;
plot(xP,yP,'r','linewidth',2)
xlim([0 9])
set(gca,'xtick',0:1:9)
xlabel(tx)
ylabel('\phi_B  [ T.m^2 ]')
grid on
box on
set(gca,'fontsize',14)

subplot(3,1,2)
yP = emf1;
plot(xP,yP,'b','linewidth',2)
hold on
yP = emf2;
plot(xP,yP,'r','linewidth',2)
xlim([0 9])
set(gca,'xtick',0:1:9)
xlabel(tx)
ylabel('emf  [ V ]')
grid on
box on
set(gca,'fontsize',14)

subplot(3,1,3)
yP = I1;
plot(xP,yP,'b','linewidth',2)
hold on
yP = I2;
plot(xP,yP,'r','linewidth',2)
xlim([0 9])
set(gca,'xtick',0:1:9)
xlabel(tx)
ylabel('I  [ A ]')
grid on
box on
set(gca,'fontsize',14)