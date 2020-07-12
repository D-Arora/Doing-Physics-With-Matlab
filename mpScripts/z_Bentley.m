% z_Bentley



clear all
close all
clc

cts = linspace(0,2,2000);

xs = sqrt(1+cts.^2) - 1;


m = 1/sqrt(2);
xD = linspace(0,1.4,500);
yD = m.*xD;
yE = (1/m) .* xD;
figure(1)

plot(xs,cts,'lineWidth',2)
hold on
plot(xD+0.4,yD+1);
plot(xD+0.4,yE+1);
grid on
axis([0 2 0 2]);
axis square
xlabel('x_s')
ylabel('c t_s');

