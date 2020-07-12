% gaussian.m


N = 500;
xmin = -2;
xmax = -xmin;
k = 4;

x = linspace(xmin,xmax,N);

U = -exp(-k.*x.^2);
F = -gradient(U);
figure(1);
set(gcf,'color',[1 1 1]);

subplot(2,1,1);
plot(x,U,'lineWidth',3);
axis off
subplot(2,1,2);
plot(x,F,'lineWidth',3);
axis off

