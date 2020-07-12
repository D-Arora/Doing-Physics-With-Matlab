% mec_tt_ball.m

tb = [0 1.81 2.57 3.24 3.78 4.2 4.66 5.09 5.46 5.78 6.18 6.5 6.77 7.1 7.45 7.7 7.99 8.25 8.52 8.81 9.03 9.30];

yb = [0 .28 .58 .89 1.23 1.49 1.82 2.14 2.47 2.74 3.06 3.4 3.65 3.97 4.33 4.59 4.89 5.22 5.53 5.85 6.13 6.42];

tb = tb .* (0.71/9.57);
yb = yb .* (2.5/7.6);


figure(10)
fs = 9;
col ='bo';
set(gcf,'units','normalized','position',[0.6 0.2 0.18 0.22]);
tx = 'time  t   [s]';
ty = 'displacement  x  [m]';
xP = tb;
yP = yb;

plot(xP,yP,col);
grid on

hold on

xP = t;
yP = x;

plot(xP,yP,'r');

xlabel(tx); ylabel(ty);