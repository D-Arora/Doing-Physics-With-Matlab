% ztest.m


clear x
clear y
clear yE
close all

y = Id(40:79);
x = 0:length(y)-1;

tH = 8;
k = log(2)/tH;
yE = y(1).*exp(-k*x);

figure (99)

plot(x,y,'o')
hold on
plot(x,yE)

% Exponetial Growth
y = Id(5:15);
x = 0:length(y)-1;

tH = 3;
k = log(2)/tH;
yE = y(1).*exp(k*x);

figure (88)

plot(x,y,'o')
hold on
plot(x,yE)