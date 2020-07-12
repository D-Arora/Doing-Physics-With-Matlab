% op_rings.m

clear all
close all
clc


wL = 6.328e-9;
k = 2*pi/wL;
ES = 1e-3;
a = 0.01;
nR = 101;
zP = 0.1;

xQ = linspace(0,a,nR);


rQP = sqrt(zP^2 + xQ.^2);


EP = (ES ./ rQP) .* exp(-1i * k * rQP);

%EP_sum = sum(EP)
EP_sum = 0;
T = zeros(1,nR);

for c = 1 :nR
    if abs(angle(EP(c))/pi) < 0.5;
       EP_sum = EP_sum + EP(c);
       T(c) = 1;
    end
end


IRR = EP_sum .* conj(EP_sum)

figure(1)
x = xQ;
y = (rQP-zP)/wL;
plot(x,y)

figure(2)
x = xQ;
y = EP;
plot(x,y)

%%
%clear all
%close all
%clc

wL = 6.328e-9;
k = 2*pi/wL;

R = [wL, 2*wL, (3/2)*wL];

z = (1./R) .* exp(-1i * k * R)

zT = sum(z)

z3 = zT .* conj(zT)
