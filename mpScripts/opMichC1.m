% opMichC1.m

clear
clc
close all

wL(1) = 588.995e-9;
wL(2) = 589.592e-9;

z1 = 0;
z2 = 3.1e-4;

zp = 0.05;

nP = 5599;
xp = linspace(-0.0050,0.005,nP);

k = 2*pi./wL;


phi(1) = 0;
phi(2) = pi;

R1 = sqrt(zp^2      + xp.^2);
R2 = sqrt((zp-z2)^2 + xp.^2);


E1 = zeros(2,nP);
E2 = zeros(2,nP);
E = zeros(1,nP);

for cn = 1 : 1
  E1(cn,:) = fun(k(1),R1,phi(1));
  E2(cn,:) = fun(k(2),R2,phi(2));
  E = E + E1(cn,:)+ E2(cn,:);
end

intensity = E.*conj(E);
hold on

figure(1)
%xP = xp;
xP = atand(xp./zp);
yP = intensity;
plot(xP,yP,'r')


function E = fun(k,R,phi)
  E = exp(1i*k*R+1i*phi)./R;
end

