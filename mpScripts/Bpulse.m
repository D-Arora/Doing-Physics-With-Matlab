%Bpulse.m


clear
close all
clc

A = 0.3;
v = 3;


Nt = 21;
tMin = 0; tMax = 5;
t = linspace(tMin,tMax,Nt);

N = 501;
zMin = -10; zMax = 10;
z = linspace(zMin,zMax,N);

B = zeros(Nt,N);

for c = 1:Nt
B(c,:) = A^3 ./ (A^2 + (v.*t(c) + z ).^2);
end

B = B./max(max(B));

figure(1)
 ylim([0,1])
for c = 1:Nt
  xP = z; yP = B(c,:);
  plot(xP,yP,'linewidth',2)
  ylim([0,1])
  pause(0.5)
  
end