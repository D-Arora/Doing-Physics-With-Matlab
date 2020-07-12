% wm_WE2D.m


close all
clear all
clc

tic

N = 1000;
Nt = 2000;
T = 0.01;  f = 1 / T;
%f = 0.1; T = 1/f;

u = zeros(N,N,Nt);
v = 20.* ones(N,N);

xMin = 0; xMax = 1;
x = linspace(xMin,xMax,N);
dx = x(2) - x(1);

k = 1/2;
dt = dx * k / max(max(v));
tMin = 0; tMax = Nt * dt;
%tMax = 10*T;
t = linspace(tMin,tMax,Nt);
%K2 = (v .* dt ./ dx) .^2;
%K1 = 2.*(1-K2);
%K = (max(max(v)) * dt / dx)^2;
u(:,1,1) = 0.1 .* sin(2*pi*t(1)/T) .* ones(N,1);
u(:,1,2) = 0.1 .* sin(2*pi*t(2)/T) .* ones(N,1);

v(:,N/2:N) = 10;

k = (dt/dx) .* v;
%K = 0.5;

for nT = 3 : Nt
for nC = 2 : N-1
for nR = 2 : N-1
 % u(nR,1,nT) = sin(2*pi*t(nT)/T); %.* ones(N,1);
  %u(nR,nC,nT) = K1(nR,nC) * u(nR,nC,nT-1) - u(nR,nC, nT-2) ...
   %             + K2(nR,nC) * (u(nR+1,nC,nT-1) + u(nR-1,nC,nT-1) ...
   %             +       u(nR,nC+1,nT-1) + u(nR,nC-1,nT-1)      );
   u(nR,1,nT) = 0.1 .* sin(2*pi*t(nT)/T); 
   
   u(nR,nC,nT) = 2*u(nR,nC,nT-1) - u(nR,nC, nT-2) ...
           + k(nR,nC)^2 * (  u(nR+1,nC,nT-1) - 2*u(nR,nC,nT-1) + u(nR-1,nC,nT-1) ...
           +          u(nR,nC+1,nT-1) - 2*u(nR,nC,nT-1) + u(nR,nC-1,nT-1)  );
   
   u(nR,N,nT) = u(nR,N,nT-2);
   u(N,nC,nT) = u(N,nC,nT-2);
   u(1,nC,nT) = u(1,nC,nT-2);
   
 end   % row nR
 end   % col nC
 end   % time nT

 % ======================================================================

figure(1)

for c = 1 : 15 : Nt
uXY = u(:,:,c);
uXY(uXY>1) = 1;
uXY(uXY<-1) = -1;

pcolor(x,x,uXY);
%colormap(gray(2));

%contourf(uXY);
shading flat;
colormap(winter);
%colormap(gray(56));
pause(0.01);

end


% =======================================================================
% figure(2)
% xP = x; yP = uXY(50,:);
% plot(xP, yP);
% % % for c = 1 :Nt
% % %  
% % %     uXY = rand(N,N);
% % %     u(:,:,c) = uXY;
% % %     
% % %     
% % % end
toc
