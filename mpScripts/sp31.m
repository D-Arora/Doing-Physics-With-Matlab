% sp31.m




%% P31 005

clear all
close all
clc

m = 9e-3;
Ls = 1.10;

mu = m/Ls;

f = [125 256 320 384 512];
L = [0.79 0.39 0.312 0.263 0.193];

% string tension
T = 4*f(2)^2*L(2)^2*mu;

Lmax = 1.00; Lmin = 0.1;
LT = linspace(Lmin,Lmax,500);

fT = 1./(2.*LT) .* sqrt(T/mu); 
 
figure(1)
  pos = [0.07 0.05 0.28 0.42];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
   

subplot(2,1,1);
  LW = 2; col = 'b';
  xP = L; yP = f;
  plot(xP,yP,'bo','linewidth',1)
  hold on
  xP = LT; yP = fT;
  plot(xP,yP,'r','linewidth',2)
  grid on
  set(gca,'fontsize',12);
  xlabel('string length L  [ m ]');
  ylabel('frequency f  [ Hz ]');
  

subplot(2,1,2);
  LW = 2; col = 'b';
  xP = 1./L; yP = f;
  plot(xP,yP,'bo','linewidth',1)
  hold on
  %xP = 1./LT; yP = fT;
  xP = [0 10]; yP = [0 1000];
  plot(xP,yP,'r','linewidth',2)
  grid on
  set(gca,'fontsize',12);
  xlabel('1/L  [ m^{-1} ]');
  ylabel('f  [ Hz ]');
  htext = text(0.5,700,'slope = rise / run = 1000/10 = 100');
  set(htext,'fontsize',12');
  
  