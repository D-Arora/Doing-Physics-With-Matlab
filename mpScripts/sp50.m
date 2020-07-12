% sp50.m

%% How do the planets move
close all
clear all
clc

  a = [5.79e10 1.08e11 1.496e11 2.28e11 7.78e11 1.43e12 2.86e12 4.52e12 5.90e12];
  T = [7.60e6  1.94e7  3.156e7  5.94e7  3.74e8  9.35e8  2.64e9  5.22e9  7.82e9];
  
  G = 6.67e-11; MS = 1.95e30;
  K = 4*pi^2/(G*MS);
  aF = linspace(0,6e12,500);
  TF = sqrt(K .* aF.^3);
  
 figure(1)
 pos = [0.07 0.05 0.28 0.25];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
 
 xP = a; yP = T;
  plot(xP,yP,'bo');
  hold on
 
  xP = aF; yP = TF;
  plot(xP,yP,'r','linewidth',2);
  grid on
  xlabel('a [ m ]'); ylabel('T  [ s ]');
  title('T = 5.51x10^{-10} a^{3/2}');
  set(gca,'fontsize',12)
  
 figure(2)
 pos = [0.37 0.05 0.28 0.25];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
 
 xP = a.^1.5; yP = T;
  plot(xP,yP,'bo');
  hold on
  xP = aF.^1.5; yP = TF;
  plot(xP,yP,'r','linewidth',2);
  grid on
  xlabel('a^{3/2} [ m^{3/2} ]'); ylabel('T  [ s ]');
  title('T = 5.51x10^{-10} a^{3/2}');
  set(gca,'fontsize',12) 
  
 figure(3)
  pos = [0.67 0.05 0.28 0.25];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
 
  xP = a.^3; yP = T.^2;
  plot(xP,yP,'bo');
  hold on
  
  xP = aF.^3; yP = TF.^2;
  plot(xP,yP,'r','linewidth',2);
  grid on
  xlabel('a^{3} [ m^{3} ]'); ylabel('T^2  [ s^2 ]');
  title('T = 5.51x10^{-10} a^{3/2}');
  set(gca,'fontsize',12) 
  
figure(4)
  pos = [0.07 0.42 0.28 0.25];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
 
  xP = log10(a); yP = log10(T);
  plot(xP,yP,'bo');
  hold on
  
  xP = log10(aF); yP = log10(TF);
  plot(xP,yP,'r','linewidth',2);
  grid on
  xlabel('log_{10}a'); ylabel('log_{10}T');
  title('T = 5.51x10^{-10} a^{3/2}');
  set(gca,'fontsize',12)   
  