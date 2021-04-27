% cov192D.m


clear
close all
clc


% Population  N2
  N  = 5;  
  N2 = N^2;
  
% location  
%  x = randi([1 N],N,N);
%  y = randi([1 N],N,N);

x = 1:N;

[xx, yy] = meshgrid(x,x);


figure(1)
  Hplot = plot(xx,yy,'ok');
  set(Hplot,'markersize',6,'markerfacecolor','k');
  

  