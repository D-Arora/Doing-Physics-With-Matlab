% Laplace10.m

close all
clc
clear

w = 2;

phi = pi;

tMax = 6;


N = 501;

t = linspace(0,pi*tMax,N);

y = cos(w*t);

z = cos(w*t - phi);
%z = 0.5.*cos(w*(t - 1.08));

figure(1)

subplot(2,1,1)
  plot(t./pi,y,'b','linewidth',2)
  hold on
  plot(t./pi,z,'r','linewidth',2)
  grid on

subplot(2,1,2)
  c = linspace(-pi,pi,501);
  plot(cos(c),sin(c))
  hold on
  plot([0 0],[-1 1],'k')
  plot([-1 1],[0 0],'k')
  plot([ 0 sin(w*t(1))],[0 cos(w*t(1))],'b','linewidth',2)
  plot([ 0 -sin(w*t(1)-phi)],[0 cos(w*t(1)-phi)],'r','linewidth',2)
  axis square

