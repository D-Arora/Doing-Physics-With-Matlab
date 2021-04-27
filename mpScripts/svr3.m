% svr3.m

close all
clear
clc

% Number of time steps [days]
  nT = 50;

% Initialize arrays
  N = zeros(nT,1);
  S = zeros(nT,1);
  D = zeros(nT,1);
  R = zeros(nT,1);
  I = zeros(nT,1);
  t = 1:nT;
  
% Initial values
  N(1) = 1001;
  S(1) = 1000;
  D(1) = 0;
  R(1) = 0;
  I(1) = 1;
  
% Parameters
  a = 2;
  b = 0.5;
  c = 5e-5;
  d = 1e-1;
  
% time evolution
for n = 2:nT
  S(n) = S(n-1) - a*S(n-1)*I(n-1)/N(n-1) + c*R(n-1);
  I(n) = I(n-1) + a*S(n-1)*I(n-1)/N(n-1) - b*I(n-1);
  R(n) = R(n-1) + b*I(n-1) - c*R(n-1);
  D(n) = D(n-1) + d*I(n);
  N(n) = S(n) + I(n) + R(n) - D(n);
  
end


figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.3 0.5]);
  set(gcf,'color','w');

  xP = t; yP = S;
    plot(xP,yP,'k','linewidth',2)
    hold on  
  yP = D;
    plot(xP,yP,'r','linewidth',2)
  yP = R;  
    plot(xP,yP,'b','linewidth',2)
  yP = I;  
    plot(xP,yP,'m','linewidth',2)
  yP = N;  
    plot(xP,yP,'c','linewidth',2)  
    
%     set(gca,'xtick',0:20:nT)
    grid on
   legend('S','D','r','I','N','orientation','horizontal','location','northoutside')
    xlabel('days')
    ylabel('s, v, r')
    set(gca,'fontsize',12)