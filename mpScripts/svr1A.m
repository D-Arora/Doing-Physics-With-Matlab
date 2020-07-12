% SVR1.m




clear 
close all
clc


b = 0.5;%0.5;
k = 1/3 ;%1/3;

%N = 7900000;
nT = 500;
tMax = 20;
t = linspace(0,tMax,nT);
dt = t(2) - t(1);


S = zeros(nT,1);
V = zeros(nT,1);
R = zeros(nT,1);


% First ime step: initial conditions
  S(1) = 10;
  V(1) = 1;
  R(1) = 0;
  N = S(1)+ V(1) + R(1);

for n = 2:nT
  S(n) = S(n-1) - dt * b*S(n-1)*V(n-1)/N;  
  V(n) = V(n-1) + dt * (b*S(n-1)*V(n-1)/N - k*V(n-1));
  R(n) = R(n-1) + dt * k*V(n-1);
%  s(n) = 1 - v(n) - r(n);       
end


% DEATHS
d = 4e-4;
D = d*N*V;



NY = [14 28 50 66 156 190 156 108 68 77 33 65 24];
t_NY = 7:7:7*13;

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.3 0.5]);
  set(gcf,'color','w');

  subplot(3,1,1)
    plot(t,S,'k','linewidth',2)
    hold on  
    plot(t,V,'r','linewidth',2)
    plot(t,R,'b','linewidth',2)
   % set(gca,'xtick',0:20:nT)
    grid on
    legend('s%','v%','r%','location','west')
    xlabel('days')
    ylabel('s, v, r')
    set(gca,'fontsize',12)
  
  subplot(3,1,2)
    plot(t,V,'r','linewidth',2)
  %  set(gca,'xtick',0:20:nT)
    grid on
    xlabel('days')
    ylabel('r%')
    set(gca,'fontsize',12)
  
%   subplot(3,1,3)
%     plot(t,D,'r','linewidth',2)
%    hold on
%    plot(t_NY,NY,'bo','linewidth',2)
%  %  set(gca,'xtick',0:20:nT)
%     
%     grid on
%     xlabel('days')
%     ylabel('DEATHS')
%     set(gca,'fontsize',12)  
    
   
