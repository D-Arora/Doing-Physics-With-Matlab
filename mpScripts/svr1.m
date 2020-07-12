% SVR1.m




clear
close all
clc


b = 0.8;
k = 0.5;
N = 1001;
nT = 500;
tMax = 100;
t = linspace(0,tMax,nT);
dt = t(2) - t(1);


s = zeros(nT,1);
v = zeros(nT,1);
r = zeros(nT,1);

% First time step: initial conditions
  s(1) = 1000/N;
  v(1) = 1/N;
  r(1) = 0;


for n = 2:nT
  s(n) = s(n-1) - dt * b*s(n-1)*v(n-1);  
  v(n) = v(n-1) + dt * (b*s(n-1)*v(n-1) - k*v(n-1));
  r(n) = r(n-1) + dt * k*v(n-1);
%  s(n) = 1 - v(n) - r(n);       
end


% % DEATHS
% d = 4e-4;
% D = d*N*v;



NY = [14 28 50 66 156 190 156 108 68 77 33 65 24];
t_NY = 7:7:7*13;

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.3 0.5]);
  set(gcf,'color','w');

  subplot(3,1,1)
    plot(t,s*100,'k','linewidth',2)
    hold on  
    plot(t,v*100,'r','linewidth',2)
    plot(t,r*100,'b','linewidth',2)
    set(gca,'xtick',0:20:nT)
    grid on
    legend('s%','v%','r%','orientation','horizontal','location','northoutside')
    xlabel('days')
    ylabel('s, v, r')
    set(gca,'fontsize',12)
  
  subplot(3,1,2)
    plot(t,v*100,'r','linewidth',2)
    set(gca,'xtick',0:20:nT)
    grid on
    xlabel('days')
    ylabel('r%')
    set(gca,'fontsize',12)
  
%   subplot(3,1,3)
%     plot(t,D,'r','linewidth',2)
%    hold on
%    plot(t_NY,NY,'bo','linewidth',2)
%    set(gca,'xtick',0:20:nT)
%     
%     grid on
%     xlabel('days')
%     ylabel('DEATHS')
%     set(gca,'fontsize',12)  
    
   
