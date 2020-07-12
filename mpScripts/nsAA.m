% ns20_001.m

clear
close all
clc


global   a b c d I dt


% CONSTANTS
 %a = 0.02; b = 0.20; c = -65; d = 8; I = 10;
 a = -0.02; b = -1; c = -60; d = 8; I = 80;

 % Time
  dt = 0.005;
  NT = 2000;
  t = 0:dt:dt*(NT-1);

v = zeros(NT,1);   u = zeros(NT,1);  
v(1) = c;
u(1) = v(1)*b;

for k = 1:NT-1
  v(k+1) = V(v(k), u(k), dt, I);
  u(k+1) = U(v(k), u(k), dt, a, b);  
  
  if v(k+1) > 30
     v(k+1) = c;
     u(k+1) = u(k+1) + d;
  end
end

% NULLCLINES
  vNull = 0.04.*v.^2 + 5.*v + 140 + I; 
  uNull = b.*v;

% GRAPHICS  ===========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.25 0.5]);
  set(gcf,'color','w');

  xP = t;
subplot(2,1,1)
  yP = v;
  plot(xP,yP,'b','linewidth',2)
  xlabel('t  [ms]')
  ylabel('v_{membrane}  [mV]')
  grid on
 % xlim([0 12])
 txt = sprintf('I_{ext} = %3.0f',I);
  title(txt)
  set(gca,'fontsize',12)
 
subplot(2,1,2)
  yP = u;
  plot(xP,yP,'b','linewidth',2)
  xlabel('t  [ms]')
  ylabel('u_{recovery}  [mV]')
  grid on
%  xlim([0 12])
%  ylim([-80 80])
  set(gca,'ytick',-80:20:80)
  set(gca,'fontsize',12)

figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.32 0.05 0.25 0.25]);
  set(gcf,'color','w');
  
  xP = v; yP = u;
  plot(xP,yP,'b','linewidth',2)
  hold on
  Hplot = plot(v(1),u(1),'og');
  set(Hplot,'markerfacecolor','g','markersize',8)
  Hplot = plot(v(end),u(end),'or');
  set(Hplot,'markerfacecolor','r','markersize',8)
  
  xP = v; yP = vNull;
  plot(xP,yP,'r','linewidth',1)
  xP = v; yP = uNull;
  plot(xP,yP,'m','linewidth',1)
  
  xlabel('v_{membrane}  [mV]')
  ylabel('u_{recovery}  [mV]')
  grid on
  ylim([-65 2.1*max(u)])
  set(gca,'ytick',-80:20:100)
  txt = sprintf('I_{ext} = %3.0f',I);
  title(txt)
  set(gca,'fontsize',12)

  % FUNCTIONS  ========================================================
  function v = V(v, u, dt, I)
   v = v + dt*(0.04*v*v + 5*v + 140 - u + I);
  end

function u = U(v, u,dt, a, b)
   u = u + dt*a*(b*v - u);
end