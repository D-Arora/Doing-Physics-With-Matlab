% chaos23.m



clear
%close all
clc

% INPUTS ==========================================================

col = [1 0 0];
%col = [1 0 0];
x0 = 0.1 ;
y0 = 0.1 ;
z0 = 0.1;
t1 = 0;
t2 = 50;

sigma = 10;
beta = 8/3;

%rho = 14;  %28
 N = 99;
 rho = linspace(0,50,N);


% SETUP =============================================================
tSpan = [t1 t2];
u0 = [x0;y0;z0];
xE = zeros(N,1); yE = xE; zE = xE;

for c = 1:N
K = [sigma,beta,rho(c)];

[t, SOL] = ode45(@(t,u) FNode(t,u,K), tSpan, u0);

x = SOL(:,1); y = SOL(:,2); z = SOL(:,3);
xE(c) = x(end); yE(c) = y(end); zE(c) = z(end);
end
rE = sqrt(xE.^2 + yE.^2 + zE.^2);

% CRITICAL POINTS  ================================================
% eta = sqrt(beta*(rho-1));
% cP1(3) = rho-1;
% cP1(2) = eta;
% cP1(1) = eta;
% cP2(3) = rho-1;
% cP2(2) = -eta;
% cP2(1) = -eta;
 
% GRAPHICS  ========================================================
  FS = 14;
figure(1)
  pos = [0.05 0.05 0.35 0.75];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
  xP = rho;
  hold on
subplot(4,1,1)
  hold on
  yP = xE;
  Hplot = plot(xP,yP,'o','color',col);
  set(Hplot,'MarkerSize',3,'markerfacecolor',col)
  ylabel('xE')
  grid on; box on
  set(gca,'FontSize',FS)
  hold on; box on
subplot(4,1,2)
  hold on
  yP = yE;
  Hplot = plot(xP,yP,'o','color',col);
  set(Hplot,'MarkerSize',3,'markerfacecolor',col)
  grid on; box on
  ylabel('yE')
  set(gca,'FontSize',FS)
  hold on
subplot(4,1,3)
  hold on
  yP = zE;
  Hplot = plot(xP,yP,'o','color',col);
  set(Hplot,'MarkerSize',3,'markerfacecolor',col)
  ylabel('zE')
  grid on;box on
  set(gca,'FontSize',FS)
  hold on
subplot(4,1,4)
  hold on
  yP = rE;
  Hplot = plot(xP,yP,'o','color',col);
  set(Hplot,'MarkerSize',3,'markerfacecolor',col)
  xlabel('\rho')
  ylabel('rE')
  grid on; box on
  set(gca,'FontSize',FS)
  hold on

%%

% FUNCTIONS ==========================================================

function du = FNode(t,u,K)
   sigma = K(1); beta = K(2); rho = K(3);
   x = u(1); y = u(2); z = u(3);
   du(1) = sigma*(y - x);
   du(2) = x*(rho - z) - y;
   du(3) = x*y - beta*z;
   du = du';
end






