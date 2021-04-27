% Laplace01.m

close all
clc
clear

global w

N = 5001;
m = 1;

A = 1;
w = 3;
k = 3;
b = 1;

tS = linspace(0,8*pi,N);

x = zeros(N,1);
y = zeros(N,1);


y = A.* cos(w.*tS);

yDot = -(A*w) .* sin(w.*tS);

K(1) = -b/m;
K(2) = -k/m;
K(3) = k*A/m;
K(4) = -b*A*w/m;

u0 = [0; A*w];

u0 = [0, 0];
tSpan = linspace(0,8*pi, N); 

RelTol = 1e-6;
options = odeset('RelTol',RelTol);


[t, SOL] = ode45(@(t,u) FNode(t,u,K), tSpan, u0, options);

x = SOL(:,1);

phi = -1.2*pi;
xP = (1/sqrt(85)).*cos(w*t-phi);

% GRAPHICS ===============================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.05 0.45 0.45]);
   set(gcf,'color','w');
   FS = 14; LW = 2;
   
   box on

   plot(t./pi,x,'r','linewidth',LW);

   hold on
   plot(tS./pi,y,'b','linewidth',LW)

   plot(t./pi,xP,'k','linewidth',LW);
   
   
   
   ylim([-2.5 2.5])
   grid on
   
   xlabel('t/\pi')
   ylabel('x & y')
   legend('Response','Input')
   set(gca,'fontsize',FS)
   
   
   


% FUNCTIONS  =============================================================   
   
function du = FNode(t,u,K)

global w

z = u(1);

zDot = u(2);

du = zeros(2,1);

% First derivative dz/dt

du(1) = zDot;

% ODE: Second derivative d2y/dt2

%du(2) = K(1)*zDot + K(2)*z + K(3)*cos(w*t) + K(4)*sin(w*t);

du(2) = -2*zDot - 2*z + cos(w*t);

end