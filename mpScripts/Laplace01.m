% Laplace01.m

close all
clc
clear

global w


% System: INPUTS  ========================================================
  N = 5001;
  tMax = 8*pi;
  
  w = 1;

  m = 1;
  b = 1;
  k = 2;

% Setup  ================================================================

  tS = linspace(0,8*pi,N);

  x = zeros(N,1);           % output

  y = zeros(N,1);           % input

  %yS = -(b*w).*sin(w*tS);
   yS = cos(w*tS);

K(1) = -b/m;
K(2) = -k/m;
K(3) = -b*w/m;

% Solve ODE  ============================================================

u0 = [0, 0];
tSpan = linspace(0,8*pi, N); 

RelTol = 1e-6;
options = odeset('RelTol',RelTol);


[t, SOL] = ode45(@(t,u) FNode(t,u,K), tSpan, u0, options);

x = SOL(:,1);


% GRAPHICS ===============================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.05 0.45 0.45]);
   set(gcf,'color','w');
   FS = 14; LW = 2;
   
   box on

   plot(t./pi,x,'r','linewidth',LW);

   hold on
   plot(tS./pi,yS,'b','linewidth',LW)

   ylim([-2.5 2.5])
   grid on
   
   xlabel('t/\pi')
   ylabel('x & y')
   legend('Response','Input')
   set(gca,'fontsize',FS)
   
   title(max(x))
   


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
du(2) = K(1)*zDot + K(2)*z + K(3)*sin(w*t) ;
%du(2) = -0.3*zDot - 2*z; % + cos(1.4*t);
end