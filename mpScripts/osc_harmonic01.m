% osc_harmonic01.m

% Spring / Mass System
% Driven harmonic oscillation with damping
% Finite difference Method for Calculations
% 180301

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm


close all
clear 
clc

% =======================================================================
% INPUTS   SI values unless stated otherwise  default values  [  ]
% =======================================================================
% mass  [0.506606]
   m = 0.506606; 
% spring constant   [20]
   k = 20; 
 
clc

disp('   ')
   f0 = sqrt(k/m)/(2*pi);       % resonance (natural) frequency
   fprintf('Resonance frequency  f0 = %3.5f Hz  \n',f0)
disp('   ')
% damping constant   [ 0 to 20 ?]
   b = input('\n   Damping constant [0 to 20]  b =  '); 
%  flagF
    flagF = input('\n   Driving force: 1 (Sinusoidal) or 2 (Impulsive) or 3 (free oscillation)  flagF =  ');
% driving frequency
   if flagF == 1
      fD = input('\n   Driving frequency [0.1 to 4]  fD = ');
   end
% max simulation time interval  [10]   
   tMax = input('\n   max simulation time interval [10]  tMax = ');                   
   
   
% =======================================================================
% SETUP
% =======================================================================
   T0 = 1/f0;                   % natural period of oscillation
   w0 = 2*pi*f0;                % natural angular frequency
   nmax = 8001;                 % max number of time steps
   tMin = 0;                    % tmin
  
   dt = (tMax-tMin)/(nmax-1);   % time step
   t = tMin : dt : tMax;        % simulation time
      
   x = zeros(nmax,1);    % displacement
   v = zeros(nmax,1);    % velocity
   a = zeros(nmax,1);    % acceleration
   %K = zeros(nmax,1);   % kinetic energy
   %U = zeros(nmax,1);   % elastic potential energy
   %E = zeros(nmax,1);   % total energy
   FD = zeros(nmax,1);   % force acting on mass


%  driving force

if flagF == 2
    FD(500:3000) = 10;     %   impulsive driving force
end
    
if flagF == 1
   wD = 2*pi*fD;
   FD = 10*sin(wD*t);     %   sinusoidal driving force
   x(2) = x(1)*sin(2*pi*dt/T0);
end

if flagF == 3
    x(1) = 1; 
    x(2) = x(1)*cos(w0*dt);
end

% Coefficients
c0 = 1 + (b/m)*(dt/2);
c1 = (2-(k/m)*dt*dt)/c0;
c11 = (2-(k/m)*dt*dt);
c2 =((b/m)*(dt/2)-1)/c0;
c3 = (dt*dt/m)/c0;

%  Finite difference Calculations
%  Position
for c = 3 : nmax
   x(c) = c1*x(c-1) + c2*x(c-2) + c3*FD(c-1);
end

% Velocity
v(1) = (x(2)-x(1))/dt;
v(nmax) = (x(nmax)- x(nmax-1))/dt;
for n = 2 : nmax-1
   v(n) = (x(n+1)- x(n-1))/(2*dt);
end

% Acceleration
a(1) = (v(2)-v(1))/dt;
a(nmax) = (v(nmax)- v(nmax-1))/dt;
for n = 2 : nmax-1
   a(n) = (v(n+1)- v(n-1))/(2*dt);
end

% Energies
K = (0.5*m).*v.^2;
U = (0.5*k).*x.^2;
E = K + U;

%s = 'damping coeff. b = '& numStr(b);
s = sprintf('b = %.1f ',b);

% ======================================================================
%  GRAPHICS
% ======================================================================
figure(1)
   pos = [0.07 0.05 0.32 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
plot(t,x,'LineWidth',2); 
title(s)
ylabel('position  x  (m)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
set(gca,'fontsize',14)

figure(2)
   pos = [0.37 0.05 0.32 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
plot(t,v,'LineWidth',3);
title(s)
ylabel('velocity  v  (m/s)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
set(gca,'fontsize',14)

figure(3)
   pos = [0.67 0.05 0.32 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
plot(t,a,'LineWidth',3);
title(s)
ylabel('acceleration  a  (m/s�)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
set(gca,'fontsize',14)

figure(4)
   pos = [0.07 0.45 0.32 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
plot(x,v,'LineWidth',3);
title(s)
ylabel('velocity  v  (m/s)','FontSize',18)
xlabel('position  x  (m)','FontSize',18)
grid on
set(gca,'fontsize',14)

figure(5)
   pos = [0.47 0.05 0.32 0.62];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
subplot(3,1,1), plot(t,x,'LineWidth',3);
title(s)
ylabel('position  x  (m)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
set(gca,'fontsize',14)

subplot(3,1,2), plot(t,v,'LineWidth',3);
%title(s)
ylabel('velocity  v  (m/s)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
set(gca,'fontsize',14)

subplot(3,1,3), plot(t,a,'LineWidth',3);
%title(s)
ylabel('accleration  a  (m/s�)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
set(gca,'fontsize',14)

figure(6)
   FS = 14;
   pos = [0.67 0.45 0.32 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
plot(t,K,'r','LineWidth',2)
hold on
plot(t,U,'b','LineWidth',2);
plot(t,E,'k','LineWidth',2);

title(s)
ylabel('energy K U E (J)','FontSize',18)
xlabel('time  t  (s)','FontSize',18)
grid on
hold off
legend('K','U','E','location','northwest')
set(gca,'fontsize',14)




