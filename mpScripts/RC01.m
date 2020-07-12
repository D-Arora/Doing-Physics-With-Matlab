% RC01.m
% ap/mp/mscripts
% 17 nov 2013

% RC Circuits: Transient Responses
% Charging a capacitor
% SI units
% R resistancce   C capacitance   E emf (battery)
% V voltage       I currents      Q charge

% calls the function simpson1d.m

clear all
close all
clc

% INPUTS -----------------------------------------------------------------
R = 1e3;   C = 1e-6;   E = 100e-3;

num = 1000;   tmin = 0;   tmax = 5e-3;

% SETUP ------------------------------------------------------------------
Vc = zeros(num,1); VcN = zeros(num,1);   VR = zeros(num,1);
I = zeros(num,1);   Q = zeros(num,1);

t = linspace(tmin,tmax,num);              % time interval
tau = R * C;                              % time constant
dt = t(2) - t(1);

% Anaytical solution of Vc
Vc = E .* (1 - exp(-t ./ tau));

% Numerical solution of Vc
for cc = 1: num-1
VcN(cc+1) = VcN(cc) + (dt/tau) * (E - VcN(cc));
end

% V and I at t = tau
Vtau = E .* (1 - exp(-1));   Itau = (E/R)*exp(-1);

% Calculations
I = (E/R) .* exp(-t./tau);
Q = C .* Vc;
VR = R .* I;
PE = E .* I;
Pc = Vc .* I;
PR = VR .* I;

% Calaculations - energy
for cc = 1:num
WE(cc) = simpson1d(PE,tmin,t(cc));
Wc(cc) = simpson1d(Pc,tmin,t(cc));
WR(cc) = simpson1d(PR,tmin,t(cc));
end

% GRAPHICS ===============================================================

figure(99)     % circuit diagram ----------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.02 0.65 0.25 0.25]);

% lines
x = [2 10]; y = [1 1];
plot(x,y,'k')
axis([0 13 0 13]);
axis equal
hold on
x = [10 10]; y = [1 3];
plot(x,y,'k')
x = [10 10]; y = [5 7];
plot(x,y,'k')
x = [10 10]; y = [8 10];
plot(x,y,'k')
x = [7 10]; y = [10 10];
plot(x,y,'k')
x = [2 5]; y = [10 10];
plot(x,y,'k')
x = [2 2]; y = [5 10];
plot(x,y,'k')
x = [2 2]; y = [1 4.5];
plot(x,y,'k')
x = [5 7]; y = [10 11];
plot(x,y,'k','lineWidth',2)

% resistor / capacitor
x = [10-0.5 10+0.5]; y = [3 3];
plot(x,y,'k','lineWidth',2)
x = [10-0.5 10+0.5]; y = [5 5];
plot(x,y,'k','lineWidth',2)
x = [10-0.5 10-0.5]; y = [3 5];
plot(x,y,'k','lineWidth',2)
x = [10+0.5 10+0.5]; y = [3 5];
plot(x,y,'k','lineWidth',2)

x = [9 11]; y = [7 7];
plot(x,y,'k','lineWidth',3)
x = [9 11]; y = [8 8];
plot(x,y,'k','lineWidth',3)

x = [1 3]; y = [5 5];
plot(x,y,'r','lineWidth',4)
x = [1.5 2.5]; y = [4.5 4.5];
plot(x,y,'k','lineWidth',2)


title_main = 'E = 100 mV     R = 10^5 \Omega     C = 10^{-8} F'; 
title(title_main);
tt1 = 'E'; tt2 = 'R';   tt3 = 'C';
text(0,5,tt1);   text(11,4,tt2);   text(11.5,7.5,tt3);   

text(1.1,0.1,'Switched close at time t = 0 ');
axis off

figure(1)     % voltage --------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.3 0.65 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'V_R   V_C   (mV)';
title_main = 'Charging a Capacitor';
tt1 = 'time constant   \tau  =';
tt2 = num2str(1e3*tau,3);   tt3 = '  ms';
tt = [tt1 tt2 tt3];

x = t.*1e3;   y = Vc.*1e3;
plot(x,y,'linewidth',2);   %  Vc
text(1.5,50,tt);
xlabel(title_x); ylabel(title_y);
title(title_main);

hold on

x = t.*1e3; y = VR.*1e3;
plot(x,y,'r','linewidth',2);   % VR

x = t.*1e3;   y = VcN.*1e3; dnum = 50;
plot(x(1:dnum:num),y(1:dnum:num),'r+','linewidth',2)   % Vc numerical

legend('V_C','V_R','V_C(num)','location','east');

x = [tau.*1e3 tau.*1e3];   y = [0 Vtau.*1e3];
plot(x,y,'-')

x = [0 tau.*1e3];   y = [Vtau.*1e3 Vtau.*1e3];
plot(x,y,'-')

x = [0 tau.*1e3];   y = [E*exp(-1).*1e3 E*exp(-1).*1e3];
plot(x,y,'r-')

figure(2)   % current ----------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.58 0.65 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'current   I   ( \muA )';
title_main = 'Charging a Capacitor:  I = I_C = I_R';
tt1 = 'time constant   \tau  =  ';
tt2 = num2str(1e3*tau,3);   tt3 = '  ms';
tt = [tt1 tt2 tt3];

x = t.*1e3;   y = I.*1e6;
plot(x,y,'linewidth',2)

title(title_main);
xlabel(title_x); ylabel(title_y);
text(1.5,50,tt);

hold on
x = [tau.*1e3 tau.*1e3];   y = [0 Itau.*1e6];
plot(x,y,'-')

x = [0 tau.*1e3];   y = [Itau.*1e6 Itau.*1e6];
plot(x,y,'-')

figure(3)   % charge -----------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.02 0.25 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'Charge   Q   (C)';
title_main = 'Charging a Capacitor';
tt1 = 'time constant   \tau  =';
tt2 = num2str(1e3*tau,3);   tt3 = '  ms';
tt = [tt1 tt2 tt3];

x = t.*1e3;   y = Q;
plot(x,y,'r','linewidth',2)

title(title_main);
xlabel(title_x); ylabel(title_y);
text(1.5,0.5,tt);

hold on
x = t.*1e3;   y = -Q;
plot(x,y,'k','linewidth',2)

legend('Q: positive plate','Q: negative plate','location','east');

figure(4)     % power -----------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.3 0.25 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'power   P   (W)';
title_main = 'Charging a Capacitor';

x = t.*1e3;   y = PE;
plot(x,y,'k','linewidth',2)

title(title_main);
xlabel(title_x); ylabel(title_y);

hold on
x = t.*1e3;   y = Pc;
plot(x,y,'b','linewidth',2)

x = t.*1e3;   y = PR;
plot(x,y,'r','linewidth',2)

legend('P_E(battery)','P_C(capacitor)','P_R(resistor)','location','east');

figure(5)     % energy ---------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.58 0.25 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'energy  W   (J)';
title_main = 'Charging a Capacitor';

x = t.*1e3;   y = WE;
plot(x,y,'k','linewidth',2)

title(title_main);
xlabel(title_x); ylabel(title_y);

hold on
x = t.*1e3;   y = Wc;
plot(x,y,'b','linewidth',2)

x = t.*1e3;   y = WR;
plot(x,y,'r','linewidth',2)

grid on
legend('W_E(battery - supplied)','W_C(capacitor - stored)','W_R(resistor - dissipated)','location','northwest');

