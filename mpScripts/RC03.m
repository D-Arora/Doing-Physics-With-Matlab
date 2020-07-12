% RC03.m
% ap/mp/mscripts
% 17 nov 2013

% RC Circuits: STEP INPUT CURRENT
% Discharging a capacitor
% SI units
% R resistancce   C capacitance   
% V voltage       I currents      Q charge

clear all
close all
clc

% INPUTS -----------------------------------------------------------------
R = 1e3;   C = 1e-6;   I0 = 1*9e-6;

num = 1000; num1 = 100;   num2 = num/2;
tmin = 0;   tmax = 10e-3;

% SETUP ------------------------------------------------------------------
V = zeros(num,1); 
I = zeros(num,1);   Ic = zeros(num,1);     IR = zeros(num,1);

t = linspace(tmin,tmax,num);              % time interval
tau = R * C;                              % time constant
dt = t(2) - t(1);

% Input current - shape of pulse 
I(num1:num2) = I0;    % single pulses
%I(100:200) = I0;   I(400:500) = I0;   I(700:800) = I0;  % pulses
%I(50:100) = I0;   I(150:200) = I0;   I(250:300) = I0;  % pulses 
%I(350:400) = I0;   I(450:500) = I0;   I(550:600) = I0;
%I(650:700) = I0;   I(750:800) = I0;   I(850:900) = I0; I(950:1000) = I0

% Calculations - finite difference method
for cc = 1:num-1
V(cc+1) = V(cc) + (dt/C)*I(cc) - (dt/tau)*V(cc);
end

IR = V ./R;   Ic = I - IR;
Q = C .* V;


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
x = [10 10]; y = [6 10];
plot(x,y,'k')
x = [5 10]; y = [10 10];
plot(x,y,'k')
x = [2 5]; y = [10 10];
plot(x,y,'k')
x = [2 2]; y = [5 10];
plot(x,y,'k')
x = [2 2]; y = [1 4.5];
plot(x,y,'k')
x = [6 6]; y = [10 13];
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

x = [1 3]; y = [5 5];
plot(x,y,'r','lineWidth',3)
x = [1 3]; y = [4.5 4.5];
plot(x,y,'k','lineWidth',3)


title_main = 'V_0 = 100 mV     R = 10^5 \Omega     C = 10^{-8} F'; 
title(title_main);
tt2 = 'C';   tt3 = 'R';   tt1 = '\downarrow  input current I';
text(0,5,tt2);   text(11.5,4,tt3);   

axis off

figure(1)     % currents --------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.3 0.65 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'currents  (uV)';
title_main = 'Response: STEP or PULSE INPUT CURRENT';
tt1 = 'time constant   \tau  =';
tt2 = num2str(1e3*tau,3);   tt3 = '  ms';
tt = [tt1 tt2 tt3];

x = t.*1e3;   y = I.*1e6;
plot(x,y,'k','linewidth',2);   %  I
xlabel(title_x); ylabel(title_y);
title(title_main);
% 
hold on
x = t.*1e3;   y = IR.*1e6;
plot(x,y,'r','linewidth',2);   %  IR

x = t.*1e3;   y = Ic.*1e6;
plot(x,y,'b','linewidth',2);   %  IR

legend('input I','I_R','I_C','location','southwest');
grid on


figure(2)   % voltage ----------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.58 0.65 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'potential  V  ( mV )';

x = t.*1e3;   y = V*1e3;
plot(x,y,'linewidth',2)

title(title_main);
xlabel(title_x); ylabel(title_y);
grid on

figure(3)   % charge -----------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.02 0.25 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'Charge   Q   (C)';

x = t.*1e3;   y = Q;
plot(x,y,'r','linewidth',2)

%title(title_main);
xlabel(title_x); ylabel(title_y);

hold on
x = t.*1e3;   y = -Q;
plot(x,y,'k','linewidth',2)

legend('+ plate','- plate','location','southeast');

figure(4)   % input current  ----------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.3 0.25 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'current I   ( \muA )';

x = t.*1e3;   y = I.*1e6;
plot(x,y,'r','linewidth',2)

%title(title_main);
xlabel(title_x); ylabel(title_y);

