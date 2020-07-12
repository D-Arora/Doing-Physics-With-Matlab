% conductance_g.m

% ap/mp/mscripts
% 23 11 2013

% Voltage clamp

clear all
close all
clc

global Vr

num = 1000;   num1 = round(0.1 * num);   num2 = round(0.6 * num);
dV = 65;

T = 6.3;
Vr = -65;

ENa = 50;  EK = -77;   EL = -75.6;      % reversal potentials [mV]

%Vt = Vr+56;

Vt = Vr .* ones(num,1);
JNa = zeros(num,1);   JK = zeros(num,1);
gKt = zeros(num,1); gNat = zeros(num,1);
nt = zeros(num,1); mt = zeros(num,1); ht = zeros(num,1);

Vt(num1:num2) = Vr + dV;


gK_max   = 36;         % max conductance K         [mmohm^-1.cm-2]
gNa_max  = 120;        % max conductance Na        [mmohm^-1.cm-1]
gL_max   = 0.3;        % max conductance leakage   [mmohm^-1.cm-2]

tmin = 0;   tmax = 20;   

t = linspace(tmin, tmax,num);
dt = t(2) - t(1);

[ An Am Ah ] = alpha(Vt,T);
[ Bn Bm Bh ] = beta(Vt,T);

[ Tn Tm Th ] = tau(Vt,T);

[ n_0 m_0 h_0 ] = N_0(T);

%n_0 = 0; m_0 = 0; h_0 = 0;

[ n_inf m_inf h_inf ] = N_inf(Vt,T);


nt(1) = An(1) / (An(1) + Bn(1));
mt(1) = Am(1) / (Am(1) + Bm(1));
ht(1) = Ah(1) / (Ah(1) + Bh(1));

gKt(1)  = gK_max .* nt(1).^4;
gNat(1) = gNa_max .* mt(1).^3 .* ht(1);

JK(1)  = gKt(1)  * (Vt(1) - EK) / 1000;   % /1000  mV --> V   [mohm.cm^-2
JNa(1) = gNat(1) * (Vt(1) - ENa)/ 1000;

for c = 1 : num-1
%   nt(c) = n_inf(c) - (n_inf(c) - n_0) .* exp(-t(c)./Tn(c));
%   mt(c) = m_inf(c) - (m_inf(c) - m_0) .* exp(-t(c)./Tm(c));
%   ht(c) = h_inf(c) - (h_inf(c) - h_0) .* exp(-t(c)./Th(c));
  
 nt(c+1) = nt(c) + dt * (An(c) *(1-nt(c)) - Bn(c) * nt(c)); 
 mt(c+1) = mt(c) + dt * (Am(c) *(1-mt(c)) - Bm(c) * mt(c)); 
 ht(c+1) = ht(c) + dt * (Ah(c) *(1-ht(c)) - Bh(c) * ht(c)); 
end

%   gKt(c+1) = gK_max .* nt(c+1).^4;
%   gNat(c+1) = gNa_max .* mt(c+1).^3 .* ht(c+1);
% 
%   IK(c+1)  = gKt(c+1)  * (Vt(c+1) - VK);
%   INa(c+1) = gNat(c+1) * (Vt(c+1) - VNa);
  
%end

gKt = gK_max .* nt.^4;
gNat = gNa_max .* mt.^3 .* ht;

JK  = gKt  .* (Vt - EK)/1000;
JNa = gNat .* (Vt - ENa)/1000;
JL = gL_max .* (Vt - EL)/1000; 
Jm = JK + JNa + JL;

%GRAPHICS ================================================================
figure(1) 
set(gcf,'units','normalized');
set(gcf,'position',[0.05 0.05 0.85 0.85]);
fs = 15;

subplot(3,2,1)
set(gca,'fontsize',fs);
xTitle = 'time   t   (ms)';
yTitle = 'Voltage clamp   V   (mV)';
gTitle1 = 'dV  =   ';  gTitle2 = num2str(dV,3); gTitle3 = ' mV';
gTitle = [gTitle1 gTitle2 gTitle3];

plot(t,Vt,'lineWidth',2)
xlabel(xTitle)
ylabel(yTitle)
title(gTitle)

subplot(3,2,2)
yTitle = 'J   (mA.cm^{-2})';
plot(t,JNa,'lineWidth',2)
set(gca,'fontsize',fs);
hold on
plot(t,JK,'r','lineWidth',2);
plot(t,JL,'m','lineWidth',2);
plot(t,Jm,'k','lineWidth',2);
legend('J_{Na}', 'J_K', 'J_L', 'J_m');
xlabel(xTitle)
ylabel(yTitle)

subplot(3,2,3)
plot(t,mt,'lineWidth',2)
yTitle = 'gate variables';
set(gca,'fontsize',fs);
hold on
plot(t,mt.^3,'r','lineWidth',2)
plot(t,ht,'k','lineWidth',2)
xlabel(xTitle)
ylabel(yTitle)
legend('m', 'm^3', 'h');

subplot(3,2,4)
plot(t,nt,'lineWidth',2)
set(gca,'fontsize',fs);
hold on
plot(t,nt.^4,'r','lineWidth',2)
xlabel(xTitle)
ylabel(yTitle)
legend('n', 'n^3');

subplot(3,2,5)
yTitle = ' g_K   (mohm^{-1}.cm^{-2})';
plot(t,gKt,'r','lineWidth',2)
set(gca,'fontsize',fs);
xlabel(xTitle)
ylabel(yTitle)

subplot(3,2,6)
yTitle = ' g_{Na}   (mohm^{-1}.cm^{-2})';
plot(t,gNat,'lineWidth',2)
set(gca,'fontsize',fs);
xlabel(xTitle)
ylabel(yTitle)





