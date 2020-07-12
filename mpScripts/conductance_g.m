% conductance_g.m

% ap/mp/mscripts
% 15 11 2013


clear all
close all
clc

global Vr

T = 6.3;
Vr = -65;

Vt = Vr+56;

gK_max = 36;
gNa_max = 120;

tmin = 0;   tmax = 10;   num = 1000;

t = linspace(tmin, tmax,num);

[ An Am Ah ] = alpha(Vt,T);
[ Bn Bm Bh ] = beta(Vt,T);

[ Tn Tm Th ] = tau(Vt,T);

[ n_0 m_0 h_0 ] = N_0(T);

[ n_inf m_inf h_inf ] = N_inf(Vt,T);

nt = n_inf - (n_inf - n_0) .* exp(-t./Tn);
mt = m_inf - (m_inf - m_0) .* exp(-t./Tm);
ht = h_inf - (h_inf - h_0) .* exp(-t./Th);


gKt = gK_max .* nt.^4;
gNat = gNa_max .* mt.^3 .* ht;


figure(4)
tp1 = 'V  =  '; tp2 = num2str(Vt,3);  tp3 = '  mV';
tp = [tp1 tp2 tp3];
subplot(2,1,1)
plot(t,nt,'lineWidth',2)
hold on
plot(t,nt.^4,'r','lineWidth',2)
xlabel('t   (ms)   (mV)')
ylabel('nt ')
title(tp);

subplot(2,1,2)
plot(t,gKt,'lineWidth',2)
xlabel('t   (ms)   (mV)')
ylabel('gK ')

figure(5)
subplot(2,1,1)
plot(t,mt,'lineWidth',2)
hold on
plot(t,mt.^3,'c','lineWidth',2),
plot(t,ht,'r','lineWidth',2)
xlabel('t   (ms)   (mV)')
ylabel('mt ')
title(tp);

subplot(2,1,2)
plot(t,gNat,'lineWidth',2)
xlabel('t   (ms)   (mV)')
ylabel('gNa ')


