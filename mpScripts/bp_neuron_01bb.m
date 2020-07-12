% bp_neuron_01bb.m

% Frequency of repetitive firing of neuron as a fuunction of step height
% for a constant current input.

clear 
close all
clc

% Data

J0 = [0.008 0.01 0.02 0.03 0.04 0.05 0.06 0.08 0.10]

f = [159 190 256 297 330 356 378 415 443]


figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.45 0.35 0.30]);
  set(gcf,'color','w');
  title_x = 'current density  J_0  [ mA.cm^{-2} ]';
  title_y = 'spike frequency  f  [ Hz ]';
  title_m = 'Constant Current Injection: Step Function';

  x = J0;   y = f;
  plot(x,y,'-o','linewidth',2);   % membrane voltage
  xlabel(title_x); ylabel(title_y);
  title(title_m)
  axis([0 0.1 0 450])
  set(gca,'fontsize',12)
  grid on
  box on
