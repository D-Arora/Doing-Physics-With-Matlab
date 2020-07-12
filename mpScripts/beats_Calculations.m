% beats_Calculation.m

% Plots of beat patterns

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170511

clear all
close all
clc


f2 = 1000;

f1 = 1050;

fBeat = abs(f1-f2);
fAvg = (f1+f2)/2;
tMax = 0.1;

t = linspace(0,tMax, 5000);


y1 = sin(2*pi*f1*t);
y2 = sin(2*pi*f2*t);
y = y1 + y2;

figure(1)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.1 0.1 0.5 0.7]);
   
   
   subplot(3,1,1);
   plot(t,y1,'b','linewidth',2);
   hold on
   plot(t,y2,'r','linewidth',2);
   axis([0 0.02 -1.1 1.1]);
   grid on
   set(gca,'fontsize',14);
   tm1 = 'f_1  =  ';
   tm2 = num2str(f1,'%3.0f\n');
   tm3 = '  Hz    ';
   tm4 = 'f_2  =  ';
   tm5 = num2str(f2,'%3.0f\n');
   tm6 = '  Hz';
   tm = [tm1 tm2 tm3 tm4 tm5 tm6];
   title(tm);
   
   
   subplot(3,1,2);
   plot(t,y,'k','linewidth',2);
   axis([0 0.02 -2.2 2.2]);
   grid on
   set(gca,'fontsize',14);
   tm1 = 'f_{avg}  =  ';
   tm2 = num2str(fAvg,'%3.0f\n');
   tm3 = '  Hz';
   tm = ([tm1 tm2 tm3]);
   title(tm);
   
   subplot(3,1,3);
   plot(t,y,'k','linewidth',2);
   axis([0 0.1 -2.2 2.2]);
   grid on
   tm1 = 'f_{beat}  =  ';
   tm2 = num2str(fBeat,'%3.0f\n');
   tm3 = '  Hz';
   tm = ([tm1 tm2 tm3]);
   title(tm);
   
   xlabel('time  t  [ s ]','fontsize',14);
   set(gca,'fontsize',14);
   
   
% figure(2)
%    plot(t,y,'b','linewidth',2);
%    axis([0 0.1 -2.2 2.2]);
%    axis off
   
   