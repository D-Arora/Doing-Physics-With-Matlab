% mw_doppler.m

% Doppler Effect

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170515


close all;
clear all;
clc;

v   = 340;        % velocity of sound
freq = 1000;      % source frequency


figure(2)
pos = [0.1 0.1 0.3 0.4];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos)
% Source moving towards a stationary observer
   v = 340;
   nF = 5000;
   vs = linspace(0,v-10,nF);
   vo = 0;
   fs = freq;
   fo = fs .* v ./ (v - vs);

   plot(vs,fo,'b','lineWidth',2);
   hold on
   plot([v v],[-1000 max(fo)],'g','lineWidth',2);

% Source moving away a stationary observer
   vs = linspace(0,600,nF);
   fo = fs .* v ./ (v + vs);
  plot(vs,fo,'k','lineWidth',2);
  fo = fs .* v ./ (v + vs);

% Observer moving towards stationary source 
  vo = linspace(0,1000,nF);
  fo = fs .* (v+vo) ./ v;
  plot(vo,fo,'r','lineWidth',2);

% Observer away stationary source 
  vo = linspace(0,600,nF);
  fo = fs .* (v-vo) ./ v;
  fo(fo<0) = 0;
  plot(vo,fo,'m','lineWidth',2); 

  
hText=text(350,4200,'v(sound)','fontsize',14);
set(hText,'color','g')

hText = text(32,3800,'souce moving','fontsize',14);
set(hText,'color','b')
hText = text(32,3400,'towards observer','fontsize',14);
set(hText,'color','b');
  
hText = text(360,1200,'souce moving','fontsize',14);
set(hText,'color','k');
hText = text(360,800,'away from observer','fontsize',14);
set(hText,'color','k');

hText = text(360,3200,'observer moving','fontsize',14);
set(hText,'color','r');
hText = text(360,2800,'towards source','fontsize',14);
set(hText,'color','r');


hText = text(40,0,'observer moving','fontsize',14);
set(hText,'color','m');
hText = text(40,-400,'away from source','fontsize',14);
set(hText,'color','m');

xlabel('speed of source or observer [ m/s ]');
ylabel('observer freq  f_o  [ Hz ]');
axis([0 600 -1000 5000])
grid on

set(gca,'fontsize',14);
