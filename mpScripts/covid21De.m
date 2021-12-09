% covid21De.m


close all
clc
clear

%save covid21D.mat covid21D

load('covid21D') 

t = 1:length(covid21D);
t = t - 301;
D = covid21D;
D = smooth(t,D,7);

dDdt = gradient(D,1);
dDdt = smooth(t,dDdt,7);

% Predicted deaths
t1 = t(end); t2 = 500;
D1 = D(end);
dD = dDdt(end)*(t2-t1);
D2 = D1 + dD;

figure(9)  
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.5 0.2 0.35 0.35])
   set(gcf,'color','w');
   FS = 14;
   

subplot(2,1,1)
   plot(t,D/1000,'b','linewidth',2)

   hold on
   xP = [t1 t2]; yP = [D1 D2]./1000;
   plot(xP,yP,'r','linewidth',2)

   xtxt = 'days';
   ytxt = 'deaths / 1000';
   xlim([-300 500])
   xlabel(xtxt)
   ylabel(ytxt)

   txt = sprintf('D_{final} = %3.0f K  deaths \n',D2/1000);
   Htxt = text(240,50,txt);
   set(Htxt,'fontsize',12,'color','r')

   text(-295,190,'6 Mar 20')
   plot([-300 -300],[0 195],'k','linewidth',1)

   text(-85,190,'1 Jan 21')
   plot([0 0],[0 195],'k','linewidth',1)

   text(270,190,'1 Jan 22')
   plot([366 366],[0 195],'k','linewidth',1)

   text(395,190,'15 May 22')
   plot([500 500],[0 195],'k','linewidth',1)



   set(gca,'fontsize',FS)
   grid on; box on;

subplot(2,1,2)
   plot(t,dDdt,'b','linewidth',2)
  
   hold on
   xP = [t(end), 500];
   yP = [dDdt(end), dDdt(end)];
   plot(xP,yP,'r','linewidth',2)

   
   xtxt = 'days';
   ytxt = 'dD /dt';
   xlim([-300 500])
   xlabel(xtxt)
   ylabel(ytxt)
    
   txt = sprintf('(dD/dt)_{final} = %3.0f deaths/day \n',dDdt(end));
   Htxt = text(180,550,txt);
   set(Htxt,'fontsize',12,'color','r')

   text(-295,1090,'6 Mar 20')
   plot([-300 -300],[0 max(dDdt)],'k','linewidth',1)

   text(-85,1090,'1 Jan 21')
   plot([0 0],[0 max(dDdt)],'k','linewidth',1)

   text(270,1090,'1 Jan 22')
   plot([366 366],[0 max(dDdt)],'k','linewidth',1)

   text(395,1090,'15 May 22')
   plot([500 500],[0 max(dDdt)],'k','linewidth',1)

   set(gca,'fontsize',FS)
   grid on; box on;