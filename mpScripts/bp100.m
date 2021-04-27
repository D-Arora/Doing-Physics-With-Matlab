% bp100.m

close all
clear
clc

Gavg(1)   = 101.29;
Gload(1)  = 15.75;
G1peaks   = 301.50;
Gt1peaks  = 1800;
Iavg(1)   = 8.73;
Iload(1)  = 0.88;
I1peaks   = 14.14;
Iitpeaks   = 1800;

Gavg(2)   = 100.86;
Gload(2)  = 15.33;
G2_peaks   = [195.76 191.83];
Gt2_peaks  = [1400 1800];
Iavg(2)   = 8.99;
Iload(2)  = 1.14;
I2_peaks   = [12.36 12.64];
I2_tpeaks   = [1400 1800];

Gavg(3)   = 100.07;
Gload(3)  = 15.54;
G3_peaks   = [195.76 194.77];
Gt3_peaks  = [1200 1800];
Iavg(3)   = 8.99;
Iload(3)  = 1.14;
I3_peaks   = [12.36 12.64];
I3_tpeaks   = [1200 1800];

Gavg(4)   = 100.98;
Gload(4)  = 15.45;
G4_peaks   = [160.00 158.34 159.63];
Gt4peaks  = [0800 1200 1800];
Iavg(4)   = 9.14;
Iload(4)  = 1.29;
I4_peaks   = [11.5 11.57 11.40];
I4_tpeaks   = [1200 1800];

Gavg(5)   = 102.17;
Gload(5)  = 16.63;
G5_peaks   = [160.00 158.34 159.63];
Gt5peaks  = [0800 1200 1800];
Iavg(5)   = 9.24;
Iload(5)  = 1.39;
I5_peaks   = [11.5 11.57 11.40];
I5_tpeaks   = [1200 1800];

% GRAPHICS  ==============================================================
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.05 0.45 0.45]);
   set(gcf,'color','w');
   FS = 14; LW = 10;
   box on
   hold on
   
     yP = 100*Gavg./Gavg(1);
   for c = 1:5
     plot([c c], [0 yP(c)],'k','linewidth',LW)
   end
   
     yP = 100*Gload./Gload(1);
   for c = 1:5
     plot([c+6 c+6], [0 yP(c)],'k','linewidth',LW)
   end
   
   yP = 100*Iavg./Iavg(1);
   for c = 1:5
     plot([c+12 c+12], [0 yP(c)],'b','linewidth',LW)
   end
   
     yP = 100*Iload./Iload(1);
   for c = 1:5
     plot([c+18 c+18], [0 yP(c)],'b','linewidth',LW)
   end
  
   text(1.5,-10,'G_{avg}','fontsize',FS)
   text(7.5,-10,'G_{load}','fontsize',FS)
   text(13.5,-10,'I_{avg}','fontsize',FS)
   text(19.5,-10,'I_{load}','fontsize',FS)
  
   set(gca,'xtick',[])
   set(gca, 'YGrid', 'on', 'XGrid', 'off')
   set(gca,'yTick',0:20:160)
   set(gca,'fontsize',FS)
   