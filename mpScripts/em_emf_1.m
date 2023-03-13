% em_emf_1.m
% Animation of a magnet moving through a coil

% Ian Cooper
% School of Physics University of Sydney
% email: cooper@physics.usyd.edu.au

clear all; close all; clc

   xmax = 40.1; xmin = -30;   % limits for X axis
   N = 200;             % steps
   Ni = round(N * 12 / 50);    % location of current peak

figure(1)  % **********************************************************
   set(gcf,'color',[1 1 1]);
   set(gcf,'units','normalized'); 
   set(gcf,'position',[0.1 0.1 0.6 0.6]);
   set(gca,'Xlim',[xmin xmax]);
   set(gca,'xcolor',[1 1 1]);
   set(gca,'ycolor',[1 1 1]);


% plot coil -----------------------------------------------------------
   Lw = 2;
   xc1 = -8; xc2 = 8;
   plot([xc1,xc2],[5 5],'k','linewidth',Lw);
   hold on
   plot([xc2,xc2],[1 10],'k','linewidth',Lw);
   plot([xc1,xc2],[10 10],'k','linewidth',Lw);
   plot([xc1,xc1],[1 10],'k','linewidth',Lw);
   plot([xc1,xc1+6],[1 1],'k','linewidth',Lw);
   plot([xc2,xc2-6],[1 1],'k','linewidth',Lw);
   x1 = linspace(xc1,xc2,20); x2 = x1 + 0.5;
   rectangle('Position',[-2,-1,4,4],'Curvature',[1 1],'linewidth',Lw)
   ht = text(-1,1,'G'); set(ht,'FontSize',18,'Color',[0 1 0],'Fontweight','bold');
   y1 = 5*ones(20,1); y2 = y1 + 5;
     for c = 1 : 19
        plot([x1(c),x2(c)],[y1(c),y2(2)],'k','linewidth',Lw);
     end
   plot(x1,y1)

% plot current graph axes -------------------------------------------- 
   xd = 20; Lw1 = 1;
   plot([-xd,xd],[-8,-8],'k','linewidth',Lw1);
   plot([-xd,-xd],[0,-16],'k','linewidth',Lw1);
   text(-21,-8,'0');
   text(17,-9,'time');
   text(-25,-3,'current');

%plot current -------------------------------------------------------
   A = 0.5; Lw2 = 2;
   xi = linspace(-20,0,N/2);
   yi = A^2 ./(A^2 + (-xi-5).^2);
   yi = 8 .* yi./max(yi) - 8; yi(end) = -8;
   plot(xi,yi,'g','linewidth',Lw2);
   xj = linspace(-0,20,N/2);
   yj = A^2 ./(A^2 + (-xj+5).^2);
   yj = -8 .* yj./max(yj) - 8;yj(1) = -8;
   plot(xj,yj,'g','linewidth',Lw2);


% plot magnet -------------------------------------------------------
   xp1 = [xmax-20,xmax-15]; yp1 = [3,3];
   hp1 = area(xp1,yp1,0);
   set(hp1,'Facecolor',[1 0 0],'EdgeColor',[1 0 0]);
   xq1 = [xmax-15,xmax-10]; yq1 = [3,3];
   hq1 = area(xq1,yq1,0);
   set(hq1,'Facecolor',[0 1 0],'EdgeColor',[0 1 0]);
% arrows
   ht = text(30,2,'\leftarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[0 0 1]);    
   ht = text(16,2,'\leftarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[0 0 1]);
   text(23,-2,'B_{magnet}');
   text(16,-1.5,'N pole','Color',[0 0 1]);
   text(30,-1.5,'S pole','Color',[0 0 1]);
   xp = [xmax-10,xmax-5]; yp = [9,9];
   hp = area(xp,yp,6);
   set(hp,'Facecolor',[1 0 0],'EdgeColor',[1 0 0]);
   xq = [xmax-5,xmax]; yq = [9,9];
   hq = area(xq,yq,6);
   set(hq,'Facecolor',[0 1 0],'EdgeColor',[0 1 0]);

% set up graph window ------------------------------------------------
   axis([xmin xmax -25 10]);
   axis equal
   set(gca,'Xlim',[xmin xmax]);
   set(gca,'Ylim',[-18 15]);
   set(gca,'xcolor',[1 1 1]);
   set(gca,'ycolor',[1 1 1]);

   M = getframe;
   %M = getframe(gcf);
   [im,map] = rgb2ind(M.cdata,256,'nodither');  %RGB to indexed images
   im(1,1,1,10) = 0;

% animation ************************************************************

% move magnet ----------------------------------------------------------

dx = (abs(xmin) + xmax - 10)/N;

for c = 1 : N
% clear plotted lines for magnet
   hp = area(xp,yp,6);
   set(hp,'Facecolor',[1 1 1],'EdgeColor',[1 1 1]);
   hq = area(xq,yq,6);
   set(hq,'Facecolor',[1 1 1],'EdgeColor',[1 1 1]); 
   text(12,12,'N pole - B_{induced} - repels magnet','Color',[1 1 1]);
   ht = text(8,11,'\rightarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[1 1 1]);  
  
% plot coil
   plot([xc1,xc2],[5 5],'k','linewidth',Lw);
   plot([xc2,xc2],[5 10],'k','linewidth',Lw);
   plot([xc1,xc2],[10 10],'k','linewidth',Lw);
   plot([xc1,xc1],[5 10],'k','linewidth',Lw);
   y1 = 5*ones(20,1); y2 = y1 + 5;
  for cc = 1 : 19
     plot([x1(cc),x2(cc)],[y1(cc),y2(2)],'k','linewidth',Lw);
  end
   plot(x1,y1,'k','linewidth',Lw);
  
% plot magnet     
   xp = xp - dx;
   xq = xq - dx;
   hp = area(xp,yp,6);
   set(hp,'Facecolor',[1 0 0],'EdgeColor',[1 0 0]);
   hq = area(xq,yq,6);
   set(hq,'Facecolor',[0 1 0],'EdgeColor',[0 1 0]);    
 
% plot coil
   plot([xc1,xc2],[5 5],'k','linewidth',Lw);
   plot([xc2,xc2],[5 10],'k','linewidth',Lw);
   plot([xc1,xc2],[10 10],'k','linewidth',Lw);
   plot([xc1,xc1],[5 10],'k','linewidth',Lw);
   y1 = 5*ones(20,1); y2 = y1 + 5;
     for cc = 1 : 19
        plot([x1(cc),x2(cc)],[y1(cc),y2(2)],'k','linewidth',Lw);
     end
   plot(x1,y1,'k','linewidth',Lw)
 
% plot point on current graph
   if c < N/2, plot(xi(c),yi(c),'go'); end;
   if c > N/2, plot(xj(c-N/2),yj(c-N/2),'go'); end;
   
   if c < N/2
   ht = text(3,1.6,'\rightarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[0 1 0]);
   else
   ht = text(3,1.6,'\rightarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[1 1 1]);    
   ht = text(3,1.6,'\leftarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[0 1 0]);    
   end

% plot Binduced arrow
   if c < N/2
   text(12,12,'N pole - B_{induced} - repels magnet');
   ht = text(8,11,'\rightarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[0 0 1]);  
   else
   text(-8,12,'N pole - B_{induced} - attracts magnet');
   ht = text(-12,11,'\leftarrow');
   set(ht,'Fontsize',32,'Fontweight','bold','Color',[0 0 1]);      
   end
   
% store image
   M = getframe; 
   %M = getframe(gcf);
   im(:,:,1,c) = rgb2ind(M.cdata,map,'nodither');
   
   pause(0.01);

   
end

%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously
ag_name = 'animation.gif';
delay = 0.1;
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);


