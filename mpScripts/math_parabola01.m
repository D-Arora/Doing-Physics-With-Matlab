%%   ellipse ===============================================================================================================
clear all
close all
clc

a = 2;
xP = -8;


xMax = 10; xMin = -xMax;
yMax = 10; yMin = -yMax;
dt = 2;  % tick spacing
N = 5597;
x0 = 0; y0 = 0;
x = linspace(xMin,xMax,N);

y = y0 + (x-x0).^2 ./(4*a);

yP = y0 + (xP-x0).^2 ./(4*a);


% distance FP
xF = x0; yF = y0+a;
dFP = sqrt((xP-xF)^2+(yP-yF)^2);

% distance DP
xD = xP; yD = y0-a;
dDP = sqrt((xP-xD)^2+(yP-yD)^2);

figure(1)   % ===========================================================
   fs = 14; 
   set(gcf,'units','normalized','position',[0.1 0.1 0.30 0.30]);
   set(gca,'fontsize',fs);
   
   % plot hyperbola
      X = x;   Y = y; col = 'b'; LW = 2;
      plot(X,Y,col,'lineWidth',LW);
      hold on
     
     % plot(X,Y,col,'lineWidth',LW);
   
  % plot X and Y axes
      X = [-xMax xMax];   Y = [0 0]; col = 'k'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      X = [0 0];   Y = [-yMax yMax]; col = 'k'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      
%  % plot y = x
%      X = x;  Y = x; col = [0.8 0.8 0.8]; LW = 1;
%      h_plot = plot(X,Y,'lineWidth',LW);
%      set(h_plot,'color',col);
%      X = -x;  Y = -x;
%      h_plot = plot(X,Y,'lineWidth',LW); 
%       set(h_plot,'color',col);
      
  % plot fixed point and vertix
     X = x0;  Y = y0+a; col = 'r';
     h_plot = plot(X,Y,'o');
     set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6); 
     
      X = x0;  Y = y0; col = 'b';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6); 
     
 % plot directrix    
     X = [xMin xMax];  Y = [y0-a y0-a]; col = 'r'; LW = 1;
     h_plot = plot(X,Y,col,'linewidth',LW');
     
   grid on; box on;
   
   % LABELS
     xlabel('x','fontsize',fs); ylabel('y','fontsize',fs);
     t1 = 'x y = a^2 / 2   y = (a^2/2) / x     a = ';
     t2 = num2str(a,2);
    
   %  t4 = num2str(b,2);
     tm = [t1 t2]; fs = 12;
   %  title(tm,'color','b','fontSize',fs);
     
     
   axis equal
   set(gca,'xLim',[-xMax xMax]);
   set(gca,'yLim',[-yMax yMax]);
   set(gca,'xTick',-xMax:dt:xMax);
   set(gca,'yTick',yMin:dt:yMax);
 
% figure(1)   % ===========================================================
%    fs = 12; LW = 1;
%    set(gcf,'units','normalized','position',[0.1 0.1 0.70 0.60]);
%    set(gca,'fontsize',fs);
%    
%    handle = subplot(1,2,1);
%    set(handle, 'OuterPosition', [0.01 0 0.5 1]) ;  
%    X = [0 100];   Y = [0 0]; col = 'k'; LW = 1;
%    plot(X,Y,col,'lineWidth',LW);
%    hold on
%    
%    
%    nD = 3;
%    
%    X = [0 100];   Y = [100 100];
%    plot(X,Y,col,'lineWidth',LW);
%    
%    X = [100 100];   Y = [100 100];
%    plot(X,Y,col,'lineWidth',LW);
%    
%    X = [100 100];   Y = [0 0];
%    plot(X,Y,col,'lineWidth',LW);
%    
%   
%   
%    X = 5; Y = 95;  tt1 = 'P(x, y) = (';  tt2 = num2str(xP,nD);
%        tt3 = ', ';  tt4 = num2str(yP,nD+2); tt5 = ')'; tt= [tt1 tt2 tt3 tt4 tt5];  text(X,Y,tt);
%   
%    
%    X =  5; Y = 90;  tt1 = 'semi-major axis   a  =  ';  tt2 = num2str(a,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    X = 55; Y = 90;  tt1 = 'major axis   2a  =  '; tt2 = num2str(2*a,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    
%    
%    X =  5; Y = 85;  tt1 = 'semi-minor axis   b  =  ';  tt2 = num2str(b,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    X = 55; Y = 85;  tt1 = 'minor axis   2b  =  ';  tt2 = num2str(2*b,nD);    tt= [tt1 tt2];  text(X,Y,tt);
%    
%    X =  5; Y = 80;  tt1 = 'focal length [OF_1  OF_2]  c  =  ';  tt2 = num2str(c,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    X = 55; Y = 80;  tt1 = 'eccentricity   e  =  ';  tt2 = num2str(e,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    
%    X =  5; Y = 75;  tt1 = 'directrices 1: x =   ';  tt2 = num2str(-D,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    X = 55; Y = 75;  tt1 = 'directrices 2: x =   ';  tt2 = num2str(D,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    
%    X =  5; Y = 70;  tt1 = 'slope tangent M_1 =   ';  tt2 = num2str(M1,nD);  tt= [tt1 tt2]; text(X,Y,tt);
%    X = 55; Y = 70;  tt1 = 'slope normal M_2 =   ';  tt2 = num2str(M2,nD); tt= [tt1 tt2]; text(X,Y,tt);
%    
%    X =  5; Y = 65;  tt1 = 'intercept tangent B_1 =   ';  tt2 = num2str(B1,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    X = 55; Y = 65;  tt1 = 'intercept normal B_2 =   ';  tt2 = num2str(B2,nD);  tt= [tt1 tt2]; text(X,Y,tt);
%    
%    X =  5; Y = 60;  tt1 = 'T  tangent cross X-axis: x = ';  tt2 = num2str(xtc,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%    X = 55; Y = 60;  tt1 = 'N normal cross X-axis: x = ';  tt2 = num2str(xnc,nD);  tt= [tt1 tt2];  text(X,Y,tt);
%        
%    X =  5;  Y = 55;  tt1 = 'distances:   PF_1 = ';   tt2 = num2str(PF1,nD);   tt= [tt1 tt2];   text(X,Y,tt);
%    X = 45; Y = 55;   tt1 = 'PF_2 = ';   tt2 = num2str(PF2,nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%    X = 65; Y = 55;   tt1 = 'PF_1 + PF_2 = ';   tt2 = num2str(PF1+PF2,nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%   
%    X =  5; Y = 50;  tt1 = 'distances:   PF_1 = ';   tt2 = num2str(PF1,nD);   tt= [tt1 tt2];   text(X,Y,tt);
%    X = 45; Y = 50;  tt1 = 'PD_1 = ';   tt2 = num2str(PD1,nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%    X = 65; Y = 50;  tt1 = 'PF_1 / PD_1 = ';   tt2 = num2str(PF1/PD1,nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%    
%    X =  5; Y = 45;  tt1 = 'distances:   PF_2 = ';   tt2 = num2str(PF2,nD);   tt= [tt1 tt2];   text(X,Y,tt);
%    X = 45; Y = 45;  tt1 = 'PD_2 = ';   tt2 = num2str(PD2,nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%    X = 65; Y = 45;  tt1 = 'PF_2 / PD_2 = ';   tt2 = num2str(PF2/PD2,nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%    
%    X =  5; Y = 40;  tt1 = 'angles:';   tt2 = '   \theta_1 = NPF_1   \theta_2 = NPF_2   ';   tt= [tt1 tt2];   text(X,Y,tt); 
%    X =  5; Y = 35;  tt1 = 'angles:';   tt2 = '   \theta_3 = NF_1P   \theta_4 = PNF_2   \theta_5 = PF_2T   ';   tt= [tt1 tt2];   text(X,Y,tt); 
%    X =  5; Y = 30;  tt1 = 'angles:';   tt2 = '   \theta_1 = \theta_4 - \theta_3       \theta_2 = \theta_5 - \theta_4     ';   tt= [tt1 tt2];   text(X,Y,tt); 
%    
%   X = 5; Y = 25;  tt1 = 'angle [degrees]   \theta_4 = ';   tt2 = num2str(theta(4),nD);   tt= [tt1 tt2];   text(X,Y,tt);
%   X =  45; Y = 25;  tt1 = '     \theta_3 = ';   tt2 = num2str(theta(3),nD);   tt= [tt1 tt2];   text(X,Y,tt);
%   X =  65; Y = 25;  tt1 = '     \theta_1 = ';   tt2 = num2str(theta(1),nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%   
%    X = 5; Y = 20;  tt1 = 'angle [degrees]   \theta_4 = ';   tt2 = num2str(theta(4),nD);   tt= [tt1 tt2];   text(X,Y,tt);
%   X =  45; Y = 20;  tt1 = '     \theta_5 = ';   tt2 = num2str(theta(5),nD);   tt= [tt1 tt2];   text(X,Y,tt);
%   X =  65; Y = 20;  tt1 = '     \theta_2 = ';   tt2 = num2str(theta(2),nD);   tt= [tt1 tt2];   text(X,Y,tt); 
%   
%   axis off
%    
%    subplot(1,2,2); % ======================================================
%       tx = 'x';  ty = 'y'; 
%    
%   % ellipse
%       X = x;   Y = y; col = 'b'; LW = 2;
%       plot(X,Y,col,'lineWidth',LW);
%       hold on
%   % tangent 
%      X = x1; Y = y1; col = 'r'; LW = 1.5;
%      plot(X,Y,col,'lineWidth',LW);
%   
%   % normal
%       X = x1; Y = y2; col = 'm'; LW = 1.5;
%       plot(X,Y,col,'lineWidth',LW);
%   
%   % focal points 
%       X = [-c c]; Y = [0 0]; col = 'r';
%       h_plot = plot(X,Y,'o');
%       set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',5);
%       
%  % directrix P T N O points     
%       X = [-D D xP xtc xnc 0]; Y = [yP yP yP 0 0 0]; col = 'k';
%       h_plot = plot(X,Y,'o');
%       set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',4); 
%  
%  
%   % Line FP
%    X = [-c xP c]; Y = [0 yP 0]; col = 'k'; LW = 1;
%      plot(X,Y,col,'lineWidth',LW);
%   
%   % Directrix
%     X = [-D -D]; Y = [yMin yMax]; col = 'b'; LW = 1;
%     plot(X,Y,col,'lineWidth',LW);
%     X = [D D]; Y = [yMin yMax]; col = 'b'; LW = 1;
%     plot(X,Y,col,'lineWidth',LW); 
%  
%   % line D1 D2
%   X = [-D D]; Y = [yP yP]; col = 'k'; LW = 0.3;
%     plot(X,Y,col,'lineWidth',LW); 
%   
%     
%   % F1 F2 P D T N
%   X = -c-0.3; Y = -0.8; L = 'F_1'; text(X,Y,L);
%   X = c-0.4; Y = -0.8; L = 'F_2';  text(X,Y,L);
%   X = 4.5; Y = 2; L = 'P';         text(X,Y,L);
%   X = 6.5; Y = 0.5; L = 'T';       text(X,Y,L);
%   X = 2; Y = 0.5; L = 'N';         text(X,Y,L);
%   X = -D-1.1; Y = yP; L = 'D_1';   text(X,Y,L);
%   X = D+0.5;  Y = yP; L = 'D_2';   text(X,Y,L);
%   X = 0;  Y = -0.5; L = 'O';       text(X,Y,L);
%   X = -3.5;  Y = 9; L = 'tangent';   text(X,Y,L,'color','r');
%   X = -3.5;  Y = -9; L = 'normal';   text(X,Y,L,'color','m');
%   
%   xlabel(tx); ylabel(ty);
%     
%    
%  
%   % axis equal
% % ax = gca;
%  
%   axis equal
% %  set(gca,'YMinorGrid','on');
% %  set(gca,'MinorGridAlpha', 0.5);
% %  set(gca,'MinorGridLineStyle','-');
% %  set(gca,'MinorGridColor',[0.8 0.8 0.8]);
% %  set(gca,'GridColor','k');
%  %set(gca,'YMinorGrid','on');
% %set(gca,'xTick',xMin:5:xMax);
%  
%  %ax.XTick = [-6*pi/6,-5*pi/6,-4*pi/6,-3*pi/6, -2*pi/6, -pi/6,0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,6*pi/6];
%  %ax.XTickLabel = {'-\pi','','-2\pi/3','','-\pi/3','','0','','\pi/3','','2\pi/3','','\pi'}; 
%  %ax.xLim = [-pi, pi];
%  %ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
% %legend('|x-3|','|x - 5|','|x-3|+|x-5|', 'y = 10');
%    
% set(gca,'xLim',[xMin xMax]);
% set(gca,'yLim',[yMin yMax]);
% set(gca,'xTick',xMin:2:xMax);
% grid on; box on;
%      

% ===========================================================================================================
figure(2)   
   fs = 14; 
   set(gcf,'units','normalized','position',[0.1 0.1 0.30 0.30]);
   set(gca,'fontsize',fs);
   
   % plot hyperbola
      X = x;   Y = y; col = 'b'; LW = 2;
      plot(X,Y,col,'lineWidth',LW);
      hold on
     
     % plot(X,Y,col,'lineWidth',LW);
   
  % plot X and Y axes
      X = [-xMax xMax];   Y = [0 0]; col = 'k'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      X = [0 0];   Y = [-yMax yMax]; col = 'k'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      
  % plot F   A   P  D
     X = x0;  Y = y0+a; col = 'r';
     h_plot = plot(X,Y,'o');
     set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6); 
     
      X = x0;  Y = y0; col = 'b';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6); 
      
      X = xP;  Y = yP; col = 'b';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6); 
      
      X = xD;  Y = yD; col = 'r';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6);  
      
     
 % plot directrix    
     X = [xMin xMax];  Y = [y0-a y0-a]; col = 'r'; LW = 1;
     h_plot = plot(X,Y,col,'linewidth',LW');
     
% plot FP
     X = [xP xF];  Y = [yP yF]; col = 'r'; LW = 1;
     h_plot = plot(X,Y,col,'linewidth',LW');

 % plot FD
     X = [xP xD];  Y = [yP yD]; col = 'r'; LW = 1;
     h_plot = plot(X,Y,col,'linewidth',LW');
     
   grid on; box on;
   
   % LABELS
     xlabel('x','fontsize',fs); ylabel('y','fontsize',fs);
     fs = 10;
     t1 = 'x_P  =  ';         t2 = num2str(xP,3);
     t3 = '   y_P  =  ';      t4 = num2str(yP,3);
     t5 = '   d_{FP}  =  ';      t6 = num2str(dFP,3);
     t7 = '   d_{DP}  =  '; t8 = num2str(dDP,3);
     tm = [t1 t2 t3 t4 t5 t6 t7 t8 ]; 
     title(tm,'color','b','fontSize',fs);
    
     X = xP+1; Y = yP; t1 = 'P'; text(X,Y,t1);  
     X = xF-1; Y = yF+1; t1 = 'F'; text(X,Y,t1);    
     X = x0; Y = y0-1; t1 = 'A'; text(X,Y,t1); 
     X = xD+0.5; Y = yD-1; t1 = 'D'; text(X,Y,t1); 
   
   axis equal
   set(gca,'xLim',[-xMax xMax]);
   set(gca,'yLim',[-yMax yMax]);
   set(gca,'xTick',-xMax:dt:xMax);
   set(gca,'yTick',yMin:dt:yMax);
 
