%%   ellipse ===============================================================================================================
clear all
close all
clc

a = sqrt(2);
b = sqrt(2);
tMin = pi;
tMax= 2*pi;

xMax = 6; xMin = -xMax;
yMax = 6; yMin = -yMax;

N = 5597;

t = linspace(tMin,tMax,N);
x = linspace(a,xMax,N);
x1 = linspace(-xMax,xMax,N);
%x = a .* cos(t);
%y = b.* sin(t);
%x = a .* sec(t);
%y = b .* tan(t); 
y = b .* sqrt(x.^2 ./ a^2 - 1);


% xP = -3.8;
% yP = sqrt(b^2 - (b/a)^2 .* xP^2);
% 
% x1 = linspace(xMin,xMax,N);
% M1 = -(b/a)^2 * (xP/yP);
% B1 = yP + (b/a)^2 * (xP^2/yP);
% 
% y1 = M1 .* x1 + B1;
% 
% M2 = -1/M1;
% B2 = yP - M2 * xP;
% 
% y2 = M2.* x1 + B2;
% 
% % focal length c
%    c = sqrt(a^2-b^2);
% % eccentricity e 
%    e = c/a;
% % directrix   D
%    D = a/e; 
% % y = 0 tangent crosses X-axis
%     xtc = -B1/M1;
% % y = 0 normal crosses X-axis
%     xnc = -B2/M2;
% % angle of tangent
%     theta1 = atand(M1);
% % angle of normal
%     theta2 = atand(M2);  
% 
%     
% % distances
%    X1 = -c; Y1 = 0; X2 = xP; Y2 = yP;
%    PF1 = sqrt((X1-X2)^2 + (Y1-Y2)^2);
%    
%    X1 = c; Y1 = 0; X2 = xP; Y2 = yP;
%    PF2 = sqrt((X1-X2)^2 + (Y1-Y2)^2);
%    
%    X1 = -D; Y1 = yP; X2 = xP; Y2 = yP;   PD1 = sqrt((X1-X2)^2 + (Y1-Y2)^2);
%      
%    X1 = D; Y1 = yP; X2 = xP; Y2 = yP;
%    PD2 = sqrt((X1-X2)^2 + (Y1-Y2)^2);
%    
%  % angles
%  theta(4) = atand(M2); if theta(4) < 0, theta(4) = 180+theta(4);end;
%  theta(3) = atand((yP-0)/(xP+c));  
%  theta(1) = theta(4)-theta(3); if theta(3) < 0, theta(3) = 180+theta(3);end;
%  
%  theta(5) = atand((yP-0)/(xP-c)); if theta(5) < 0, theta(5) = 180+theta(5);end;
%  theta(2) = theta(5)-theta(4);
 
figure(1)   % ===========================================================
   fs = 14; 
   set(gcf,'units','normalized','position',[0.1 0.1 0.30 0.30]);
   set(gca,'fontsize',fs);
   
   % plot hyperbola
      X = x;   Y = y; col = 'b'; LW = 2;
      plot(X,Y,col,'lineWidth',LW);
      hold on
      X = x;   Y = -y; col = 'b'; LW = 2;
      plot(X,Y,col,'lineWidth',LW);
      X = -x;   Y = -y; col = 'b'; LW = 2;
      plot(X,Y,col,'lineWidth',LW);
      X = -x;   Y = y; col = 'b'; LW = 2;
      plot(X,Y,col,'lineWidth',LW);
   
  % plot X and Y axes
      X = [-xMax xMax];   Y = [0 0]; col = 'k'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      X = [0 0];   Y = [-yMax yMax]; col = 'k'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      
  % plot asymptotes
      X = x1;   Y = (b/a) .* x1; col = 'r'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);
      X = x1;   Y = -(b/a) .* x1; col = 'r'; LW = 1;
      plot(X,Y,col,'lineWidth',LW);    
      
   grid on; box on;
   
   % LABELS
     xlabel('x','fontsize',fs); ylabel('y','fontsize',fs);
     t1 = 'x^2 / a^2 - y^2 / b^2 = 1     a = ';
     t2 = num2str(a,2);
     t3 = '     b = ';
     t4 = num2str(b,2);
     tm = [t1 t2 t3 t4]; fs = 12;
    % title(tm,'color','b','fontSize',fs);
     
     
   axis equal
   set(gca,'xLim',[xMin xMax]);
   set(gca,'yLim',[yMin yMax]);
   set(gca,'xTick',xMin:2:xMax);
   set(gca,'yTick',yMin:2:yMax);
 
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