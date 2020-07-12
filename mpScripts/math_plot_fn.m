function math_plot_fn
% math_plot_graph
close all
clear all
clc

% define X values
  N = 2000;
  xMin = 5; xMax = 30;
  x = linspace(xMin,xMax,N);

% define function
  y = fn(x);
  
  yMin = -0.5; yMax = 0.5;
  
  
 % plot function
    tx = ' x';
    ty = ' y';
    tm = 'y  =  x^3 - 2x^2 - x + 2';
    fs = 14; 
   set(gcf,'units','normalized','position',[0.1 0.1 0.30 0.30]);
   set(gca,'fontsize',fs);
    col = 'b'; LW = 3;
   X = x;   Y = y;
   plot(X,Y,col,'linewidth',LW)
   xlabel(tx); ylabel(ty);
   %title(tm);
   %set(gca,'yLim',[yMin yMax]);
   
   hold on
   X = x;   Y = -y;
   plot(X,Y,col,'linewidth',LW)
   
   X = -x;   Y = -y;
   plot(X,Y,col,'linewidth',LW)
   
   X = -x;   Y = y;
   plot(X,Y,col,'linewidth',LW)
   
   % tangent 
   M1 = 4/(5*sqrt(3)); B1 = -2/sqrt(3);
   xT = linspace(-xMax,xMax, 200);
   yT = M1 .* xT + B1;
   
   xN = linspace(-2,xMax, 200);
   M2 = -1/M1; B2 = 14.5*sqrt(3);
   yN = M2 .* xN + B2;
   
   X = xT;   Y = yT;  col = 'r'; LW = 1;
   plot(X,Y,col,'linewidth',LW)
   
   X = xN;   Y = yN;  col = 'r'; LW = 1;
   plot(X,Y,col,'linewidth',LW)
   
   a = 5; b = 2; c = sqrt(a^2+b^2);
   xP = 10; yP = fn(xP);
   
   X = [a^2/c a^2/c];   Y = [-30 30];  col = 'k'; LW = 1;
   plot(X,Y,col,'linewidth',LW)
   
   X = [-a^2/c -a^2/c];   Y = [-30 30];  col = 'k'; LW = 1;
   plot(X,Y,col,'linewidth',LW)
   
   % asymptotoes
   X = xT;   Y = (b/a).*xT;  col = 'k'; LW = 1;
   plot(X,Y,col,'linewidth',LW)
   plot(X,-Y,col,'linewidth',LW)
    
% plot points
     hold on
       X = [-a a -c c];  Y = [0 0 0 0]; col = 'b';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',3);

       X = [xP 0 0];  Y = [yP B1 B2]; col = 'r';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',5);
%      
      X = [-B1/M1 -B2/M2];  Y = [0 0]; col = 'm';
      h_plot = plot(X,Y,'o');
      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',5);
%      
%      X = -1/sqrt(2);  Y = fn(X); col = 'r';
%      h_plot = plot(X,Y,'o');
%      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6);
%      
%      X = -1/sqrt(2);  Y = -fn(X); col = 'r';
%      h_plot = plot(X,Y,'o');
%      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6);
%      
%      X = 1/sqrt(2);  Y = fn(X); col = 'r';
%      h_plot = plot(X,Y,'o');
%      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6);
%      
%      X = 1/sqrt(2);  Y = -fn(X); col = 'r';
%      h_plot = plot(X,Y,'o');
%      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6);
%      
set(gca,'xLim',[-30 30]);
set(gca,'yLim',[-30 30]);
% set(gca,'yTick',[yMin :0.1: yMax]);
 % fn(3)   
axis square 
box on; grid on; 
     
function y = fn(x)
    %a(5) = 3; a(4)= -8; a(3) = -30; a(2) = 72; a(1) = 27;
    %y = a(5).*x.^4 + a(4).*x.^3 + a(3).*x.^2 + a(2).*x + a(1);    
  % y = x.^2 .* (1-x.^2);
  
  y = sqrt(4.*x.^2-100) ./5;
   