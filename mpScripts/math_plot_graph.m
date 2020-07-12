% math_plot_graph
close all
clear all
clc


% INPUTS ================================================================
% X AXIS
   N = 2000;
   xMin = 0.0; xMax = 1.00;
   x = linspace(xMin,xMax,N);
   xL = [xMin xMax];            % X limits
   xTicks = 0:0.1:xMax;
  tT = 'time t = 0.000 s';
   
% Y AXIS
   yL = [-15.20 15.2];                    % Y limits
   
   wL = 2; A = 15;
   y = A.*sin(2*pi*x/wL);

% Labels  X Y Ttile   
  xT = 'current  I  [ A ]';  
  yT = 'voltage  V  [ V ]'; 
  
  
   %a(5) = 3; a(4)= -8; a(3) = -30; a(2) = 72; a(1) = 27;
   %y = a(5).*x.^4 + a(4).*x.^3 + a(3).*x.^2 + a(2).*x + a(1);
   %yMin = -105; yMax = 105;
   %A = 100; T = 8;
   %y = -100 .* sin(2*pi*x/T);

   
 %% GRAPHICS =============================================================
    figure(1);    
    set(gcf,'units','normalized','position',[0.1 0.1 0.30 0.30]);
    fs = 14;
    col = 'b'; LW = 2;
    yP = [0 2 4.1 5.9 8 9.8];
    xP = [0 0.1 0.2 0.3 0.4 0.5];
    plot(xP,yP,'bo');
    hold on
    xP = [0 0.5];
    yP = [0 10];
    plot(xP,yP,col,'linewidth',LW)
   
   % yP = [0 -1];
   % col = 'r';
   % plot(xP,yP,col,'linewidth',LW)
    xlabel(xT);
    ylabel(yT);
    %title(tT,'fontweight','normal');
%     set(gca,'yLim',yL);
%     set(gca,'xLim',xL);
%     set(gca,'xTick',xTicks);
    box on; grid on;
%   legend('A','B');
    
    set(gca,'fontsize',fs);
% % plot points
%      X = 1;  Y = 1; col = 'b';
%      h_plot = plot(X,Y,'o');
%      set(h_plot,'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',6); 
   