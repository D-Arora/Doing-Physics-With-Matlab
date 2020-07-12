close all
clear all
clc


x= 1:10
y1 = x.^2;
y2 = 2.*x+1;



figure(1)                      % open figure window for plotting
set(gcf,'Name','Sinusoidal Function')
set(gcf,'NumberTitle','off')
set(gcf,'Color',[0.8 0.8 0])
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1 0.1 0.6 0.7])
set(gcf,'PaperType','A4')

%plot titles
%tm = 'Sinusodal Function';     % Main title
%tx = 'time t';                 % X-axis title
%ty = 'position d'              % Y-axis title

%plot

subplot(2,1,2)
%p1 = plot(x,y2,'LineWidth',2);        % XY plot
set(gca,'Position',[0.15 0.112 0.7 0.2]);
set(gca,'XLim',[0 10])
set(gca,'YLim',[0 10])
text(2,8,'slope ')
axis off




subplot(2,1,1)
plot(x,y1,'LineWidth',2)        % XY plot
set(gca,'Position',[0.15 0.4 0.7 0.5])
grid on
%title(tm)
%xlabel(tx)
%ylabel(ty)

%grid on
