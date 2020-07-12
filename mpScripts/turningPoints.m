function [indexMin indexMax] = turningPoints(xData, yData) 

size = length(xData);    %Get the length of the dataset

a1 = yData(1,1); a2 = yData(1,2);
if a2 > a1, flag = 1; else flag = 2; end;
v = 0;

% find max
for x = 2:size-1
    a1 = yData(1,x);     %Get two adjacent samples of the dataset
    a2 = yData(1,x+1);

    if flag == 1 && a2 > a1; x = x+1; end;
    if (flag == 1 && a2 < a1); v = v + 1; indexMax(v) = x; x = x+1; end;
    if a2 <= a1, flag = 0; end;
    if a2 > a1, flag = 1; end;
end

a1 = yData(1,1); a2 = yData(1,2);
if a2 < a1, flag = 1; else flag = 2; end;
v = 0;

% find min
for x = 2:size-1
    a1 = yData(1,x);     %Get two adjacent samples of the dataset
    a2 = yData(1,x+1);

    if flag == 1 && a2 < a1; x = x+1;
    end;
    if (flag == 1 && a2 > a1); v = v + 1; indexMin(v) = x; x = x+1; end;
    if a2 >= a1, flag = 0; end;
    if a2 < a1, flag = 1; end;
end

figure(99)
fs = 10;
set(gca,'fontsize',fs);
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.15 0.15]);
plot(xData,yData,'lineWidth',2)
hold on

hp1 = plot(xData(indexMax),yData(indexMax),'o');
set(hp1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);

hp1 = plot(xData(indexMin),yData(indexMin),'s');
set(hp1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',4);


% disp('  ');
% disp('Max values - indices / xData / yData');
% indexMax;
% xData(indexMax);
% yData(indexMax);
% 
% disp('  ');
% disp('  ');
% disp('Min values - indices / xData / yData');
% indexMin;
% xData(indexMin);
% yData(indexMin);

 
