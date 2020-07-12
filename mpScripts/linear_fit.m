function [a1, a2, Ea1, Ea2, r] = linear_fit(xyData,xmin, xmax, flag)

% Ian Cooper
% Doing Physics with Matlab

% Function to determine the coefficients a1 and a2 
%    for a linear fit to experimental data
% Flag = 1: Linear relationship y = mx + b = a1 + a2 x  a1 = b and a2 = m
% Flag = 2: Power relationship y = a1 x^a2
% Flag = 3: Exponential relationship y = a1 exp(a2)

% Inputs to function
%    (x,y) data in a matrix (n x 2) with column 1 for x data & column 2 for y data
%    Graph min and max x values are xmin & xmax
%    Flag to select type of relationship

% Outputs
%    a1 and a2 for the coefficients as determined by the value of the flag
%    Ea1 and Ea2 the uncertainties for a1 and a2
%    Correclation coefficent r

close all                    % close all Plot Windows

% Input Graph Titles --------------------------------------------------

switch flag

    case 1
     tm = '';                      % Title for plot
     tx1 = 'extension  e (mm)';    % X-axis label        
     ty1 = 'load  F  (N)';         % Y-axis label
%      disp
%      disp('Enter the main title for the plot and labels for the axes');     
%      tm = input('Title for plot   ');
%      tx1 = input('X-axis label    ');
%      ty1 = input('Y-axis label    ');
    
    case 2
     tm = '';                      % Title for plot
     tx1 = 'load  m  (kg)';    % X-axis (linear scale) label        
     ty1 = 'period  T (s)  (N)';         % Y-axis (linear scale) label 
%      disp
%      disp('Enter the main title for the plot and labels for the axes');      
%      tm = input('Title for plot   ');
%      tx1 = input('X-axis (linear scale)label    ');
%      ty1 = input('Y-axis (label scale) label    ');
    
     tx2 = 'log_{10}(m)';         % X-axis (log10 scale) label        
     ty2 = 'log_{10}(T)';          % Y-axis (log10 scale) label 
%    tx2 = input('X-axis (log10 scale)label    ');
%    ty2 = input('Y-axis (log10 scale)label    ');
    
    case 3    
     tm = '';                      % Title for plot
     tx1 = 'time  t (s)';    % X-axis (linear scale) label        
     ty1 = 'amplitude  A (mm)';         % Y-axis (linear scale) label 
%      disp
%      disp('Enter the main title for the plot and labels for the axes');      
%      tm = input('Title for plot   ');
%      tx1 = input('X-axis (linear scale)label    ');
%      ty1 = input('Y-axis (label scale) label    ');
    
     tx2 = 'time  t (s)';         % X-axis (log10 scale) label        
     ty2 = 'log_{e}(A)';          % Y-axis (log10 scale) label 
%    tx2 = input('X-axis (log10 scale)label    ');
%    ty2 = input('Y-axis (log10 scale)label    ');

end


% Process input data --------------------------------------------------
switch flag

case 1                       % linear 
x = xyData(:,1);
y = xyData(:,2);

case 2                       % power 
x = log10(xyData(:,1));
y = log10(xyData(:,2));

case 3                       % exponentail 
x = xyData(:,1);
y = log(xyData(:,2));
end 

n = length(x);               % number of data points

% Calculations for slope and intercept y = m x + b --------------------
 sx = sum(x);
 sy = sum(y);
 sxx = sum(x .* x);
 syy = sum(y .* y);
 sxy = sum(x .* y);
 sx2 = sx * sx;
 sy2 = sy * sy;
 
 m = (n * sxy - sx * sy)/(n * sxx - sx2);
 
 b = (1/n) * (sy - m * sx);
 
 s = sqrt((syy - sy2/n - m * (sxy - sx * sy/n))/(n-2));
 
 Em = s * sqrt(n / (n * sxx - sx2));
 
 Eb = s * sqrt(sxx / (n * sxx - sx2));
 
 r = (sxy - sx * sy /n) / sqrt((sxx - sx2/n) * (syy - sy2/n));
 
 
% Calculations for coeffficients a1 and a2 & fitted functions 
%  and text labelling for graphs --------------------------------------
nx = 400;
xf = linspace(xmin, xmax, nx);

t31 = 'correlation coefficient  r  =  ';
t32 = num2str(r, '%0.4g');
text3 = [t31 t32];          % correlation coeff r


% Graphical output ---------------------------------------------------

switch flag
case 1
     a1 = b; Ea1 = Eb;
     a2 = m; Ea2 = Em;
     yf = m .* xf + b;
     
     %labels  
     t11 = 'Intercept     a_1 = ';
     t12 = num2str(b, '%0.4g');
     t13 = '       E(a_1) = ';
     t14 = num2str(Eb, '%0.2g');
     text1 = [t11 t12 t13 t14];
     t21 = 'slope          a_2 = ';
     t22 = num2str(m, '%0.4g');
     t23 = '     E(a_2) = ';
     t24 = num2str(Em, '%0.2g');
     text2 = [t21 t22 t23 t24];
        
     figure(1)
     set(gcf,'Units','Normalized')
     set(gcf,'Position',[0.1 0.1 0.5 0.6])
     set(gcf,'PaperType','A4')
     
     subplot(2,1,2)          % SUMMARY
     set(gca,'Position',[0.15 0.02 0.75 0.2]);
     set(gca,'XLim',[0 10])
     set(gca,'YLim',[0 10])
     text(1,10,'y = a_1 + a_2 x','fontSize',12);
     text(1,7.5,text1,'fontSize',12);
     text(1,5,text2,'fontSize',12);
     text(1,2.5,text3,'fontSize',12);
     axis off
         
     subplot(2,1,1)          % LINEAR PLOT 
     plot(xf,yf,'LineWidth',2);
     hold on
     plot(x,y,'o');
     set(gca,'Position',[0.15 0.4 0.75 0.55])
     set(gca,'FontSize',12)
     grid on
     xlabel(tx1,'FontSize',12);
     ylabel(ty1,'FontSize',12);
     title(tm);

   
case 2
     a1 = 10^b; Ea1 = 10^b * Eb;
     a2 = m; Ea2 = Em; 
     yf = a1 .* xf.^a2;
     %labels  
     t11 = 'coefficient     a_1 = ';
     t12 = num2str(a1, '%0.4g');
     t13 = '       E(a_1) = ';
     t14 = num2str(Ea1, '%0.2g');
     text1 = [t11 t12 t13 t14];        % a1
     t21 = 'coefficient     a_2 = ';
     t22 = num2str(a2, '%0.4g');
     t23 = '     E(a_2) = ';
     t24 = num2str(Ea2, '%0.2g');
     text2 = [t21 t22 t23 t24];        % a2
    
     figure(1)               
     set(gcf,'Units','Normalized')
     set(gcf,'Position',[0.1 0.1 0.5 0.6])
     set(gcf,'PaperType','A4')
     
     subplot(2,2,3)          % SUMMARY        
     set(gca,'Position',[0.15 0.02 0.75 0.2]);
     set(gca,'XLim',[0 10])
     set(gca,'YLim',[0 10])
     text(1,10,'y = a1 * x^{a2}','fontSize',12);
     text(1,7.5,text1,'fontSize',12);
     text(1,5,text2,'fontSize',12);
     text(1,2.5,text3,'fontSize',12);
     axis off
         
     subplot(2,2,2)          % LINEAR PLOT        
     plot(xf,yf,'LineWidth',2);
     hold on
     plot(xyData(:,1),xyData(:,2),'o');
     set(gca,'Position',[0.60 0.4 0.35 0.55])
     set(gca,'FontSize',12)
     set(gca,'XLim',[xmin xmax]);
     grid on
     xlabel(tx1,'FontSize',12);
     ylabel(ty1,'FontSize',12);
     title(tm);
     
     subplot(2,2,1)          % LOG LOG PLOT
     plot(log10(xf),log10(yf),'LineWidth',2);
     hold on
     plot(x,y,'o');
     set(gca,'Position',[0.1 0.4 0.35 0.55])
     set(gca,'FontSize',12)
     grid on
     xlabel(tx2,'FontSize',12);
     ylabel(ty2,'FontSize',12);
     title(tm);
     
         
case 3
    a1 = exp(b); Ea1 = exp(b) * Eb;
    a2 = m; Ea2 = Em; 
    
    yf = a1 .* exp(a2 .* xf);
 
    %labels  
     t11 = 'coefficient     a_1 = ';
     t12 = num2str(a1, '%0.4g');
     t13 = '       E(a_1) = ';
     t14 = num2str(Ea1, '%0.2g');
     text1 = [t11 t12 t13 t14];        % a1
     t21 = 'coefficient     a_2 = ';
     t22 = num2str(a2, '%0.4g');
     t23 = '     E(a_2) = ';
     t24 = num2str(Ea2, '%0.2g');
     text2 = [t21 t22 t23 t24];        % a2
    
     figure(1)               
     set(gcf,'Units','Normalized')
     set(gcf,'Position',[0.1 0.1 0.5 0.6])
     set(gcf,'PaperType','A4')
     
     subplot(2,2,3)          % SUMMARY        
     set(gca,'Position',[0.15 0.02 0.75 0.2]);
     set(gca,'XLim',[0 10])
     set(gca,'YLim',[0 10])
     text(1,10,'y = a1 * exp(a2 * x)','fontSize',12);
     text(1,7.5,text1,'fontSize',12);
     text(1,5,text2,'fontSize',12);
     text(1,2.5,text3,'fontSize',12);
     axis off
         
     subplot(2,2,2)          % LINEAR PLOT        
     plot(xf,yf,'LineWidth',2);
     hold on
     plot(xyData(:,1),xyData(:,2),'o');
     set(gca,'Position',[0.55 0.4 0.40 0.55])
     set(gca,'FontSize',12)
     set(gca,'XLim',[xmin xmax]);
     grid on
     xlabel(tx1,'FontSize',12);
     ylabel(ty1,'FontSize',12);
     title(tm);
     
     subplot(2,2,1)          % LOG LINEAR PLOT
     plot(xf,log(yf),'LineWidth',2);
     hold on
     plot(x,y,'o');
     set(gca,'Position',[0.1 0.4 0.35 0.55])
     set(gca,'FontSize',12)
     grid on
     xlabel(tx2,'FontSize',12);
     ylabel(ty2,'FontSize',12);
     title(tm); 
    
end
 
 % text output to Command Window --------------------------------------
 switch flag
     
 case 1
    disp('y = m x + b') 
    fprintf('n =  %0.4g \n', n)
    fprintf('slope m =  %0.4g \n', m)
    fprintf('intercept b =  %0.4g \n', b)
    fprintf('Em =  %0.4g \n', Em)
    fprintf('Eb =  %0.4g \n', Eb)
    fprintf('correlation r =  %0.4g \n', r)
    
case 2
    disp('y = a1 x ^(a2)') 
    fprintf('n =  %0.4g \n', n)
    fprintf('a1 =  %0.4g \n', a1)
    fprintf('a2 =  %0.4g \n', a2)
    fprintf('Ea1 =  %0.4g \n', Ea1)
    fprintf('Ea2 =  %0.4g \n', Ea2)
    fprintf('correlation r =  %0.4g \n', r)
 
case 3
    disp('y = a1 exp(a2 * x)') 
    fprintf('n =  %0.4g \n', n)
    fprintf('a1 =  %0.4g \n', a1)
    fprintf('a2 =  %0.4g \n', a2)
    fprintf('Ea1 =  %0.4g \n', Ea1)
    fprintf('Ea2 =  %0.4g \n', Ea2)
    fprintf('correlation r =  %0.4g \n', r) 
 
end

 
 
 
 