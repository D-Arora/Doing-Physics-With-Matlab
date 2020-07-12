%close all
%clear all
%clc

xS = linspace(0,26e-3,500);

aS1 = a(1);
aS2 = a(2);
aS3 = a(3);
aS4 = -a(4);
%yS = aS1  + aS2.* exp(-aS3.*xS);
%n = length(yS);
% aS1 = -1.1544;
% aS2 = 1550.2;
% aS3 = 0.0003205;
% aS4 = 0.0036732;
 yS = aS1  + aS2.* exp(-((xS-aS3)./aS4)); %+ rand(length(10),1);
% 
% dxS = 0.0005.* xS;
%5dyS = 0.0005 .* yS;
 
 
% n = 40;
% y(n) = 6;
% y(1) = 2055;
% y(n2) = 45;
% a(1,1) = y(n);
%     a(2,1) = y(1);
%     n2 = round(n/2);
%    a(4,1) = (x(n2) - x(1))/log((y(n2)-a(1,1))/(y(1)-a(1,1)));
%     a(3,1) = log((y(1)-a(1,1))/(y(n2)-a(1,1)))/(x(n2)-x(1));
% a(1) = 6; a(2) = 2000; a(3) = 200;
% eqType = 10;
% yS = fitFunction(a, xS, eqType) ;   
figure(88)
plot(xS,yS);
grid on
hold on
plot(xf,yf,'r');
%  wData(:,1) = xS;
%  wData(:,2) = yS;
%  wData(:,3) = dxS;
%  wData(:,4) = dyS;
%  
% weighted2


