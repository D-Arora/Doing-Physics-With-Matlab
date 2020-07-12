%ztest1.m

%clear
%close all
%clc
k = 1.6e-5;
t = linspace(0,14e4,501);
s = 3.5e4.*(1-exp(-k.*t));



%zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
%zx = zci(s)    
hold on

figure(99)
%  set(gcf,'units','normalized');
%   set(gcf,'position',[0.02 0.02 0.2 0.2]);
%   set(gcf,'color','w');
% plot([-3 10],[0 0],'k','linewidth',2)
% hold on

 plot(t,s,'b','linewidth',1)

 
% t = linspace(-3,10,501);
% s = 4.*exp(-(t-3).^2./4);
% hold on
% plot(t,s,'r','linewidth',3)
%  axis off
