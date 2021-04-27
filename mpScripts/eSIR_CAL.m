%eSIR_CAL.m

clc
close all
clear

% save covid.mat covid

% DATA: Id Rd Dd / Number of data days for covid19 
%   load('covid.mat')
%   Ndays  = length(covid); %find(covid19 == 0,1) - 1; 
%   
%   Ndays = 280;
%   
%   Idtot  = covid(1:Ndays,1);
%   Id     = covid(1:Ndays,2);
%   Dd     = covid(1:Ndays,3);  

load('covid.mat')
  Ndays  = length(covid); %find(covid19 == 0,1) - 1;  
  Idtot  = covid(1:Ndays,4);
  Id     = covid(1:Ndays,5);
  Dd     = covid(1:Ndays,6);

% load('covid19.mat')
%   Ndays  = length(covid19); %find(covid19 == 0,1) - 1;  
%   Idtot  = covid19(1:Ndays,1);
%   Id     = covid19(1:Ndays,2);
%   Dd     = covid19(1:Ndays,3);

  
  Cd = Idtot - Id - Dd;
  Rd = Cd + Dd;
  
  
  t = 1:Ndays;
  
  b = 0.017;
   R(1) = Rd(1);
  for c = 2 : Ndays
   R(c) = R(c-1) + b*Id(c);   
  end
  
 
  errR = sqrt(sum((R'-Rd).^2))
 
  m = 80;
  N = 130;
  x1 = t(N);
  y1 = Rd(N);
  
  b1 = y1 - m*x1;
  
  y = m*t + b1;
  
  
  b = m/Id(N)
%   


Re = 1 + gradient(Id,1)./(b.*Id);

a = -gradient(Re,1)./(Re.*Id);

  
% GRAPHICS ==========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.5 0.05 0.3 0.6]);
   set(gcf,'color','w');
   FS = 14;
   
   
% subplot(2,1,1)
   xP = t;  yP = Rd;
   plot(xP,yP,'b')
   
   yP = y;
   hold on
   plot(xP,yP,'r')
   
   grid on
   xlabel('t')
   ylabel('R_d')
   set(gca,'fontsize',FS)
   
subplot(2,1,2)
   xP = t;  yP = Re;
   plot(xP,yP)
   
   hold on
   
%    yP = y;
%    plot(xP,yP,'r')
   grid on
   
   xlabel('t')
   ylabel('I_d')
   set(gca,'fontsize',FS)
   