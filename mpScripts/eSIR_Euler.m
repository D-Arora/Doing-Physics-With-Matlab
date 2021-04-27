% eSIR_Euler

clear 
clc
close all

 a = 0.27;
 b = 0.027;
 f = 18e5;
 
 % Initialize matrices
  n = 50000;
  S  = zeros(n,1);       % Susceptible population
  I = zeros(n,1);        % Active infected population
  R = zeros(n,1);        % Removals form infected population  
  
  tMax = 1000;
  t = linspace(0,tMax,n);
  h = t(2) - t(1);
 
  S(1) = 1;
  I(1) = 2e-9;
  R(1) = 0;
  
 % Euler
 for c = 1 : n-1
     S(c+1) = S(c) - h*a*S(c)*I(c);
     I(c+1) = I(c) + h*I(c)*(a*S(c) - b);
     R(c+1) = R(c) + h*b*I(c);
     
     if c > 800  && c < 2601;  S(c+1) = 0.2; end
     if c > 2600  && c < 3201;  S(c+1) = 0.3; end
 end
 
 % Graphics
  figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.02 0.01 0.4 0.8]);
   set(gcf,'color','w');
   FS = 14;
   xP = t;
   
 subplot(4,1,1)
   yP = S;
   plot(xP,yP,'b','linewidth',2)
   grid on
   set(gca,'fontsize',FS)
   ylabel('S');
 
 subplot(4,1,2)
   yP = I+R;
   plot(xP,yP,'b','linewidth',2)
   grid on
   set(gca,'fontsize',FS)
   ylabel('I_{tot}'); 
   
 subplot(4,1,3)
   yP = I;
   plot(xP,yP,'b','linewidth',2)
   grid on
   set(gca,'fontsize',FS)
   ylabel('I');  
   
subplot(4,1,4)
   yP = R;
   plot(xP,yP,'b','linewidth',2)
   grid on
   set(gca,'fontsize',FS)
   ylabel('R');     
   xlabel('t')
   