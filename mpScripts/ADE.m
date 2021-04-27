% ADA.m

% 10 07 2020 Matlab 2019b
% AD  Alzheimer's Dissease  


clear
close all
clc

global a g w b d

tic

% CONSTANTS  ==========================================================

% alpha
  a(1)  = 1e-5;  a(2)  = 1e-3;  a(3)  = 1e-2;  a(4)  = 1e-4;
  a(5)  = 1e-2;  a(6)  = 1e-2;  a(7)  = 1e-4;  a(8)  = 1e-2;
  a(9)  = 1e-2;  a(10) = 1e-2;  a(11) = 1e-2;  a(12) = 1e-4;
  a(13) = 1e-2;  a(14) = 1e-4;  a(15) = 1;     a(16) = 1e-2;
  a(17) = 1;   % alpha_r

% gamma
  g = 1e-2;   % 1e-3;
 % omega
  w = 0.01 ;%1e-2;
% beta 
  b = 0.06;
% delta
  d = 0.03;

% Time
  N = 9999; tMax = 30;
  t = linspace(0,tMax,N);
  h = t(2) - t(1);
  
% VARIABLES and INITIAL VALUES  =======================================

% neuron pop: survived
  Ns = zeros(N,1); Ns(1) = 1e4;
% neuron pop: dead
  Nd = zeros(N,1); Nd(1) = 1e2;  
% Astroglia: quiescenent
  Aq = zeros(N,1); Aq(1) = 1e5;
% Astroglia: proliferation
  Ap = zeros(N,1); Ap(1) = 1e3;
% Microglia: normal - resting / antiflammatory
  M2 = zeros(N,1); M2(1) = 1e5;
% Microglia: reactive - active pro-inflamatory
  M1 = zeros(N,1); M1(1) = 1e3;  
% Microglia: neutral
  M3 = zeros(N,1); M3(1) = 1e3;  
% Amyloid-beta peptide: number of molecules  
  Ab = zeros(N,1); Ab(1) = 1e3;
  

  for n = 1: N-1
  
   Ns(n+1) = Ns(n) + h*(a(1)*Aq(n) - a(2)*Ap(n) - a(3)*M1(n));
 
   Nd(n+1) = Nd(n) + h*(- a(1)*Aq(n) + a(2)*Ap(n) + a(3)*M1(n));
  
   Aq(n+1) = Aq(n) + h*(a(4)*M2(n) - a(5)*M1(n));
   
   Ap(n+1) = Ap(n) - h*(a(4)*M2(n) + a(5)*M1(n));
   
   M2(n+1) = M2(n) + h*((a(6)+a(11))*Ns(n) - a(10)*Nd(n) + (a(7)+a(12))*Aq(n) ...
             - a(9)*M1(n) + a(14)*M2(n) - (a(8)+a(13))*Ab(n)); 

   M1(n+1) = M1(n) +h*(-(a(6)+a(11))*Ns(n)+ a(10)*Nd(n) - (a(7)+a(12))*Aq(n) ...
             + a(9)*M1(n) - a(14)*M2(n) + (a(8)+a(13))*Ab(n)); 
         
   Ab(n+1) = Ab(n) +h*(a(15)*Ns(n) - a(16)*M2(n) - a(17)*Ab(n) + g*M3(n));
   
%   M3(n+1) = M3(n) + h*(w - b*M1(n) + d*M3(n));
   
   %  0.00001.*x(3)-0.001.*x(4)-0.01.*x(6)
% -0.00001.*x(3)+0.001.*x(4)+0.01.*x(6)
%  0.0001.*x(5)-0.01.*x(6)
% -0.0001.*x(5)+0.01.*x(6)
%  (0.01+0.01).*x(1)-0.01.*x(2)+(0.0001+0.0001).*x(3)-0.01.*x(6)+0.0001.*x(5)-(0.01+0.01).*x(7)
% -(0.01+0.01).*x(1)+0.01.*x(2)-(0.0001+0.0001).*x(3)+0.01.*x(6)-0.0001.*x(5)+(0.01+0.01).*x(7)
% 1.*x(1)-0.01.*x(5)-1.*x(7)+0.01.*x(8)
% 0.01-0.06.*x(6)+0.03.*x(8)]);
   
  end


% GRAPHICS  ===========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.02 0.10 0.5 0.6]);
   set(gcf,'color','w');
   
   FS = 14;
   xP = t;
   subplot(2,2,1)
   
   yyaxis left
   yP = Nd/(Ns(1));
   plot(xP,yP,'r','linewidth',2)
   ylabel('N_d')
   set(gca,'YColor','r')
   
   yyaxis right
   yP = Ns/Ns(1);
   plot(xP,yP,'b','linewidth',2)
   ylabel('N_d')
   xlabel('time  [years]')
   ylabel('N_s') 
   set(gca,'fontsize',FS,'fontname','times')
   set(gca,'YColor','b');
   box on

   figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.30 0.10 0.25 0.3]);
   set(gcf,'color','w');
   
   FS = 14;
   
   xP = t; yP = Nd;
   plot(xP,yP,'b','linewidth',2)
   hold on
   yP = Ap;
   plot(xP,yP,'r','linewidth',2)
   
   grid on
   box on
  
   xlabel('time  [years]')
   ylabel('Populations  N_d,  M_1,  A_p') 
   legend('N_d', 'A_p','location','east')
   set(gca,'fontsize',FS,'fontname','times')   
   

toc   
   
