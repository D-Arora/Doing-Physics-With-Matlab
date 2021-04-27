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
  g = 1e-3;
% omega
  w = 1e-2;
% beta 
  b = 0.06;
% delta
  d = 0.03;

% Time
  N = 999; tMax = 20;
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
  M3 = zeros(N,1); M3(1) = 1e5;  
% Amyloid-beta peptide: number of molecules  
  Ab = zeros(N,1); Ab(1) = 1e3;
  

% Compute y and yDdot values as a function of time =====================
for n = 1 : N-1
   k1Ns = Ns_DOT(t(n),Ns(n),Aq(n),Ap(n),M2(n),M1(n),Ab(n),M3(n));
   k1Nd = -k1Ns;
   k1Aq = Aq_DOT(t(n),Ns(n),Aq(n),Ap(n),M2(n),M1(n),Ab(n),M3(n));
   k1Ap = -k1Aq;
   k1M2 = M2_DOT(t(n),Ns(n),Aq(n),Ap(n),M2(n),M1(n),Ab(n),M3(n),Nd(n));
   k1M1 = -k1M2;
   k1Ab = Ab_DOT(t(n),Ns(n),Aq(n),Ap(n),M2(n),M1(n),Ab(n),M3(n));
   k1M3 = M3_DOT(t(n),Ns(n),Aq(n),Ap(n),M2(n),M1(n),Ab(n),M3(n));
   
   
   k2Ns = Ns_DOT(t(n)+h/2,Ns(n)+k1Ns/2,Aq(n)+k1Aq/2,Ap(n)+k1Ap/2 ...
          ,M2(n)+k1M2/2,M1(n)+k1M1/2,Ab(n)+k1Ab/2,M3(n)+k1M3/2);
   k2Nd = -k2Ns;   
      
   k2Aq = Aq_DOT(t(n)+h/2,Ns(n)+k1Ns/2,Aq(n)+k1Aq/2,Ap(n)+k1Ap/2 ...
          ,M2(n)+k1M2/2,M1(n)+k1M1/2,Ab(n)+k1Ab/2,M3(n)+k1M3/2);
   k2Ap = -k2Aq;    
           
   k2M2 = M2_DOT(t(n)+h/2,Ns(n)+k1Ns/2,Aq(n)+k1Aq/2,Ap(n)+k1Ap/2 ...
          ,M2(n)+k1M2/2,M1(n)+k1M1/2,Ab(n)+k1Ab/2,M3(n)+k1M3/2,Nd(n)+k1Nd/2);
   k2M1 = -k2M2;   
       
   k2Ab = Ab_DOT(t(n)+h/2,Ns(n)+k1Ns/2,Aq(n)+k1Aq/2,Ap(n)+k1Ap/2 ...
          ,M2(n)+k1M2/2,M1(n)+k1M1/2,Ab(n)+k1Ab/2,M3(n)+k1M3/2);
       
   k2M3 = M3_DOT(t(n)+h/2,Ns(n)+k1Ns/2,Aq(n)+k1Aq/2,Ap(n)+k1Ap/2 ...
          ,M2(n)+k1M2/2,M1(n)+k1M1/2,Ab(n)+k1Ab/2,M3(n)+k1M3/2); 
      
      
   k3Ns = Ns_DOT(t(n)+h/2,Ns(n)+k2Ns/2,Aq(n)+k2Aq/2,Ap(n)+k2Ap/2 ...
          ,M2(n)+k2M2/2,M1(n)+k2M1/2,Ab(n)+k2Ab/2,M3(n)+k2M3/2);
   k3Nd = -k3Ns;   
      
   k3Aq = Aq_DOT(t(n)+h/2,Ns(n)+k2Ns/2,Aq(n)+k2Aq/2,Ap(n)+k2Ap/2 ...
          ,M2(n)+k2M2/2,M1(n)+k2M1/2,Ab(n)+k2Ab/2,M3(n)+k2M3/2);
   k3Ap = -k3Aq;   
      
   k3M2 = M2_DOT(t(n)+h/2,Ns(n)+k2Ns/2,Aq(n)+k2Aq/2,Ap(n)+k2Ap/2 ...
          ,M2(n)+k2M2/2,M1(n)+k2M1/2,Ab(n)+k2Ab/2,M3(n)+k2M3/2,Nd(n)+k2Nd/2 );
   k3M1 = -k3M2;   
        
   k3Ab = Ab_DOT(t(n)+h/2,Ns(n)+k2Ns/2,Aq(n)+k2Aq/2,Ap(n)+k2Ap/2 ...
          ,M2(n)+k2M2/2,M1(n)+k2M1/2,Ab(n)+k2Ab/2,M3(n)+k2M3/2);
   
   k3M3 = M3_DOT(t(n)+h/2,Ns(n)+k2Ns/2,Aq(n)+k2Aq/2,Ap(n)+k1Ap/2 ...
          ,M2(n)+k2M2/2,M1(n)+k2M1/2,Ab(n)+k2Ab/2,M3(n)+k1M3/2);   
      
  
   k4Ns = Ns_DOT(t(n)+h,Ns(n)+k3Ns,Aq(n)+k3Aq,Ap(n)+k3Ap ...
          ,M2(n)+k3M2,M1(n)+k3M1,Ab(n)+k3Ab,M3(n)+k3M3);
   k4Nd = -k4Ns;
      
   k4Aq = Aq_DOT(t(n)+h,Ns(n)+k3Ns,Aq(n)+k3Aq,Ap(n)+k3Ap ...
          ,M2(n)+k3M2,M1(n)+k3M1,Ab(n)+k3Ab,M3(n)+k3M3);
   k4Ap = -k4Aq;
      
   k4M2 = M2_DOT(t(n)+h,Ns(n)+k3Ns,Aq(n)+k3Aq,Ap(n)+k3Ap ...
          ,M2(n)+k3M2,M1(n)+k3M1,Ab(n)+k3Ab,M3(n)+k3M3,Nd(n)+k3Nd);
   k4M1 = -k4M2;   
           
   k4Ab = Ab_DOT(t(n)+h,Ns(n)+k3Ns,Aq(n)+k3Aq,Ap(n)+k3Ap ...
          ,M2(n)+k3M2,M1(n)+k3M1,Ab(n)+k3Ab,M3(n)+k3M3);
      
   k4M3 = M3_DOT(t(n)+h,Ns(n)+k3Ns,Aq(n)+k3Aq,Ap(n)+k3Ap ...
          ,M2(n)+k3M2,M1(n)+k3M1,Ab(n)+k3Ab,M3(n)+k3M3);  
      
      
   Ns(n+1) = Ns(n) + h*(k1Ns + 2*k2Ns + 2*k3Ns + k4Ns)/6;
   Aq(n+1) = Aq(n) + h*(k1Aq + 2*k2Aq + 2*k3Aq + k4Aq)/6;
   M2(n+1) = M2(n) + h*(k1M2 + 2*k2M2 + 2*k3M2 + k4M2)/6;
   Ab(n+1) = Ab(n) + h*(k1Ab + 2*k2Ab + 2*k3Ab + k4Ab)/6;
   M3(n+1) = M3(n) + h*(k1M3 + 2*k2M3 + 2*k3M3 + k4M3)/6;
   
   Nd(n+1) = Nd(n) - (Ns(n+1) - Ns(n));
   Ap(n+1) = Ap(n) - (Aq(n+1) - Aq(n));
   M1(n+1) = M1(1) - (M2(n+1) - M2(n));
    
    
end  


% GRAPHICS  ===========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.02 0.10 0.3 0.3]);
   set(gcf,'color','w');
   
   FS = 14;
   
   xP = t; yP = Ns;
   plot(xP,yP,'b','linewidth',2)
   hold on
   yP = M1;
   plot(xP,yP,'r','linewidth',2)
   yP = Ab;
   plot(xP,yP,'g','linewidth',2)
   grid on
   box on
  
   xlabel('time  [years]')
   ylabel('Populations N_s, M_1, A\beta') 
   legend('Ns', 'M1', 'A\beta')
   set(gca,'fontsize',FS,'fontname','times')
  
figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.32 0.10 0.3 0.3]);
   set(gcf,'color','w');
   
   FS = 14;
   
   xP = t; yP = Nd;
   plot(xP,yP,'b','linewidth',2)
   hold on
   yP = Ap;
   plot(xP,yP,'r','linewidth',2)
   yP = Ab;
 %  plot(xP,yP,'g','linewidth',2)
   grid on
   box on
  
   xlabel('time  [years]')
   ylabel('Populations N_s, M_1, A\beta') 
   legend('Nd', 'Ap')
   set(gca,'fontsize',FS,'fontname','times')   
   

toc   
   
% FUNCTIONS  ==========================================================

function f = Ns_DOT(t,Ns,Aq,Ap,M2,M1,Ab,M3)
   global a g w b d
   f = a(1)*Aq - a(2)*Ap - a(3)*M1;
end

function f = Aq_DOT(t,Ns,Aq,Ap,M2,M1,Ab,M3)
   global a g w b d
   f = a(4)*M2 - a(5)*M1;
end

function f = M2_DOT(t,Ns,Aq,Ap,M2,M1,Ab,M3,Nd)
   global a g w b d
   f = (a(6)+a(11))*Ns - a(10)*Nd + (a(7)+a(12))*Aq - a(9)*M1 + a(14)*M2 - (a(8)+a(13))*Ab; 
end

function f = Ab_DOT(t,Ns,Aq,Ap,M2,M1,Ab,M3)
   global a g w b d
   f = -a(15)*Ns - a(16)*M2 + a(17)*Ab + g*M3;
end

function f = M3_DOT(t,Ns,Aq,Ap,M2,M1,Ab,M3)
   global a g w b d
   f = w - b*M1 + d*M3;
end

