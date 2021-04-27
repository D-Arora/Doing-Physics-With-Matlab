% bp001.m

clear
close all
clc



% INPUTS  =============================================================  
R0 = 2.1;
EG0 = 1.0e-3;
SI = 3.06e-3;
a = 1e4;
Imax = 0.28;
kI = 0.01;
kSTO = 0.04;
kGUT = 0.12;

G_basal = 70.86;

I_basal = 9.36; 



NT = 3600*24;


% SETUP  ==============================================================
tMax = 24*60;
t = linspace(0,tMax,NT);
dt = t(2) - t(1);

gGUT = zeros(NT,1);
G    = zeros(NT,1);
I    = zeros(NT,1);




% FOOD INTAKE  =======================================================

% Stomach
    Q1 = 200;
    k1 = 0.01;
    g1 = zeros(NT,1);
    R1 = 8*3600:NT;
    g1(R1) = Q1.*exp(-k1 .* (t(R1)-t(R1(1))));

    Q2 = 0;
    k2 = 0.036;
    g2 = zeros(NT,1);
    R2 = 10*3600:NT;
    g2(R2) = Q2.*exp(-k2 .* (t(R2)-t(R2(1))));

    
   gSTO = k1*g1+ k2.*g2;
    
    

%

% COMPUTATIONS  =======================================================
 
G(1) = G_basal;
I(1) = I_basal;



for n = 1 : NT-1

% GLUCOSE -------------------------------------------------------------
% Glucose Intake: GUT
 %   gGUT(n+1) = gGUT(n) + dt * (k1*g1(n) - kGUT*gGUT(n));
  gGUT(n+1) = gGUT(n) + dt * (k1*g1(n)+k2*g2(n) - kGUT*gGUT(n));  
    
 %   gGUT(n+1) = gGUT(n) + dt * (gSTO(n) - kGUT*gGUT(n));
    G_gut = dt*kGUT*gGUT(n+1);

% Glucose production
   G_prod = dt*R0;

% Glucose decay   
   G_decay = -dt*EG0*G(n);
   
% Glucose metabolism by insulin
  G_insulin =  -dt*SI*I(n)*G(n);
 
    

% Glucose total  
  G(n+1)    = G(n) + G_gut + G_prod + G_decay + G_insulin;

% INSULIN  ------------------------------------------------------------    
% Insulin  pancrease
  I_pancrease = dt * Imax*(G(n)^2/(a + G(n)^2)) ;

% I_decay
  I_decay =  - dt*kI*I(n);
    
  I(n+1)    = I(n) + I_pancrease + I_decay;  
end





% GRAPHICS  ===========================================================

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.0 0.10 0.35 0.7]);
  set(gcf,'color','w');
  FS = 14;
  
subplot(4,1,1)
xP = t./60; yP = gSTO;
plot(xP,yP,'linewidth',2)
grid on
xlim([0 24])
set(gca,'xtick',0:4:24)
%set(gca,'ytick',0:20:100)
ylabel('g_{STO}')
set(gca,'fontsize',FS)


subplot(4,1,2)
xP = t./60; yP = gGUT;
plot(xP,yP,'linewidth',2)
grid on
xlim([0 24])
set(gca,'xtick',0:4:24)
ylabel('g_{GUT}')
set(gca,'fontsize',FS)


subplot(4,1,3)
xP = t./60; yP = G;
plot(xP,yP,'linewidth',2)
grid on
xlim([0 24])
set(gca,'xtick',0:4:24)
%set(gca,'ytick',0:10:max(y))
ylabel('glucose')
set(gca,'fontsize',FS)


subplot(4,1,4)
xP = t./60; yP = I;
plot(xP,yP,'linewidth',2)
grid on
xlim([0 24])
%ylim([0 20])
set(gca,'xtick',0:4:24)
%ylabel('insulin')
xlabel('time  [h]')
set(gca,'fontsize',FS)



% Sigmoidal function

  f = linspace(0,200,2000);
  S = zeros(2000,1);
  for n = 1 : 2000
        S(n) = sigm(f(n),0,1e3);
  end   

figure(2)
   set(gcf,'units','normalized');
  set(gcf,'position',[0.5 0.10 0.25 0.35]);
  set(gcf,'color','w');
  FS = 14;
  
  plot(f,S)
  

% FUNCTIONS  ==========================================================

function S = sigm(f,f1,f2)
 
  F = (f-f1)^2;
  S = F / (f2 + F);

end

