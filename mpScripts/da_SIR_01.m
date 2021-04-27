%da_SIR_01.m

clear
close all
clc

global a b
% Model Parameters
  a = 0.2;
  b = 0.05;
  num = 5000;
  tMax = 300;
  
% Initialize matrices
  S  = zeros(num,1);       % Susceptible population
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population  
  C = zeros(num,1);        % ReCovered population
  D = zeros(num,1);        % Dead population
  
  S(1) = 1;
  I(1) = 0.001;
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
   kS1 = SDOT(t(n), S(n), I(n), R(n));
   kI1 = IDOT(t(n), S(n), I(n), R(n));
   kR1 = RDOT(t(n), S(n), I(n), R(n));
   
   kS2 = SDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kI2 = IDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kR2 = RDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   
   kS3 = SDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kI3 = IDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kR3 = RDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   
   kS4 = SDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kI4 = IDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kR4 = RDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   
   S(n+1) = S(n) + h*(kS1+2*kS2+2*kS3+kS4)/6;
   I(n+1) = I(n) + h*(kI1+2*kI2+2*kI3+kI4)/6;
   R(n+1) = R(n) + h*(kR1+2*kR2+2*kR3+kR4)/6;
   
 %  if n == 2500; S(n+1) = 0.8; end
end
  
dIdt = a.*S.*I - b.*I;

Iinc = I(1).*exp((a-b).*t);

Iindex = find(t>50,1);
I0 = 0.82;
tdec = t(Iindex:end);
Idec = I0.*exp(-b.*(tdec - tdec(1)));


% GRAPHICS  ===========================================================


figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.10 0.10 0.22 0.32]);
   set(gcf,'color','w');
   FS = 10;
   
   
subplot(2,1,1)
   yyaxis left
   plot(t,S,'b','linewidth',2)
   hold on
   plot(t,R,'r-','linewidth',2)
   grid on; box on;
   xlabel('Days elapsed','fontname','times')
 %  title('S and R populations','Fontweight','normal')
   ylabel('S, R','fontname','times')
   ylim([0 1.1])
   set(gca,'ytick',0:0.2:1)
   text(180,0.87,'R','color','r','fontname','times','fontsize',FS)
   text(12,0.87,'S','color','b','fontname','times','fontsize',FS)
   text(38,0.9,'R_e','color','k','fontname','times','fontsize',FS)
   ytickformat('%,0.1f')
   set(gca,'YColor',[0 0 1]);
   set(gca,'fontsize',FS)
   
   yyaxis right
   plot(t,S.*a./b,'k','linewidth',2)
   set(gca,'YColor',[0 0 0]);
   ylabel('R_e','fontname','times')
   set(gca,'fontsize',FS)
   set(gca,'fontname','times')
   
   
 subplot(2,1,2)
   yyaxis left
   plot(t,I,'b','linewidth',2)
   hold on
   plot(t,Iinc,'k-')
   plot(tdec,Idec,'k-')
   
   grid on; box on;
   xlabel('Days elapsed','fontname','times')
 %  title('Active infections','Fontweight','normal')
   ylabel('Active infections','fontname','times')
   ylim([0 1.1])
   set(gca,'ytick',0:0.2:1)
   set(gca,'YColor',[0 0 1]);
   ytickformat('%,0.1f')
   text(6,0.2,'e^{(a-b) t}','fontname','times','fontsize',12)
   text(60,0.7,'e^{- b t}','fontname','times','fontsize',12)
   text(170,0.5,'dI/dt','color','r','fontname','times','fontsize',12)
   text(170,0.1,'I','color','b','fontname','times','fontsize',12)
   set(gca,'fontsize',FS)  
   
   yyaxis right
   plot(t,dIdt,'r','linewidth',2)
   ylabel('dI/dt','fontname','times')
   set(gca,'fontsize',12) 
   set(gca,'YColor',[1 0 0]);
   set(gca,'fontsize',FS)
   set(gca,'fontname','times')
% subplot(3,1,3)   
%    plot(t,R,'b','linewidth',2)
%    grid on; box on;
%    title('Removed Population','Fontweight','normal')
%    xlabel('days elapsed')
%    ylabel('R','fontname','times')
%    set(gca,'fontsize',FS)
%    ylim([0 1.1])
%    set(gca,'ytick',0:0.2:1)
%    set(gca,'fontsize',FS)
   
figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.50 0.10 0.25 0.25]);
   set(gcf,'color','w');
   
   plot(t,S,'b','linewidth',2)
   hold on
   plot(t,I,'r','linewidth',2)
   plot(t,R,'k','linewidth',2)
   plot(t,I+R,'m','linewidth',2)
   grid on; box on;
   xlabel('t')
   ylabel('S  I  R')
   legend('S','I','R','I_{tot}','location','east')
   set(gca,'fontsize',12)
   ylim([0 1.1])
   set(gca,'ytick',0:0.2:1)
   
      

% FUNCTIONS  ===========================================================

function FS = SDOT(t,S,I,R)
  global a
    FS = -a*S*I;
end


function FI = IDOT(t,S,I,R)
  global a b
    FI = a*S*I - b*I;
end

function FR = RDOT(t,S,I,R)
  global b
    FR = b*I;
end

