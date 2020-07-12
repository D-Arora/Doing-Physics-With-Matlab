% m41_simulations.m


%%  m41EV.docx       electric field due to two charges at two points
clear all
close all
clc

N = 2;
% Charge  Q = [10, 0, 0, 0, 0] 
   Q = [20, -20] .* 1e-6;
   
% X & Y components of position of charges  [0, 0, 0, 0, 0]
   xC = [-0.5,  0.5];
   yC = [0,  0];

  x = [0 0.25];
  y = [0.5 0.5];
  
  eps0 = 8.854e-12;
  kC = 1/(4*pi*eps0);
  V = zeros(N,N);
  Rx = zeros(N,N); Ry = zeros(N,N); R = zeros(N,N); R3 = zeros(N,N);
  V = zeros(N,N); Ex = zeros(N,N); Ey = zeros(N,N); 
  theta = zeros(N,N);


% CALCULATION: POTENTIAL & ELECTRIC FIELD ================================

for nC = 1 : N
 for n = 1 : N
   Rx(nC,n) = x(n) - xC(nC);
   Ry(nC,n) = y(n) - yC(nC);
   R(nC,n)  = sqrt(Rx(nC,n).^2 + Ry(nC,n).^2);
   theta(nC,n) = atan2d(Ry(nC,n),Rx(nC,n));
   V(nC,n) = kC .* Q(nC) ./ R(nC,n);
   R3(nC,n) = R(nC,n).^3;
   Ex(nC,n) = kC .* Q(nC) .* Rx(nC,n) ./ R3(nC,n);
   Ey(nC,n) = kC .* Q(nC) .* Ry(nC,n) ./ R3(nC,n);
end
end
Rx
Ry
R
theta
V
Ex
Ey
Vtot = sum(V)
Extot = sum(Ex)
Eytot = sum(Ey)
Etot = sqrt(Extot.^2 + Eytot.^2)
thetaE = atan2d(Eytot,Extot)


%% m41Voltage.docx   Example 1 solution
clear all
close all
clc
format shorte
 
q = 5.23e-6
E = 1.56e3
m = 5e-3
s = 100e-3

FE = q * E
a = FE/m
v = sqrt(2*a*s)
EKB = 0.5*m*v^2
WBA = FE*s
dEP = - WBA
dV=dEP/q


%% Coulomb's Law     m41B.docx   ========================================
clear all
close all
clc

% INPUTS
r = 250e-3;
QA = 2.23e-6;
QB = 4.45e-6;
theta = 180;

% CALCULATIONS
[F, Fx, Fy] = coulomb(QA,QB,r,theta);

% DISPLAY RESULTS
disp('INPUTS:')
textD = ['  r  =  ',num2str(r),'  m'];
disp(textD);
textD = ['  QA =  ',num2str(QA),'  C'];
disp(textD)
textD = ['  QB =  ',num2str(QB),'  C'];
disp(textD)
textD = ['  angle: theta =  ',num2str(theta),'  deg'];
disp(textD);
disp('OUTPUTS:')
textD = ['  FE =  ',num2str(F),'  N'];
disp(textD);
textD = ['  Fx =  ',num2str(Fx),'  N'];
disp(textD);
textD = ['  Fy =  ',num2str(Fy),'  N'];
disp(textD);

% FUNCTION ===============================================================

function [F, Fx, Fy] = coulomb(QA,QB,r, theta)
  k = 9e9;
  F = k * QA*QB / r^2;

  Fx = F * cosd(theta);
  Fy = F * sind(theta);
end

% %% Coulomb's Law   Graphs    m41B.docx =================================
% clear all
% close all
% clc
% 
% n = 5;
% e = 1.602e-19;
% QA = n*e;
% QB = n*e;
% k = 9e9;
% rMin = 0.1e-9;
% rMax = 10e-9;
% NR = 500;
% r = linspace(rMin,rMax,NR);
% FE1 = (k*QA*QB) ./ r.^2;
% FE2 = 5.*FE1;
% 
% figure(1)
% fs = 14;
% set(gcf,'Units','Normalized');
% set(gcf,'Position',[0.2 0.2 0.2 0.25]);
% set(gcf,'color',[1 1 1]);
% xP = r.*1e9; yP = FE1.*1e9;
% plot(xP,yP,'lineWidth',2);
% hold on
% xP = r.*1e9; yP = FE2.*1e9;
% plot(xP,yP,'lineWidth',2);
% xlabel('r  [ nm ]','fontsize',fs);
% ylabel('F_E  [ nN ]','fontsize',fs);
% set(gca,'xTick',0:2:10);
% set(gca,'yTick',0:5:20);
% set(gca,'yLim',[0 20]);
% legend('n_A = 5  n_B = 5','n_A = 5  n_B = 25');
% grid on
% box on
% set(gca,'fontsize',fs);
% 
% %  F function of q
% rQ = [2 6] .* 1e-9;
% n = 0:25;
% QA = e;
% QB = n .* e;
% FQ1 = (k*QA) .* QB ./ rQ(1)^2;
% FQ2 = (k*QA) .* QB ./ rQ(2)^2;
% 
% figure(2)
% fs = 14;
% set(gcf,'Units','Normalized');
% set(gcf,'Position',[0.5 0.2 0.2 0.25]);
% set(gcf,'color',[1 1 1]);
% xP = QB; yP = FQ1;
% plot(xP,yP,'bo');
% hold on
% xP = QB; yP = FQ2;
% plot(xP,yP,'ro');
% xlabel('Q_B  [ C ]','fontsize',fs);
% ylabel('F_E  [ N ]','fontsize',fs);
% %set(gca,'xTick',0::10);
% %set(gca,'yTick',0:5:20);
% set(gca,'xLim',[0 4.1e-18]);
% legend('r_1 = 2  nm','r_2 = 6  nm','location','northwest');
% grid on
% box on
% set(gca,'fontsize',fs);
