% QMG23C.m
% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG23C.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230519   Matlab R2021b

% Schrodinger equation: Harmonic Oscillator
% Hermite Polynomials solutions for wavefunction

clear; close all;clc

% INPUTS =========================================================
   NX = 801;                   % x grid points
   xMin = -0.25; xMax = 0.25;  % x range [nm]
% Stationary state quantum number n   
   n = 1;                   
% Potential well   depth [eV]  /  width  [nm]
   U0 = -400; x1 = 0.5;

   NE = 13;    % number of energy levels displayed in Figure(4)

% CONSTANTS ========================================================
  h    = 6.62607015e-34;
  hbar = 1.05457182e-34;
  me    = 9.1093837e-31;
  e    = 1.60217663e-19;
  Ls = 1e-9; Es = e;      % scaling factors: position [nm]  energy [eV]
 % Cs = -hbar^2/(2*me*Ls^2*Es);

% SETUP ===========================================================   
   x = linspace(xMin,xMax,NX);
   k = -8*(U0*Es)/(x1*Ls)^2;      % effective spring constant [N/m]
   w = sqrt(k/me);                % fundamental frequency [rad/s]

   xi = sqrt(me*w/hbar).*(x*Ls);       % scaled x 
   C1 = (me*w/(pi*hbar))^0.25;
   
% Wavefunction and probabilty density   
   psi = zeros(NE,NX); E = zeros(NE,1); probD = zeros(NE,NX);
   for c = 1:NE
     C2 = sqrt(1/(2^(c-1)*factorial(c-1)));
     psi(c,:) = (C1*C2).* hermiteH(c-1,xi).*exp(-xi.^2/2);
     E(c) = hbar*w*(c-1/2)/e;
   end

% Probability density
  for c = 1:NE
    probD(c,:) = psi(c,:).*psi(c,:);
  end

% Check normalized wavefunction
  check = zeros(NE,1);
  for c = 1:NE
    fn = probD(c,:); %psi(c,:).*psi(c,:);
    check(c) = Ls*simpson1d(fn,x(1),x(NX));
  end



% Potential energy / total energy / kinetic energy
     U = (0.5*k).*(x.*Ls).^2;
     Ue = U./Es;
     M = 4;
     EM = E(4).*ones(1,NX);
     K = EM - Ue;
% K > 0
     Kind = find(K>0); K1 = Kind(1); K2 = Kind(end);
     
% Prob of electron in forbidden region for stationary state n = 2 (c = 3);   
  fn = probD(3,1:K1-1);
  xS = x(1); xF = x(K1-1); 
  probF = 2*Ls*simpson1d(fn,xS,xF);

% Hermite polynomials ============================================= 
   NP = 5;
   Hn = zeros(NX,NP);
   xH = linspace(-2,2,NX);
   for c = 1:NP
     Hn(:,c) = hermiteH(c-1,xH);
   end


% GRAPHICS  =======================================================
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.05 0.30 0.30]);
   set(gcf,'color','w');
   FS = 14;
hold on
   xP = xH;
   for c = 1:NP
     yP = Hn(:,c);
     plot(xP,yP,'LineWidth',2);
   end
   grid on
   ylim([-40 40])
   xlabel('x  [ nm ]','FontSize',14);
   ylabel('H_n(x)','FontSize',14);
   title('Hermite Poylnomials','FontWeight','normal')
   legend('H_0(x)', 'H_1(x)', 'H_2(x)', 'H_3(x)', 'H_4(x)', ...
       'Location', 'north','Orientation','horizontal','box','off')
   box on
   set(gca,'fontsize',FS)

figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.4 0.05 0.30 0.60]);
   set(gcf,'color','w');
   FS = 14;

   xP = x;
   for c = 1:6
       subplot(3,2,c)   
       yP = psi(c,:);
       plot(xP,yP,'b','LineWidth',2);
       xticks(-0.25:0.25:0.25)
       grid on
       xlabel('x  [ nm ]')
       ylabel('\psi(x)');
       txt = sprintf('n = %0.0f  \n',c-1);
       title(txt,'FontWeight','normal')
       set(gca,'fontsize',FS)
   end

figure(3)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.6 0.25 0.25]);
   set(gcf,'color','w');
   FS = 14;

   xP = 0:NE-1; yP = E(1:NE);
       Hplot = plot(xP,yP,'bo');
       set(Hplot,'markersize',6,'MarkerFaceColor','b')
       grid on
       xlabel('n')
       ylabel('E_n','FontSize',14);
       txt = sprintf('E_0 = %0.4f eV  \n',E(1));
       text(3,12,txt,'FontSize',14)
       txt = sprintf('\\DeltaE = %0.4f eV  \n',E(2)-E(1));
       text(8,210,txt,'FontSize',14)
       set(gca,'fontsize',FS)
  
 figure(4)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.65 0.6 0.25 0.25]);
   set(gcf,'color','w');
   FS = 14;

   xP = x; yP = Ue;
       plot(xP,yP,'b','LineWidth',2);
       hold on
       for c = 1:NE
           xP = [xMin xMax]; yP = [E(c), E(c)];
           plot(xP,yP,'r','linewidth',2)
       end
    %   xlim([-0.2, 0.2])
       grid on
       xlabel('x')
       ylabel('U & E_n  [ eV ]','FontSize',14);
       set(gca,'fontsize',FS)

figure(5)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.7 0.05 0.25 0.50]);
   set(gcf,'color','w');
   FS = 14;
  
   xP = x;
subplot(2,1,1)   
   yP = Ue;
   plot(xP,yP,'b','LineWidth',2);
   hold on
   plot(xP,EM,'r','linewidth',2)
   plot(xP,K,'k','linewidth',2)
   plot([xMin xMax], [0 0],'m','linewidth',1)
   plot([x(K1) x(K1)],[-400 400],'b')
   plot([x(K2) x(K2)],[-400 400],'b')
   grid on
   xlabel('x')
   ylabel('U, E_n, K  [ eV ]','FontSize',14);
   set(gca,'fontsize',FS)
subplot(2,1,2)   
   yP = psi(3,:);
   plot(xP,yP,'b','LineWidth',2);
   hold on
   plot([x(K1) x(K1)],[-2e5 2e5],'b')
   plot([x(K2) x(K2)],[-2e5 2e5],'b')
   grid on
   xlabel('x')
   ylabel('\psi_2','FontSize',14);
   set(gca,'fontsize',FS) 

 figure(6)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.55 0.6 0.25 0.25]);
   set(gcf,'color','w');
   FS = 14;

   xP = x; yP = probD(3,:);
   plot(xP,yP,'b','LineWidth',2);
   hold on
   xP = x(1:K1); yP = probD(3,1:K1);
   area(xP,yP,'FaceColor','r')
   xP = x(K2:NX);yP = probD(3,K2:NX);
   area(xP,yP,'FaceColor','r')
   grid on; box on
   xlabel('x  [ nm ]')
   ylabel('\psi(x)','FontSize',14);
   txt = sprintf('Prob. electron in fobidden region = %0.4f  \n',probF);
   title(txt,'FontWeight','normal')
   set(gca,'fontsize',FS)

  