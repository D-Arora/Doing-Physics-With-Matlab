% QMG23B.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG23F.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230517   Matlab R2021b

% ANHARMONIC OSCILLATOR: MORSE POTENTIAL

clear; close all;clc
% INPUTS  ================================================================
   NX = 801;                % x grid points
   xMin = 0; xMax = 0.5;    % x range [nm]
   x1 = 0.1;                % width of truncated potential well [nm]
   U0 = -600;               % well depth [eV]
   S = 0.12;                % Width of well
   M = 5;                   % Stationary states M 
  
% CONSTANTS ========================================================
  h    = 6.62607015e-34;
  hbar = 1.05457182e-34;
  me    = 9.1093837e-31;
  e    = 1.60217663e-19;
  Ls = 1e-9; Es = e;      % scaling factors: position [nm]  energy [eV]
  Cs = -hbar^2/(2*me*Ls^2*Es);

% SETUP  ===========================================================
  x = linspace(xMin,xMax, NX);

% Potential well 
  U_matrix = zeros(NX-2);
  U = U0.*(1-(exp((x1-x)/S)-1).^2);  

% Make potential energy matrix
  dx = (x(2)-x(1));
  dx2 = dx^2;
  for cn = 1:(NX-2)
    U_matrix(cn,cn) = U(cn+1);
  end

% KINETIC ENERGY and HAMILTONIAN MATRICES ============================
% Second Derivative Matrix 
    off = ones(NX-3,1);                 
    SD_matrix = (-2*eye(NX-2) + diag(off,1) + diag(off,-1))/dx2;
% KE Matrix
    K_matrix = Cs * SD_matrix;            
% Hamiltonian Matrix
   H_matrix = K_matrix + U_matrix;

% EIGENVALUES and EIGENFUNCTIONS =====================================
   [e_funct, e_values] = eig(H_matrix);

% All Eigenvalues 1, 2 , ... n  where E_N < 0
   flagE = 0;
   n = 1;
while flagE == 0
    E(n) = e_values(n,n);
    if E(n) > 0, flagE = 1; end 
    n = n + 1;
end  
    E(n-1) = [];
    n = n-2;

% Corresponding Eigenfunctions 1, 2, ... ,n: Normalizing the wavefunction
    psi = zeros(NX,n); area = zeros(1,n);
for cn = 1 : n
    psi(:,cn) = [0; e_funct(:,cn); 0];
    area(cn) = simpson1d((psi(:,cn) .* psi(:,cn))',xMin,xMax);
  if psi(5,cn) < 0, psi(:,cn) = -psi(:,cn); end  % curve starts positive
end % for
    psi = psi./sqrt(area);
    probD = psi.* psi;
    PROB = simpson1d(probD',xMin,xMax);

% SINGLE QUANTUM STATE M and EXPECTATION VALUES  ==================
   K = E(n) - U;          % kinetic energy  (eV)
   EB = -E(M);            % binding energy  (eV)
   bra = psi(:,M)';
   Psi = psi(:,M)';

% Probability
  ket = Psi;
  braket = bra .* ket;
  Prob = simpson1d(braket,xMin,xMax);

% Position x
  ket = x .* Psi;
  braket = bra .* ket;
  xavg = simpson1d(braket,xMin,xMax);

% Position x^2
  ket = (x.^2) .* Psi;
  braket = bra .* ket;
  x2avg = simpson1d(braket,xMin,xMax);

% Momentum ip      change length units from nm to m
  ket = gradient(Psi,x)* hbar / Ls;
  braket = bra .* ket;
  ipavg = simpson1d(braket,xMin,xMax);

% Momentum ip^2    chnage length unit from nm to m
  psi1 = gradient(Psi,x);
  ket = -gradient(psi1,x) * hbar^2 / Ls^2;
  braket = bra .* ket;
  ip2avg = simpson1d(braket,xMin,xMax);

% Potential energy U
  ket = U .* Psi;
  braket = bra .* ket;
  Uavg = simpson1d(braket,xMin,xMax);

% Kinetic energy K
  psi1 = gradient(Psi,x);
  ket = - (hbar^2/(2*me))* gradient(psi1,x) / (Ls^2*Es);
  braket = bra .* ket;
  Kavg = simpson1d(braket,xMin,xMax);

% total energy
  psi1 = gradient(Psi,x);
  ket = - (hbar^2/(2*me))* gradient(psi1,x) / (Ls^2*Es) + (U.*Psi);	
  braket = bra .* ket;
  Eavg = simpson1d(braket,xMin,xMax);

% Calculate uncertainites 
  deltax = sqrt(x2avg - xavg^2) * Ls;         % length units nm to m
  deltaip = sqrt(ip2avg + ipavg^2);            % unit N.s
  dxdp = (deltax * deltaip)/hbar;              % N.s.m = J.s   (SI unit)



% OUTPUT =============================================================
% Display eigenvalues in Command Window
  disp('  ');
  fprintf('No. bound states found =  %0.0f   \n',n);
  disp('   ');
  disp('Quantum State / Eigenvalues  En / Theory ET  (eV)');

for cn = 1 : n
    fprintf('  %0.0f   ',cn);
    fprintf('   %0.5g  \n ',E(cn));
end    

% Display results of calculations in Command Window
  disp('  ');
  disp('Single stationary state computations  ')
  fprintf('Stationary State n = M = %0.0f \n',M)
  fprintf('Total energy EM = %0.4f  eV \n',E(M))
  fprintf('Total Probability = %0.6g   \n',Prob);
  fprintf('<x> = %0.4f   nm\n',xavg);
  fprintf('<x^2> = %0.4f   nm^2\n',x2avg);
  fprintf('<ip> = %0.4e   N.s\n',ipavg);
  fprintf('<ip^2> = %0.4e   N^2.s^2\n',ip2avg);
  fprintf('deltax,dx = %0.4e  m \n',deltax);
  fprintf('deltap,dp = %0.4e   N.s\n',deltaip);
  fprintf('(dx dp)/hbar = %0.2f   \n',dxdp);
  fprintf('Uncertainty Principle (dx dp)/hbar >= 0.5\n');
  disp('  ')
  fprintf('<U> = %0.4f   eV\n',Uavg);
  fprintf('<K> = %0.4f   eV\n',Kavg);
  fprintf('<E> = %0.4f   eV\n',Eavg);
  fprintf('<K> + <U> = %0.4f   eV\n',Kavg+Uavg);
  disp('  ')
  

  % GRAPHICS  ======================================================

figure(1)  % 1111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.05 0.25 0.25]);
   set(gcf,'color','w');
   FS = 14;
   plot(x,U,'b','LineWidth',2);
    hold on
   plot([xavg xavg],[-600 400],'m','LineWidth',2)
   xlabel('x  [ nm ]')
   ylabel('|\psi)|^2')
   grid on; box on
   xlabel('position x  [ nm ]','FontSize',14);
   ylabel('U [ eV ]','FontSize',14);
   txt = sprintf('U_0 = %0.0f   x_1 = %0.2f   S = %0.2f',U0,x1,S);
   title(txt,'FontWeight','normal')
   set(gca,'fontsize',FS)
   
figure(2)  % 2222222222222222222222222222222222222222222222222222
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.4 0.25 0.5]);
   set(gcf,'color','w');
   FS = 14;
subplot(2,1,1)
   plot(x,U,'b','LineWidth',2);
   hold on
   cnmax = length(E);
   for cn = 1 : cnmax
       ys(1) = E(cn);
       ys(2) = ys(1);
       plot([xMin xMax],ys,'r','LineWidth',2);
  end 
   box on; grid on
   xlabel('position x  [ nm ]','FontSize',14);
   ylabel('U & E [ eV ]','FontSize',14);
   set(gca,'fontsize',FS)
subplot(2,1,2)
   Hplot = plot(1:length(E),E,'bo');
   set(Hplot,'MarkerFaceColor','b')
   box on; grid on
   xlabel('n','FontSize',FS);
   ylabel('E_n','FontSize',FS);
   set(gca,'fontsize',FS)
 
figure(3)  % 333333333333333333333333333333333333333333333333333
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.305 0.05 0.30 0.70]);
  set(gcf,'Color',[1 1 1]);
  FS = 14;
  nMax = 6;
  for cn = 1:nMax
    subplot(nMax,2,2*cn-1);
    y1 = psi(:,cn) ./ (max(psi(:,cn)-min(psi(:,cn))));
    y2 = 1 + 2 * U ./ (max(U) - min(U));
    plot(x,y1,'m','lineWidth',2)
    hold on
    plot(x,y2,'b','lineWidth',1)
    axis off
    title_m = ['\psi   n = ', num2str(cn)] ;
    title(title_m,'Fontsize',FS);
    
    subplot(nMax,2,2*cn);
    y1 = probD(:,cn) ./ max(probD(:,cn));
    y2 = 1 + 2 * U ./ (max(U) - min(U));
    plot(x,y1,'m','lineWidth',2)
    hold on
    plot(x,y2,'b','lineWidth',1)
    title_m = ['\psi^2   n = ', num2str(cn)] ;
    title(title_m,'Fontsize',14);
    axis off
  end
  
% Graphical display for quantum state M  % 44444444444444444444444444
figure(4) 
   set(gcf,'units','normalized');
   set(gcf,'position',[0.65 0.05 0.25 0.60]);
   set(gcf,'color','w');
   FS = 14;
subplot(3,1,1)
   plot(x,psi(:,M),'b-','linewidth',2)
   xlabel('x  [ nm ]')
   ylabel('\psi')
   txt = sprintf('n = %0.0f  ',M);
   title(txt,'FontWeight','normal')
   grid on
   set(gca,'FontSize',FS);
subplot(3,1,2)
   plot(x,psi(:,M).^2,'r-','linewidth',2)
   hold on
   plot([xavg xavg],[0 max(psi(:,M).^2)],'m','LineWidth',2)
   xlabel('x  [ nm ]')
   ylabel('|\psi)|^2')
   grid on
   set(gca,'FontSize',FS);   
subplot(3,1,3)
% Each point represents the location of the particle after
%    a measuerment is made on the system 
   hold on
   num1 = 10000;
   axis off

for c = 1 : num1
  xIndex = ceil(1+(NX-2)*rand(1,1));
  yIndex = rand(1,1);
  pIndex = max(probD(:,n))*rand(1,1);
   if pIndex < probD(xIndex,M)
   plot(x(xIndex),yIndex,'s','MarkerSize',2,'MarkerFaceColor','m','MarkerEdgeColor','none');
   end
end
   set(gca, 'Xlim',[xMin xMax]);
   set(gca, 'Ylim',[0,1]);
   set(gca,'FontSize',14)

   


