% qmN001.m

close all
clc
clear


% *********************************************************************
% INPUTS
% *********************************************************************
% Number of grid points (must be an ODD number)
   N = 999;
% x Range  [nm]   
   xMin = -0.10;  xMax = 0.10;
% Percentage width of well  [19 < x < 90]
   w = 90;
% Well depth
  U0 = 400;
% Graphical display for the state with quantum qn
   qn = 2;
%  Orthonormal startes  qn1  qn2
   qn1 = 1;  qn2 = 2;
   
   
% *********************************************************************
% SETUP
% *********************************************************************
  x = linspace(xMin,xMax,N);
  dx = x(2)-x(1);
  dx2 = dx^2;

% Elementary charge [C]  
    e = 1.6021766208e-19;
% Planck constant [ J.s] / Reduced Planck's constant
    h = 6.626070040e-34;            
    hbar = 1.054571800e-34; 
% Electron mass [kg] 
    me = 9.10938356e-31; 
% Scaling constants
    Es = e;
    Ls = 1e-9;
    Cs = -hbar^2/(2*me) / (Ls^2*Es); 


% *********************************************************************
% POTENTIAL ENERGY MATRIX
% *********************************************************************
% Square Well dimensions: boundary wB NB and wW Nw 
%     L = abs(xMin) + abs(xMax);   
%     w = (w/100)*L;  wB = (L-w)/2;
%     Nw = (w/L)*N; NB = round((wB/L)*N,0);
%     NR = (NB :(N-NB))';
%     
%     U = zeros(N,1);
%     U_matrix = zeros(N-2,N-2);
%     U(NR) = -U0;
%     
%   for cn = 1:(N-2)
%       U_matrix(cn,cn) = U(cn+1);
%   end
  
 % Potential well 
   x1 = 0.04;                 % width of truncated potential well [nm]
   U0 = -400;                % well depth [eV]
  num = N;
  U = zeros(num,1);
  U_matrix = zeros(num-2);
    
   for cn = 1 : num
   if abs(x(cn))<=x1/2,   U(cn) = -(4*U0/(x1*x1))*x(cn)^2+U0; end
   end  

% Make potential energy matrix
  dx = (x(2)-x(1));
  dx2 = dx^2;
  for cn = 1:(num-2)
    U_matrix(cn,cn) = U(cn+1);
  end



% *********************************************************************
% KINETIC ENERGY and HAMILTONIAN MATRICES
% *********************************************************************  
% Second Derivative Matrix 
    off = ones(N-3,1);                 
    SD_matrix = (-2*eye(N-2) + diag(off,1) + diag(off,-1))/dx2;
% KE Matrix
    K_matrix = Cs * SD_matrix;            
% Hamiltonian Matrix
   H_matrix = K_matrix + U_matrix;

   
% *********************************************************************
% EIGENVALUES and EIGENFUNCTIONS
% *********************************************************************    
   [e_funct, e_values] = eig(H_matrix);

% All Eigenvalues 1, 2 , ... n  where E_N < 0
   flag = 0;
   n = 1;
while flag == 0
    E(n) = e_values(n,n);
    if E(n) > 0, flag = 1; end 
    n = n + 1;
end  
    E(n-1) = [];
    n = n-2;

% Corresponding Eigenfunctions 1, 2, ... ,n: Normalizing the wavefunction
    psi = zeros(N,n); area = zeros(1,n);
for cn = 1 : n
    psi(:,cn) = [0; e_funct(:,cn); 0];
    area(cn) = Simpson1d((psi(:,cn) .* psi(:,cn))',xMin,xMax);
  if psi(5,cn) < 0, psi(:,cn) = -psi(:,cn); end  % curve starts positive
end % for
    psi = psi./sqrt(area);
    probD = psi.* psi;

    PROB = simpson1d(probD',xMin,xMax);
    
% *********************************************************************
% SINGLE QUANTUM STATE and EXPECTATION VALUES  qn  
% *********************************************************************
   K = E(n) -U;          % kinetic energy  (eV)
   EB = -E(qn);            % binding energy  (eV)

  bra = psi(:,qn)';
  Psi = psi(:,qn)';

% probability
  ket = Psi;
  braket = bra .* ket;
  Prob = simpson1d(braket,xMin,xMax);

% position x
  ket = x .* Psi;
  braket = bra .* ket;
  xavg = simpson1d(braket,xMin,xMax);

% position x^2
  ket = (x.^2) .* Psi;
  braket = bra .* ket;
  x2avg = simpson1d(braket,xMin,xMax);

% momentum ip      change length units from nm to m
  ket = gradient(Psi,x)* hbar / Ls;
  braket = bra .* ket;
  ipavg = simpson1d(braket,xMin,xMax);

% momentum ip^2    chnage length unit from nm to m
  psi1 = gradient(Psi,x);
  ket = -gradient(psi1,x) * hbar^2 / Ls^2;
  braket = bra .* ket;
  ip2avg = simpson1d(braket,xMin,xMax);

%potential energy U
  ket = U' .* Psi;
  braket = bra .* ket;
  Uavg = simpson1d(braket,xMin,xMax);

% kinetic energy K
  psi1 = gradient(Psi,x);
  ket = - (hbar^2/(2*me))* gradient(psi1,x) / (Ls^2*Es);
  braket = bra .* ket;
  Kavg = simpson1d(braket,xMin,xMax);

% total energy
  psi1 = gradient(Psi,x);
  ket = - (hbar^2/(2*me))* gradient(psi1,x) / (Ls^2*Es) + (U.*Psi')';	
  braket = bra .* ket;
  Eavg = simpson1d(braket,xMin,xMax);

% Calculate uncertainites 
deltax = sqrt(x2avg - xavg^2) * Ls;         % length units nm to m
deltaip = sqrt(ip2avg + ipavg^2);            % unit N.s
dxdp = (deltax * deltaip)/hbar;              % N.s.m = J.s   (SI unit)


% Two quantum numbers for two states qn1  qn2
%    to test normalization or that two different 
%    eigenvectors are orthogonal to each other orthogonal
%    Evaluate integral -------------------------------------------
     bra = psi(:,qn1)';
     ket = psi(:,qn2)';
     braket = bra .* ket;
     integral = round(simpson1d(braket,xMin,xMax));
     
     
% *********************************************************************
% Display eigenvalues in Command Window
% *********************************************************************
  disp('   ');
  disp('================================================================  ');
  disp('  ');
  fprintf('No. bound states found =  %0.0g   \n',n);
  disp('   ');
  disp('Quantum State / Eigenvalues  En  (eV)');

for cn = 1 : n
    fprintf('  %0.0f   ',cn);
    fprintf('   %0.5g   \n',E(cn));
end

% Display results of calculations in Command Window
  disp('  ');
  disp('  ');
  fprintf('Quantum number, n  =  %0.0g   \n',qn);
  fprintf('Energy, E =  %0.6g  eV \n',E(qn));
  fprintf('Total Probability = %0.6g   \n',Prob);
  fprintf('<x> = %0.6g   nm\n',xavg);
  fprintf('<x^2> = %0.6g   nm^2\n',x2avg);
  fprintf('<ip> = %0.6g   N.s\n',ipavg);
  fprintf('<ip^2> = %0.6g   N^2.s^2\n',ip2avg);
  fprintf('<U> = %0.6g   eV\n',Uavg);
  fprintf('<K> = %0.6g   eV\n',Kavg);
  fprintf('<E> = %0.6g   eV\n',Eavg);
  fprintf('<K> + <U> = %0.6g   eV\n',Kavg+Uavg);
  disp('  ')
  fprintf('deltax = %0.6g  m \n',deltax);
  fprintf('delta| = %0.6g   N.s\n',deltaip);
  fprintf('(dx dp)/hbar = %0.6g   \n',dxdp);
  disp('  ')
  fprintf('Integral =  %0.1f  \n',integral);
  if integral < 0.5
    disp('Wavefunctions are orthogonal');
  else
    disp('Wavefunction is normalized');
  end 
   
% **********************************************************************
% GRAPHICS  
% **********************************************************************
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.05 0.2 0.25]);
   set(gcf,'color','w');
   FS = 14;
   
   plot(x,U,'b','LineWidth',2);
   set(gca,'fontsize',14)
   xlabel('position x (nm)','FontSize',14);
   ylabel('energy U, E_n (eV)','FontSize',14);
   title('SQUARE WELL','FontSize',12,'fontweight','normal');
   hold on

  cnmax = length(E);
for cn = 1 : cnmax
  ys(1) = E(cn);
  ys(2) = ys(1);
  plot([xMin xMax],ys,'r','LineWidth',2);
end %for   
  axis([xMin-eps xMax min(U)-50 max(U)+50]);

% Plots first 5 wavefunctions & probability density functions
if n < 6
    nMax = n;
else
    nMax = 5;
end


figure(2)
  clf
  set(gcf,'Units','Normalized');
  set(gcf,'Position',[0.26 0.05 0.25 0.35]);
 % set(gcf,'NumberTitle','off');
 % set(gcf,'Name','Eigenvectors & Prob. densities');
  set(gcf,'Color',[1 1 1]);
  
  %nMax = 8;
  for cn = 1:nMax
    subplot(nMax,2,2*cn-1);
    y1 = psi(:,cn) ./ (max(psi(:,cn)-min(psi(:,cn))));
    y2 = 1 + 2 * U ./ (max(U) - min(U));
    plot(x,y1,'m','lineWidth',2)
    hold on
    plot(x,y2,'b','lineWidth',1)
    %plotyy(x,psi(:,cn),x,U);
    axis off
    %title('\psi cn);
    title_m = ['\psi   n = ', num2str(cn)] ;
    title(title_m,'Fontsize',14);
    
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

  
% Graphical display for quantum state qn
figure(3)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.05 0.4 0.2 0.25]);
   set(gcf,'color','w');
   FS = 14;
   left_color = [0 0 1];
   right_color = [1 0 1];
   
   yyaxis left
   plot(x,U,'b-','linewidth',2)
   hold on
   plot(x,K,'k-','linewidth',2);
   plot([xMin xMax], [-EB -EB],'r-','linewidth',2)
   legend('U','K','E','Orientation','horizontal','location','north','AutoUpdate','off')
   
   ylim([-1.1*U0 1.5*U0])
   txt = sprintf('qn = %2.0f  \n',qn);
   title(txt,'fontweight','normal')
   grid on
   set(gca,'YColor',left_color);
   
   yyaxis right
   plot(x,psi(:,qn),'m','linewidth',2)
   hold on
   plot([xMin xMax],[0 0],'m-','linewidth',0.5)
   ylabel('   \psi','fontsize',16,'rotation',0)
   set(gca,'ytick',[])
   set(gca,'YColor',right_color);
  
   
figure(4)
   set(gcf,'Name','Schrodinger Equation: Bound States');
   set(gcf,'NumberTitle','off');
   set(gcf','Units','normalized')
   set(gcf,'Position',[0.55 0.05 0.25 0.40]);
   set(gcf,'Color',[1 1 1]);
   eft_color = [0 0 1];
   right_color = [1 0 1];
   
   axes('position',[0.1 0.48 0.75 0.35]);
   
   yyaxis left
   plot(x,U,'b-','linewidth',2)
   hold on
    plot([xMin xMax], [-EB -EB],'r-','linewidth',2)
   legend('U','E','Orientation','horizontal','location','north','AutoUpdate','off')
   
   ylim([-1.1*U0 10])
   xlabel('x  [nm]')
   txt = sprintf('qn = %2.0f  \n',qn);
   title(txt,'fontweight','normal')
   grid on
   set(gca,'YColor',left_color);
   set(gca,'fontsize',14)
    
   yyaxis right
   plot(x,probD(:,qn),'m','linewidth',2)
   hold on
   plot([xMin xMax],[0 0],'m-','linewidth',0.5)
   ylabel('     \psi^2','fontsize',16,'rotation',0)
   set(gca,'ytick',[])
   set(gca,'YColor',right_color);
   

% Each point represents the location of the particle after
%    a measeurment is made on the system 
   axes('position',[0.15 0.05 0.75 0.25]);
   hold on
   num1 = 10000;
   axis off

for c = 1 : num1
  xIndex = ceil(1+(N-2)*rand(1,1));
  yIndex = rand(1,1);
  pIndex = max(probD(:,n))*rand(1,1);
   if pIndex < probD(xIndex,qn)
   plot(x(xIndex),yIndex,'s','MarkerSize',2,'MarkerFaceColor','m','MarkerEdgeColor','none');
   end
end
   set(gca, 'Xlim',[xMin xMax]);
   set(gca, 'Ylim',[0,1]);

hold off
  
figure(5)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.28 0.5 0.25 0.30]);
   set(gcf,'color','w');
   FS = 12;
   
   bra = psi(:,qn)';
   ket = x .* psi(:,qn)';
   braket = bra .* ket;
    
subplot(3,1,1)
   plot(x,Psi,'m','LineWidth',2);
   grid on
   ylabel('bra  \psi');
   tm1 = '<x>  =  ';
   tm2 = num2str(xavg,'%2.2f');
   tm3 = '  nm';
   tm = [tm1 tm2 tm3];
   title(tm);  
   set(gca,'fontsize',FS)
   
subplot(3,1,2)
   plot(x,ket,'b','LineWidth',2);
   grid on
   ylabel('ket  x \psi ');
   set(gca,'fontsize',FS)
   
subplot(3,1,3);
   plot(x,braket,'b');
   fill(x,braket,'r')
   grid on
   ylabel('bra ket  \psi x \psi');
   hold on
   plot([xavg xavg],[min(min(braket)) max(max(braket))],'lineWidth',3);
   xlabel('position x (nm)');
   hold off
   set(gca,'fontsize',FS)
   
   
% **********************************************************************
% FUNCTIONS  
% **********************************************************************

function integral = Simpson1d(f,a,b)
  % [1D] integration - Simpson's 1/3 rule
  %      f function    a = lower bound    b = upper bound
  %      Must have odd number of data points
  %      Simpson's coefficients   1 4 2 4 ... 2 4 1

    numS = length(f);               % number of data points
  if mod(numS,2) == 1
      sc = 2*ones(numS,1);
      sc(2:2:numS-1) = 4;
      sc(1) = 1; sc(numS) = 1;
      h = (b-a)/(numS-1);
      integral = (h/3) * f * sc;
  else
      integral = 'Length of function must be an ODD number'; 
  end
end



