% qp_hydrogen.m
% Ian Cooper, School of Physics, The University of Sydney
% calls simpson1d.m to find area under curve

clear 
close all
clc


% INPUTS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('   ');
disp('Solution of the [1D] Radial Schrodinger Equation for hydrogen like atoms')
disp('    ')
r_max = input('max radial distance (default 10e-10 m), r_max =  ');
disp('   ');
L = input('orbital quantum number (default 0), L =  ');
disp('   ');
m_L = input('magnetic quantum number (default 0), m_L =  ');
disp('   ');
Z = input('nuclear charge (default 1), Z =  ');

disp('    ');
disp(' Please wait: calculating ....');

num = 1201;                             % number of data points default 1201


% CONSTANTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
eps0 = 8.854e-12;      % F/m


% POTENTIAL WELL +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% potential energy in electron volts (eV)
r_min = 1e-15;    % default 1e-15
r = linspace(r_min,r_max, num);
dr = r(2)-r(1);
dr2 = dr^2;

% Coulomb term
K = -Z*e/(4*pi*eps0);
U_c = K./r;

% Angular momnetum term
U_L = (hbar^2*L*(L+1)/(2*me*e))./r.^2;

% Effective potential energy
U = U_c + U_L;

U(U < -2000) = -2000;     % default value 1000
U(U > 1000) = 1000;

U_matrix = zeros(num-2,num-2);
for cn = 1:(num-2)
    U_matrix(cn,cn) = U(cn+1);
end



tic;
% SOLVE SCHRODINGER EQUATION +++++++++++++++++++++++++++++++++++++++++++++
% Make Second Derivative Matrix 
off     = ones(num-3,1);                 
SD_matrix = (-2*eye(num-2) + diag(off,1) + diag(off,-1))/dr2;

% Make KE Matrix
K_matrix = -hbar^2/(2*me*e) * SD_matrix;            

% Make Hamiltonian Matrix
H_matrix = K_matrix + U_matrix;

% Find Eignevalues E_n and Eigenfunctions psi_N 
[e_funct, e_values] = eig(H_matrix);

% All Eigenvalues 1, 2 , ... n  where E_N < 0
flag = 0;
n = 1;

while flag == 0
    E(n) = e_values(n,n);
    if E(n) > 0, flag = 1; end % if
    n = n + 1;
end  % while
E(n-1) = [];
n = n-2;

psi  = zeros(num,n);
prob = zeros(num,n);

% Corresponding Eigenfunctions 1, 2, ... ,n: Normalizing the wavefunction
for cn = 1 : n
psi(:,cn) = [0; e_funct(:,cn); 0];
area = simpson1d((psi(:,cn).* psi(:,cn))',r_min,r_max);
psi(:,cn) = psi(:,cn)/sqrt(area);
prob(:,cn) = psi(:,cn) .* psi(:,cn);
end % for


% CALCULATE EXPECTATION VALUES +++++++++++++++++++++++++++++++++++++++++++

disp('  ')
disp('Enter Principal Quantum Number for calculation of')
disp('   expectation values and graphical display')
pqn = input('Enter Principal Quantum Number (n > L),  n  =  ');
qn = pqn-L;


bra = psi(:,qn)';
Psi = psi(:,qn)';
Prob_density = Psi .* Psi;

% probability
ket = Psi;
braket = bra .* ket;
Prob = simpson1d(braket,r_min,r_max);

% position r
ket = r .* Psi;
braket = bra .* ket;
ravg = simpson1d(braket,r_min,r_max);
% position of prob_density peak
r_peak = r(Prob_density == max(Prob_density));

%position r^2
ket = (r.^2) .* Psi;
braket = bra .* ket;
r2avg = simpson1d(braket,r_min,r_max);

% momentum ip
ket = gradient(Psi,r)*hbar;
braket = bra .* ket;
ipavg = simpson1d(braket,r_min,r_max);

% momentum ip^2
psi1 = gradient(Psi,r);
ket = -gradient(psi1,r)*hbar^2;
braket = bra .* ket;
ip2avg = simpson1d(braket,r_min,r_max);

%potential energy U
ket = U .* Psi;
braket = bra .* ket;
Uavg = simpson1d(braket,r_min,r_max);

%kinetic energy K
psi1 = gradient(Psi,r);
ket = - (hbar^2/(2*me*e))*gradient(psi1,r);
braket = bra .* ket;
Kavg = simpson1d(braket,r_min,r_max);

%total energy
psi1 = gradient(Psi,r);
ket = - (hbar^2/(2*me*e))*gradient(psi1,r) + (U.*Psi);	
braket = bra .* ket;
Eavg = simpson1d(braket,r_min,r_max);

% Calculate uncertainites 
deltar = sqrt(r2avg - ravg^2);
deltaip = sqrt(ip2avg + ipavg^2);
drdp = (deltar * deltaip)/hbar;

E_max = max(abs(E));


% DISPLAY VALUES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
t_execution = toc;
clc
disp('---------------------------------------------------------------------');
fprintf('No. bound states found =  %0.0g   \n',n);
disp('  ');
disp('Quantum State / Eigenvalues  En  (eV)');
for cn = 1 : n
    fprintf('  %0.0f   ',L+cn);
    fprintf('   %0.5g   \n',E(cn));
end

disp('  ');
fprintf('Principal Quantum Number, n  =  %0.0g   \n',pqn);
fprintf('Orbital Quantum Number, L  =  %0.0g   \n',L);
fprintf('Magnetic Quantum Number, m_L  =  %0.0g   \n',m_L);
fprintf('Nuclear charge, Z  =  %0.0g   \n',Z);
disp('   ')
fprintf('Energy, E =  %0.6g   \n',E(qn));
fprintf('Total Probability = %0.6g   \n',Prob);
disp('  ')
fprintf('r_peak = %0.6g   m\n',r_peak);
fprintf('<r> = %0.6g   m\n',ravg);
fprintf('<r^2> = %0.6g   m^2\n',r2avg);

disp('  ')
fprintf('<ip> = %0.6g   N.s\n',ipavg);
fprintf('<ip^2> = %0.6g   N^2.s^2\n',ip2avg);
disp('  ')
fprintf('<U> = %0.6g   eV\n',Uavg);
fprintf('<K> = %0.6g   eV\n',Kavg);
fprintf('<E> = %0.6g   eV\n',Eavg);
fprintf('<K> + <U> = %0.6g   eV\n',Kavg+Uavg);
disp('  ')
fprintf('deltar = %0.6g   \n',deltar);
fprintf('deltaip = %0.6g   \n',deltaip);
fprintf('(dr dk)/hbar = %0.6g   \n',drdp);
disp('  ')
fprintf('execution time = %0.3g s  \n',t_execution);


% GRAPHICS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% plot potential energies
figure(1);   
  set(gcf,'Name','Potentials','NumberTitle','off')
  set(gcf,'Units','normalized','Position',[0.44 0.62 0.2 0.3]);
  plot(r.*1e9,U_c,'b','linewidth',2);
  hold on;
  plot(r.*1e9,U_L,'r','linewidth',2);
  plot(r.*1e9,U,'k','linewidth',2);
  legend('coul', 'Ang Mon', 'eff');
  grid on;
  axis([1e9.*r_min 1e9.*r_max -5*E_max 5*E_max]);
  xlabel('radial position  r  (nm)');
  ylabel('potential energy  U  (eV)');
  title('Potential Energy Functions','fontweight','normal');
  set(gca,'fontsize',12')
   
   
% Plot Potential Energy and eigenvalues -----------------------------
figure(2)
  set(gcf,'Name','Energy','NumberTitle','off')
  set(gcf,'Units','normalized','Position',[0.44 0.22 0.2 0.3]);
  plot(1e9.*r,U,'b','lineWidth',2);
  xlabel('radial position r (nm)');
  ylabel('energy U, E_n (eV)');
  title('Energy Spectrum','fontweight','normal');
  hold on
  for cn = 1 : n
    h_plot = plot([1e9.*r_min 1e9.*r_max], [E(cn) E(cn)]);
    set(h_plot, 'lineWidth',0.1, 'color', [1 0 0]);
  end
  axis([1e9.*r_min 1e9.*4*r_peak -5*E_max 10]);
  grid on
  set(gca,'fontsize',12);


% Plot Wave function and Probability function ----------------------------
figure(3)
fs = 12;
set(gcf,'Name','Wave Function','NumberTitle','off')
set(gcf,'Units','normalized','Position',[0.64 0.58 0.32 0.35]);

subplot(1,2,1)
   if psi(3,qn) > 0
       sf = 1;
   else
       sf = -1;
   end
   plot(1e9.*r,sf .* psi(:,qn),'LineWidth',2);
   grid on; hold on
   plot(1e9.*[r_peak r_peak],[0, max(sf .* psi(:,qn))],'r','linewidth',2)
   tm1 = 'r_{peak} =  ';
   tm2 = num2str(1e9.*r_peak,4);
   tm3 = '  nm';
   tm = [tm1 tm2 tm3];
   h_title = title(tm);
   set(h_title,'color',[1 0 0],'fontsize',fs,'fontweight','normal');
   ylabel('wave function,   \psi ','fontsize',fs);
   xlabel('position x (nm)','fontsize',fs);
   axis([1e9.*r_min 1e9.*4*r_peak -1.1*max(abs(psi(:,qn))) 1.1*max(abs(psi(:,qn)))]);
   set(gca,'fontsize',fs);
   set(gca,'Position',[0.066 0.12 0.38 0.71])

   subplot(1,2,2);
   plot(1e9.*r,prob(:,qn),'LineWidth',2);
   grid on
   xlabel('position x (nm)','fontsize',fs);
   ylabel('prob. density,  \psi^2 (m^{-1}) ','fontsize',14);
   hold on
   plot(1e9.*[r_peak r_peak],[0, max(Prob_density)],'r','linewidth',2)
   plot(1e9.*[ravg ravg],[0, max(Prob_density)],'m','linewidth',2)
   tm1 = 'r_{avg} =  ';
   tm2 = num2str(1e9.*ravg);
   tm3 = '  nm';
   tm = [tm1 tm2 tm3];
   title(tm);
   h_title = title(tm);
   set(h_title,'color',[1 0 1],'fontsize',fs,'fontweight','normal');
   axis([1e9.*r_min 1e9.*4*r_peak 0 1.1*max(abs(prob(:,qn)))]);
   set(gca,'fontsize',fs);
  
   
  
% Graphical output for <x>
figure(4)
set(gcf,'Name','<r>','NumberTitle','off')
set(gcf,'Units','normalized','Position',[0.64 0.12 0.32 0.35]);
fs = 12; 
subplot(3,1,1)
   plot(1e9.*r,-Psi,'LineWidth',2);
   grid on
   ylabel('bra = \psi','fontsize',fs);
   tm1 = '<r>  =  ';
   tm2 = num2str(1e9.*ravg,4);
   tm3 = '  nm';
   tm = [tm1 tm2 tm3];
   title(tm,'fontsize',fs,'fontweight','normal'); 
   axis([1e9.*r_min 1e9.*4*r_peak -1.1*max(abs(psi(:,qn))) 1.1*max(abs(psi(:,qn)))]);
   set(gca,'fontsize',fs)
   
subplot(3,1,2)
   plot(1e9.*r,ket,'LineWidth',2);
   grid on
   ylabel('ket ','fontsize',fs);
   axis([1e9.*r_min 1e9.*4*r_peak -1.1*max(abs(ket)) 1.1*max(abs(ket))]);
   set(gca,'fontsize',fs)
   
subplot(3,1,3);
   plot(1e9.*r,braket);
   fill(1e9.*r,braket,'r')
   grid on
   ylabel('braket ','fontsize',fs);
   hold on
   plot([ravg ravg],[min(braket) max(braket)],'lineWidth',3);
   xlabel('position r (nm)','fontsize',fs);
   hold off
   axis([1e9.*r_min 1e9.*4*r_peak -1.1*max(abs(braket)) 1.1*max(abs(braket))]);
   set(gca,'fontsize',fs)
   
   
% Probability Clouds ++++++++++++++++++++++++++++++++++++++++++++++++++++++
nump = 201;
maxp = r_max;

xxx = Prob_density;
EB = -E(pqn-L);

dp = 2*maxp/(nump-1);
xp = -maxp : dp : maxp; 
yp = -maxp : dp : maxp;

k = 0;

probcloud = zeros(nump,nump);


for c1 = 1 : nump
for c2 = 1 : nump
   
   rp = sqrt(xp(c1)^2 + yp(c2)^2);
   
if rp < maxp   
   
if xp(c1) >= 0 && yp(c2) >= 0
   theta = atan(abs(yp(c2)/(xp(c1)+eps)));
end
if xp(c1) <= 0 && yp(c2) >= 0
   theta = pi - atan(abs(yp(c2)/(xp(c1)+eps)));
end
if xp(c1) <= 0 && yp(c2) <= 0
   theta = pi + atan(abs(yp(c2)/(xp(c1)+eps)));
end
if xp(c1) >= 0 && yp(c2) <= 0
   theta = 2*pi - atan(abs(yp(c2)/(xp(c1)+eps)));
end

  leg01 = legendre(L,cos(theta));
  leg02 = leg01(m_L+1,:);

  N = round(1 + (num-1)*(rp/r_max));
  probcloud(c1,c2) = xxx(N)*leg02^2;    

end   %if   

end   %for c1
end   %for c2

probcloud = probcloud./max(max(probcloud));
probS = probcloud .^0.5;

%%
s = sprintf('n = %.0g    l = %.0g   m_l = %.0g     E_B = %0.4g  eV',pqn,L,m_L,EB);

figure(82);
  set(gcf,'Name','Cloud1','NumberTitle','off')
  set(gcf,'Units','normalized','Position',[0.02 0.62 0.2 0.3]);
  pcolor(xp.*1e9,yp.*1e9,probS);
  colorbar off;
  shading interp
  axis square
  title(s);
  xlabel('x  (nm)  ');
  ylabel('y   (nm) ');
  set(gca,'fontsize',12)
  
figure(83);
  set(gcf,'Name','Cloud2','NumberTitle','off')
  set(gcf,'Units','normalized','Position',[0.02 0.22 0.2 0.3]);
  pcolor(xp.*1e9,yp.*1e9,probS);
  colormap copper
  colorbar off;
  shading flat;
  title(s);
  axis square
  xlabel('x  (nm)  ');
  ylabel('y   (nm) ');
  set(gca,'fontsize',12)
  

figure(84)
  set(gcf,'Name','Cloud3','NumberTitle','off')
  set(gcf,'Units','normalized','Position',[0.23 0.62 0.2 0.3]);
  surf(xp,yp,probcloud);
  colorbar off;
  colormap jet
  shading interp;
  title(s);
  xlabel('x  (ao)  ');
  ylabel('y   (ao) ');
  axis off
  view(-86,70)
  set(gca,'fontsize',12)
  
figure(99)
  set(gcf,'Units','normalized','Position',[0.23 0.22 0.2 0.3]);
  fs = 12;
  set(gca,'fontsize',fs);
  subplot(1,2,1)
  xx = xp .* 1e9;  yy = xp .*1e9;
  pcolor(xx,yy,probS);
  colormap jet
  colorbar off;
  shading interp
  axis square
% title(s);
  xlabel('x  (nm)  ','fontsize',fs);
  ylabel('y   (nm) ','fontsize',fs);
% axis off
  set(gca,'fontsize',fs);

subplot(1,2,2)
surf(xp,yp,probcloud);
colorbar off;
shading interp;
s = sprintf('Electron density:  n = %.0g    L = %.0g   mL = %.0g     EB = %0.4g  eV',pqn,L,m_L,EB);
title(s,'fontsize',fs);
%xlabel('x  (ao)  ');
%ylabel('y   (ao) ');
axis off
set(gca,'fontsize',fs);
view(-86,70)


