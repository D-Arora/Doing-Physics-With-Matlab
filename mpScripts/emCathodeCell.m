% emCathodeCell.m

% CATHODE CELL: All particle simulations of cathodic arc plasmas 
% CATHODE SPOTS

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/emCathodeCell.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b
% 191227

clear
close all
clc

tic

% INPUTS ==============================================================

% Number of particles (must be even number)   [100]
  N = 100; 
% Number of electrons (first set of particles) and indices
  n1 = (1/2)*N;     n1Ind = 1:n1;
% Number of positive ions (second set of particles)
  n2 = (1/2)*N;     n2Ind = n1+1:N;  

% Number of time steps  [1e5]
  nt = 2e6;   
% Time step  [ < 1e-18  s]
  dt = 1e-18;
% Time interval between plot update in animation  [2e2]; 
  ng = 1e3;
% Time indices for time evoultion plots [plot 1000 points]
  tInd = 1:round(nt/5000):nt;

  % Atomic mass of positive ion   [Ti+ 47.9]
  A = 47.9;
% Initial cell size LxLxL   [1e-9 m] 
  L = 1e-9;
% Electron temperature  [K]
  TE = 3e4;
% Applied external electric field   [0 V/m]
  Ex = 0; 
% Soft potential: min dist between any two particles  [1e-11 m]
  rdMin = 1e-11 .*ones(N,N);   



% CHARGE: electrons 1  / positive ions 2    [C] 
  txtIon = 'Ti^{+}';
% Charge on positive ions  [1]
  nQ = 1;
  
% Fundamental charge [C]
  e = 1.60218e-19;
  q1 = -e;
  q2 = nQ*e;
  q = q1 .* ones(1,N);
  q(n2Ind) = q2;
  qq = meshgrid(q,q);


  
% SETUP ===============================================================

% CONSTANTS

% Electron mass [kg]
  mE = 9.10939e-31; 
% Nucleon mass [kg]
  mN = 1.674e-27; 
% Permittivity of free space  [SI units]
  eps0 = 8.85419e-12;  
% Coulomb constant  [SI units]
  k = 1/(4*pi*eps0);
% Boltzmaan constant [SI units]
  kB = 1.38066e-23;
  

% MASS: electrons / positive ions    [kg] 
  m = mE .* ones(1,N);  
  m(n2Ind) = A*mN;
  
% Calculation constants
  A1 = (k .* q) ./m; 
  A2 = A1 .* dt^2;
  
% Time grid
  t = 0 : dt : dt * (nt-1); 
  
  
% Electric field: acceleartion (Force per unit mass)  [kg/m   kg.s^2/m]
  Fx = (q .*Ex) ./ m;          
  FxT = Fx .* dt^2;
  
% Displacement components of displacement  [m]
  x = zeros(nt,N);
  y = zeros(nt,N);
  z = zeros(nt,N);
  
  rng('shuffle');
% Initial positions of particles (time step 1)
  x(1,:) = -L + (2*L) .* rand(1,N);
  y(1,:) = -L + (2*L) .* rand(1,N);
  z(1,:) = -L + (2*L) .* rand(1,N);
  
% VELOCITY components of particles  [m/s]
  vx = zeros(nt,N);
  vy = zeros(nt,N);
  vz = zeros(nt,N);
  
% ENERGIES {eV]
% kinetic energy   
  K = zeros(nt,N);
% PE: electron - electron interaction
  U1= zeros(nt,1);
% PE: ion - ion interaction
  U2 = zeros(nt,1);
% PE: electron - ion interaction
  U12 = zeros(nt,1);
  
% =====================================================================
% Maxwellian Distribution: hot electrons   cold protons
% Initial velocities of particles (time step 1]
%  vx(1,n1Ind) = -1e6 + 2e6 .* rand(1,n1);
%  vy(1,n1Ind) = -1e6 + 2e6 .* rand(1,n1);
%  vz(1,n1Ind) = -1e6 + 2e6 .* rand(1,n1);
% Grid points / scale factor / constant / flag
  num = 500;
  scale = 4;
  C = -0.5*mE/(kB*TE);
  flagM = 1;
% Most probable speed / max speed / sppeds 
  vProb = sqrt(2*kB*TE/mE);
  vMax = scale* vProb;
  vm = linspace(0,vMax,num);
  f = vm.^2 .* exp(C .* vm.^2);
  fMax = max(f);
  
  vM = zeros(1,n1);
  for c = 1 : n1
      while flagM == 1
          vR = vMax * rand;
          fC = (vR.^2 .* exp(C .* vR.^2)) ./fMax; 
           fR = rand;
            if fR < fC
              vM(c) = vR;
              flagM = 0;
            end
      end
  flagM = 1;
  end

  theta = -pi + (2*pi) .* rand(1,n1);
  phi   = -pi/2 + pi   .* rand(1,n1);
  
  vx(1,n1Ind) = vM .* cos(theta);
  vy(1,n1Ind) = vM .* sin(theta);
  vz(1,n1Ind) = vM .* sin(phi);
  
clear theta phi vM vR fC f vm 
  
% TIME EVOLUTION OF SYSTEM  ===========================================

% START  time step 1 
%  Position arrays
    xx = meshgrid(x(1,:),x(1,:));
    yy = meshgrid(y(1,:),y(1,:));
    zz = meshgrid(z(1,:),z(1,:));
    
%  Displacement between particle  [m]  
    xd = xx - xx';
    yd = yy - yy';
    zd = zz - zz';
    rd = sqrt(xd.^2 + yd.^2 + zd.^2);
    rd = rd + rdMin;
    rd3 = rd.^3;
    
% Potential energies [eV]    
    [U1(1), U2(1), U12(1)] = PE(rd,n1,n2);
    
% Time step 2
   Sx = (qq.*xd) ./rd3; 
   Sy = (qq.*yd) ./rd3; 
   Sz = (qq.*zd) ./rd3; 
   
% Acceleration components
   ax = -A1 .* sum(Sx') + Fx;
   ay = -A1 .* sum(Sy');
   az = -A1 .* sum(Sz');

   x(2,:) = x(1,:) + vx(1,:) .* dt + (0.5*dt^2) .* ax;
   y(2,:) = y(1,:) + vy(1,:) .* dt + (0.5*dt^2) .* ay;
   z(2,:) = z(1,:) + vz(1,:) .* dt + (0.5*dt^2) .* az;

clear ax ay az   
   
% Potential energies [eV]    
    U1(2) = U1(1); U2(2) = U2(1); U12(2)= U12(1);
   
   
% TIME STEPS > 2 
for c = 3 : nt
   xx = meshgrid(x(c-1,:),x(c-1,:));
   yy = meshgrid(y(c-1,:),y(c-1,:));
   zz = meshgrid(z(c-1,:),z(c-1,:));
   xd = xx - xx';
   yd = yy - yy';
   zd = zz - zz';
   rd = sqrt(xd.^2 + yd.^2 + zd.^2); 
   rd = rd + rdMin;  
   rd3 = rd.^3;
   Sx = (qq.*xd) ./rd3;
   Sy = (qq.*yd) ./rd3; 
   Sz = (qq.*zd) ./rd3;
   SSx = -A2 .* sum(Sx');
   SSy = -A2 .* sum(Sy');
   SSz = -A2 .* sum(Sz');    
    
   x(c,:) = 2.*x(c-1,:) - x(c-2,:) + SSx + FxT;
   y(c,:) = 2.*y(c-1,:) - y(c-2,:) + SSy;
   z(c,:) = 2.*z(c-1,:) - z(c-2,:) + SSz; 
   
   [U1(c), U2(c), U12(c)] = PE(rd,n1,n2);
end

clear Sx Sy Sz SSx SSy SSz A1 A2 qq xx yy zz xd yd zd rd rd3

% Distance of particles from Origin (0, 0)
  R1 = zeros(nt,n1); R1avg = zeros(nt,1);
  R2 = zeros(nt,n2); R2avg = zeros(nt,1);
for c = 1 : nt
  R1(c,n1Ind) = sqrt(x(c,n1Ind).^2 +  y(c,n1Ind).^2 + z(c,n1Ind).^2);
  R1avg(c) = mean(R1(c,:));
  
  R2(c,n2Ind) = sqrt(x(c,n2Ind).^2 +  y(c,n2Ind).^2 + z(c,n2Ind).^2);
  R2avg(c) = mean(R2(c,:));
end

% Final max distances from Origin  [nm]
  R1_max = max(R1(end,n1Ind))*1e9;
  R2_max = max(R2(end,n2Ind))*1e9;
  
% Velocities and kinetic energies   [J]
  for c = 2 : nt-1
    vx(c,:) = (x(c+1,:)- x(c-1,:))./(2*dt);
    vy(c,:) = (y(c+1,:)- y(c-1,:))./(2*dt);
    vz(c,:) = (z(c+1,:)- z(c-1,:))./(2*dt);
    K(c,:) = (0.5.*m) .* (vx(c,:).^2 + vy(c,:).^2 + vz(c,:).^2);
  end
    K(1,:) = (0.5.*m) .* (vx(1,:).^2 + vy(1,:).^2 + vz(1,:).^2);
    K(nt,:) = K(nt-1,:);
    
% Total and max final kinetic energies:  Start (S)   Finish (F)  [eV]    
    K1 = sum(K(:,n1Ind),2)./e;
    K2 = sum(K(:,n2Ind),2)./e;
    Ktot = K1 + K2;
    
    K_EmaxS = max(K(1,n1Ind))/e;
    K_PmaxS = max(K(1,n2Ind))/e;
    
    K_EmaxF = max(K(end,n1Ind))/e;
    K_PmaxF = max(K(end,n2Ind))/e;
    
% Potential Energies [eV]
    U1  = (k*q1*q1/e).* U1;
    U2  = (k*q2*q2/e).* U2;
    U12 = (k*q1*q2/e).* U12;

% Total energy  [eV]
    Utot = U1 + U2 + U12;
    E = Ktot + Utot;
  
    
% GRAPHICS ============================================================

%%
figure(1)
 s = 1e9;    % conversion m to nm
 SLim = 10; 
 
 set(gcf,'units','normalized');
 set(gcf,'position',[0.02 0.55 0.23 0.3]);
 set(gcf,'color','w');
 
 % Max time steps for initial trajectories
   cMax = 5e4;
 for c = 1:1e3:cMax
   plot3(s.*x(c,n1Ind), s.*y(c,n1Ind),s.*z(c,n1Ind),'k.')
   hold on
   Hplot = plot3(s.*x(c,n2Ind), s.*y(c,n2Ind),s.*z(c,n2Ind),'ro');
   set(Hplot,'markersize',5,'markerfacecolor','r')
   box on
   grid on
   
   %axis equal
   xlim([-SLim SLim]);
   ylim([-SLim SLim]);
   zlim([-SLim SLim]);
   set(gca,'xtick',-SLim:10:SLim)
   set(gca,'ytick',-SLim:10:SLim)
   set(gca,'ztick',-SLim:10:SLim)
   axis([-SLim SLim -SLim SLim -SLim SLim])
   
   tm = sprintf('Initial trajectories in first %3.2f pc',cMax*dt*1e12);
   title(tm) 
  % pause(0.01)
 end
  
   xlabel('{\itx}  [nm]','FontSize',12)
   ylabel('{\ity}  [nm]','FontSize',12)
   zlabel('{\itz}  [nm]','FontSize',12)
   set(gca,'Fontsize',12);
   
%%   
 figure(2)
   s = 1e9;    % conversion m to nm
   SLim = 100; 
   set(gcf,'units','normalized');
   set(gcf,'position',[0.02 0.1 0.23 0.3]);
   set(gcf,'color','w');
 
  plot3(s.*x(tInd,n1Ind), s.*y(tInd,n1Ind), s.*z(tInd,n1Ind),'k.')
  hold on
  plot3(s.*x(tInd,n2Ind), s.*y(tInd,n2Ind), s.*z(tInd,n2Ind),'r.')
  Hplot = plot3(s.*x(end,n1Ind), s.*y(end,n1Ind), s.*z(end,n1Ind),'ko');
    set(Hplot,'markersize',5,'markerfacecolor','k')
  Hplot = plot3(s.*x(end,n2Ind), s.*y(end,n2Ind), s.*z(end,n2Ind),'ro');
    set(Hplot,'markersize',5,'markerfacecolor','r')
     
   box on
   grid on
   axis equal
   xlim([-SLim SLim]);
   ylim([-SLim SLim]);
   zlim([-SLim SLim]);
   set(gca,'xtick',-SLim:50:SLim)
   set(gca,'ytick',-SLim:50:SLim)
   set(gca,'ztick',-SLim:50:SLim)
   %view(0,0)
 
   title('Particle Trajectories','fontweight','normal');
   xlabel('{\itx}  [nm]','FontSize',12)
   ylabel('{\ity}  [nm]','FontSize',12)
   zlabel('{\itz}  [nm]','FontSize',12)
   set(gca,'Fontsize',12);
 
%%   
 figure(3)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.26 0.45 0.23 0.40]);
   set(gcf,'color','w');
  
   s = 1e9;    % conversion m to nm
   st = 1e12;  % conversion s to ps
   subplot(2,1,1)
   plot(st.*t(tInd), s.*R1avg(tInd),'k','linewidth',2)
   box on
   grid on
   ylabel('{\itR_{e avg}}  [ nm ]','FontSize',12)
   set(gca,'Fontsize',12);
   title('Average Particle Displacements','fontweight','normal');
   
   subplot(2,1,2)
   plot(st.*t(tInd), s.*R2avg(tInd),'r','linewidth',2)
   box on
   grid on
   xlabel('{\it} t  [ ps ]','FontSize',12)
   ylabel('{\itR_{p avg}}  [ nm ]','FontSize',12)
   set(gca,'Fontsize',12);
   
%%   
 figure(4)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.50 0.45 0.23 0.4]);
   set(gcf,'color','w');
  
   subplot(2,1,1)
   plot(st.*t(tInd), K1(tInd)./(N/2),'k','linewidth',2)
   box on
   grid on
   ylabel('{\it KE_e}  [ eV ]','FontSize',12)
   set(gca,'Fontsize',12);  
   title('Average Particle KE','fontweight','normal')
   
   subplot(2,1,2)
   plot(st.*t(tInd), K2(tInd)./(N/2),'r','linewidth',2)
   box on
   grid on
   xlabel('{\it t}  [ ps ]','FontSize',12)
   ylabel('{\it KE_p}  [ eV ]','FontSize',12)
   set(gca,'Fontsize',12);
   
   

 figure(5)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.74 0.45 0.23 0.40]);
   set(gcf,'color','w');
  
   plot(st.*t(tInd), U1(tInd),'k','linewidth',2)
   hold on
   plot(st.*t(tInd), U2(tInd),'r','linewidth',2)
   plot(st.*t(tInd), U12(tInd),'b','linewidth',2)
   plot(st.*t(tInd), Utot(tInd),'m','linewidth',2)
   plot(st.*t(tInd), E(tInd),'g','linewidth',2)
   
   box on
   grid on
   HL = legend('U_E','U_P','U_{EP}','U_{tot}','E');
   set(HL,'orientation','horizontal','location','northoutside')
   xlabel('{\it t} [ ps ]','FontSize',12)
   ylabel('{\it U and E}  [ eV ]','FontSize',12)
   set(gca,'Fontsize',12);  
  

figure(6)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.50 0.10 0.23 0.4]);
   set(gcf,'color','w');
  
  subplot(2,1,1)
   Hplot = plot(s*R1(end,n1Ind), K(end,n1Ind)./e,'ko');
   set(Hplot,'markersize',4,'markerfacecolor','k')
   hold on
  % plot(s*R1(1,n1Ind), K(1,n1Ind)./e,'k+')
   box on
   grid on
   ylabel('{\it KE_e}  [ eV ]','FontSize',12)
   set(gca,'Fontsize',12);  
  title('Final Kinetic Energies','fontweight','normal')
   
  subplot(2,1,2)
   Hplot = plot(s*R2(end,n2Ind), K(end,n2Ind)./e,'ro');
   set(Hplot,'markersize',4,'markerfacecolor','r')
   hold on
   %plot(s*R2(1,n1Ind), K(1,n2Ind)./e,'r+')
   box on
   grid on
   xlabel('{\it R}  [ nm ]','FontSize',12)
   ylabel('{\it KE_p}  [ eV ]','FontSize',12)
   set(gca,'Fontsize',12);   

   
% DISPLAY RESULTS -----------------------------------------------------

%%
figure(7)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.26 0.10 0.35 0.60]);
   set(gcf,'color','w');

   xlim([0 160])
   ylim([0 250])
   h = 250; dh = -18;
   
   txt = txtIon;
   text(2,h,txt,'fontsize',12)
   txt = sprintf('Cell size  L = %3.0f  nm', L*s);
   text(18,h,txt,'fontsize',12)
    txt = sprintf('No. time steps nt = %3.0f', nt);
   text(72,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('Time increment dt = %3.2e s', dt);
   text(2,h,txt,'fontsize',12)
   txt = sprintf('Simultation time = %3.2f ps', 1e12*t(end));
   text(82,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('No. electrons = %3.0f', n1);
   text(2,h,txt,'fontsize',12)
   txt = sprintf('No. ions = %3.0f', n2);
   text(52,h,txt,'fontsize',12)
   txt = sprintf('Ion density = %3.2e  /m^3', n2/L^3);
   text(97,h,txt,'fontsize',12)
   h = h+dh;
   
   
   
   txt = sprintf('INITIAL ENERGY VALUES PER PARTICLE  [ eV ]');
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('    U_E = %3.1f    U_P = %3.1f    U_{EP} = %3.1f', ...
       U1(1)/n1, U2(1)/n2, U12(1)/N);
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('    K_{Eavg} = %3.1f    K_{Pavg} = %3.1f    K_{avg} = %3.1f     E = %3.1f', ...
       K1(1)/n1, K2(1)/n2, Ktot(1)/N, E(1)/N );
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
    txt = sprintf('    K_{Emax} = %3.1f    K_{Pmax} = %3.1f', ...
       K_EmaxS, K_PmaxS);
   text(2,h,txt,'fontsize',12)
   h = h+dh; 
   
  txt = sprintf('FINAL ENERGY VALUES PER PARTICLE  [ eV ]');
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('    U_E = %3.1f    U_P = %3.1f    U_{EP} = %3.1f', ...
       U1(end)/n1, U2(end)/n2, U12(end)/N);
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('    K_{Eavg} = %3.1f    K_{Pavg} = %3.1f    K_{avg} = %3.1f     E = %3.1f', ...
       K1(end)/n1, K2(end)/n2, Ktot(end)/N, E(1)/N );
   text(2,h,txt,'fontsize',12)
   h = h+dh; 
   
   txt = sprintf('    K_{Emax} = %3.1f    K_{Pmax} = %3.1f', ...
       K_EmaxF, K_PmaxF);
   text(2,h,txt,'fontsize',12)
   h = h+dh; 
   
   txt = sprintf('FINAL DISPLACEMENTS FROM ORIGIN   [ nm ]');
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('    R_{Eavg} = %3.1f    R_{Pavg} = %3.1f', ...
       1e9*R1avg(end), 1e9*R2avg(end));
   text(2,h,txt,'fontsize',12)
   h = h+dh;
   
   txt = sprintf('    R_{Emax} = %3.1f    R_{Pmax} = %3.1f', ...
       R1_max, R2_max);
   text(2,h,txt,'fontsize',12)
   
   h = h+dh;
   
   txt = sprintf('Execution time %3.2f  min',round(toc/60));
   text(2,h,txt,'fontsize',12)
   
  axis off
  
   
%%   
   
 
   
   
toc

% =====================================================================
% FUNCTIONS
% =====================================================================

function [U1, U2, U12] = PE(rd,n1,n2) 
  

  M = 1./rd;

  U1  = triu(M(1:n1,1:1:n1)); 
  U2  = triu(M(n1+1:n1+n2,n1+1:n1+n2)); 
  U12 = M(1:n1,n1+1:n1+n2); 

  U1 = U1 - U1(1,1).*eye(n1,n1);
  U2 = U2 - U2(1,1).*eye(n2,n2);
  
  U1  = sum(U1,'all');
  U2  = sum(U2,'all');
  U12 = sum(U12,'all');

end

