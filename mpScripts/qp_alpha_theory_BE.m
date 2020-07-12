% qp_alpha_theory_BE.m

% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% ALPHA DECAY: relationship between half-lives and kinetic energy for
%  polonium (1) and uranium (2)
%   Experimental values /  Theoetical values / Simulation values
% Comment / uncomment lines for theoretical calculation of half-lives for
%     experimental energies T1 and T2

% 20190325   Matlab R2018b

close all
clear 
clc


% INPUT PARAMETERS ====================================================

% flagE = 0  graphical output for T = Tmin to Tmax;
% flagE = 1  display polonium data
% flagE = 2  display uranium data
flagE = 0;

% Range of kinetic enegies of emitted alpha particle  [4 to 10 Mev}
Tmin = 4;
Tmax = 10;

% Radius factor for calculation of radius of alpha and daughter  [fm]
%   Polonium data (1) and uranium data (2)   
Rf = [1.07 1.12];

% Number of grid point for calculations
N = 10; 


% Constants ===========================================================
   e = 1.602176565e-19;    % elementary charge
   hbar = 1.054571726e-34; % hbar Planck's constant
   mA = 6.64465675e-27;    % alpha particle mass 
   eps0 = 8.854187e-12;    % permittivity of free space


% Data: energy T [J]  /  half-life h [s] ============================== 
%   Po Z = 84 
      T1 = [6.020 6.537 6.778 7.37 7.687 8.376 8.785 5.304 4.883 5.11]'.*e.*1e6;
 %  Experimental data      
      h1  = [186 1.53 0.15 1.80e-3 1.62e-4 3.70e-6 3.0e-7 1.20e7 3.22e9 9.15e7]'; 
%   Simulation #1 data      
      h11 = [420 3.5e-1 2.4e-2 4.5e-3 4.6e-4 6.4e-6 9.0e-8 3.2 1.4e6 5.6e10 5.4e6]';

%   U  Z = 92
      T2 = [4.151 4.445 4.596 4.722 4.729 5.236 5.818 6.410 7.060 7.402 7.875 8.780]'.*e.*1e6;  
%   Experimental data
        h2 = [1.41e17 7.39e14 2.22e16 7.75e12 5.02e12 2.23e9 1.80e6 5.46e2 66 0.26 84e-3 18e-6]';
 %   Simulation data #1
        h21 = [4.7e16 2.8e14 2.3e14 3.7e13 3.0e13 1.0e10 4.4e6 9.9e2 3.4 0.41 3.3e-3 4.4e-5]';
     
 

% CALCULATIONS ========================================================
% Atomic Number  /  Mass number: polonium Po (1) and uranium U (2)
   Z = [84 92];
   A = [218 238];

% Radius of alpha plus daughter [m]
   R0 = 1e-15.*Rf.*(4^(1/3) + (A-4).^(1/3));    
   
% Range of kinetic energies of emitted alpha particle [J]
  T = 1e6*e*linspace(Tmin,Tmax,N)';
  if flagE == 1; T = T1; end
  if flagE == 2; T = T2; end
     
% Reduced mass [kg]
   mR = mA.*(A-4)./A;

% Velocity of alpha particle when escapes from nucleus [m/s]
   vA = sqrt(2.*T./mR);     

% Frequency of alpha particle inside nucleus striking potential barrier
   f = vA ./ (2.*R0);    % [1/s]    

% Transmission probability of an alpha particle escaping nucleus
   K1 = (4*e/hbar).*(mR./(pi*eps0)).^0.5 .* R0.^0.5 .* (Z-2).^0.5; 
   K2 = (e^2/(hbar * eps0)) .* (mR./2).^0.5 .* (Z-2) ;
   P = exp(K1 - K2.* T.^(-0.5));
% decay constant  [1/s]   
   gamma = f.* P;

% Half-life [s]   
   t_half = log(2) ./ gamma;


 tH = TP(K1,K2,f,gamma,T); 
   
%  TABLE OUTPUT ======================================================= 
  if flagE == 1
   E_alpha = T1./(e.*1e6);
   t_alpha = h1;
   table(E_alpha, t_alpha)
    disp('  ')
   E_alpha = T./(e.*1e6);
   t_alpha = t_half;
   table(E_alpha, t_alpha(:,1))
   disp('polonium:  energy  [MeV]  /  half-life [s]')
  end

  if flagE == 2
   E_alpha = T2./(e.*1e6);
   t_alpha = h2;
   table(E_alpha, t_alpha)
   disp('  ')
   E_alpha = T./(e.*1e6);
   t_alpha = t_half;
   table(E_alpha, t_alpha(:,2))
  end

   

if flagE == 0
% GRAPHICS ============================================================
figure(1)
  pos = [0.05 0.1 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  
% Polonium theoretical values
  xP = T.^(-0.5); yP = log(t_half(:,1));
  yP = log(tH(:,1));
    plot(xP,yP,'b','linewidth',2);
  hold on
% Uranium theoretical values  
  xP = T.^(-0.5); yP = log(t_half(:,2));
    plot(xP,yP,'r','linewidth',2);
    
% Polonium experimental values
 xP = T1.^(-0.5); yP = log(h1);
  plot(xP,yP,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',8);
 xP = T1.^(-0.5); yP = log(h11);
 % plot(xP,yP,'o','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1],'MarkerSize',6);  
  
% Uranium experimental values
 xP = T2.^(-0.5); yP = log(h2);
  plot(xP,yP,'o','MarkerFaceColor',[1 0 0],'MarkerSize',8,'MarkerEdgeColor',[1 0 0])
 xP = T2.^(-0.5); yP = log(h21);
%   plot(xP,yP,'o','MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[1 1 0],'MarkerSize',6);  

  legend('Po','U','location','northwest');
  xlabel('1 / \surdT')
  ylabel('log_{e}( t_{1/2} )','fontsize',14)
  set(gca,'fontsize',12)
  grid on

end

function tH = TP(K1,K2,f,gamma,T)
  % transmission pobability
    P = exp(K1 - K2.* T.^(-0.5));
% decay constant  [1/s]   
   gamma = f.* P;
% Half-life [s]   
   tH = log(2) ./ gamma;
end

