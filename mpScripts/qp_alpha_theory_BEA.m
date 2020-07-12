% qp_alpha_theory_BE.m

% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

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


% Po Z = 84   A = 218 217 216 215 214 213 212 210 209 208
% Po Experimental kinetic energies T [J] 
     T1 = [6.002 6.537 6.778 7.37 7.687 8.376 8.785 5.304 4.883 5.11]'.*e.*1e6;
% Po Experimental half-lives [s]     
     h1  = [186 1.53 0.15 1.80e-3 1.62e-4 3.70e-6 3.0e-7 1.20e7 3.22e9 9.15e7]'; 
% Po Numerical simulations half-lives  [s]      
     h11 = [28e2 3.6 1.4e-1 43e-3 41e-4 56e-6 15e-7 0.63e7 0.59e9 8.4e7]';
    
% U  Z = 92   A = 238 236 234 233 232 230 228 227 226 225 223
% U  Experimental kinetic energies T [J] 
     T2 = [4.151 4.445 4.596 4.722 4.729 5.236 5.818 6.410 7.060 7.402 7.875 8.780]'.*e.*1e6;  
% U  Experimental half-lives [s]
     h2 = [1.41e17 7.39e14 2.22e16 7.75e12 5.02e12 2.23e9 1.80e6 5.46e2 66 0.26 84e-3 18e-6]';
% U  Numerical simulation data half-lives  [s]
     h21 = [0.003e17 0.25e14 0.14e14 2.1e12 1.7e12 0.63e9 0.33e6 1.2e2 0.4 0.26e-1 0.26e-3 4.8e-6]';
 

% CALCULATIONS ========================================================
% Atomic Number  /  Mass number: polonium Po (1) and uranium U (2)
   Z = [84 92];
   A = [218 238];

% Radius of alpha plus daughter [m]
   R0 = 1e-15.*Rf.*(4^(1/3) + (A-4).^(1/3));    
   
% Range of kinetic energies of emitted alpha particle [J]
  T = 1e6*e*linspace(Tmin,Tmax,N)';
      
% Reduced mass [kg]
   mR = mA.*(A-4)./A;

% Transmission probability of an alpha particle escaping nucleus
    K1 = (4*e/hbar).*(mR./(pi*eps0)).^0.5 .* R0.^0.5 .* (Z-2).^0.5; 
    K2 = (e^2/(hbar * eps0)) .* (mR./2).^0.5 .* (Z-2) ;
%    P = exp(K1 - K2.* T.^(-0.5));
% % decay constant  [1/s]   
%    gamma = f.* P;

% Half-life [s]   
%    t_half = log(2) ./ gamma;


 tH  = caltH(R0,mR,K1,K2,T); 
 tH1 = caltH(R0,mR,K1,K2,T1); 
 tH2 = caltH(R0,mR,K1,K2,T2);
  
%  TABLE OUTPUT ======================================================= 
T_Po =  T1./(e*1e6); tH_Po = tH1(:,1);  
table(T_Po,tH_Po)
disp('Po  Theoretical values')
   
T_U =  T2./(e*1e6); tH_U = tH2(:,2);  
table(T_U,tH_U) 
disp('U  Theoretical values')  
 

   


% GRAPHICS ============================================================
figure(1)
  pos = [0.05 0.1 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  
% Polonium theoretical values
  xP = T.^(-0.5); yP = log(tH(:,1));
    plot(xP,yP,'b','linewidth',2);
  hold on
% Uranium theoretical values  
  xP = T.^(-0.5); yP = log(tH(:,2));
    plot(xP,yP,'r','linewidth',2);
    
% Polonium experimental values
 xP = T1.^(-0.5); yP = log(h1);
  plot(xP,yP,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',8);
 xP = T1.^(-0.5); yP = log(h11);
  plot(xP,yP,'o','MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1],'MarkerSize',6);  
  
% Uranium experimental values
 xP = T2.^(-0.5); yP = log(h2);
  plot(xP,yP,'o','MarkerFaceColor',[1 0 0],'MarkerSize',8,'MarkerEdgeColor',[1 0 0])
 xP = T2.^(-0.5); yP = log(h21);
  plot(xP,yP,'o','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1],'MarkerSize',6);  

  legend('Po','U','location','northwest');
  xlabel('1 / \surdT')
  ylabel('log_{e}( t_{1/2} )','fontsize',14)
  set(gca,'fontsize',12)
  grid on


function tH = caltH(R0,mR,K1,K2,T)
   % Velocity of alpha particle when escapes from nucleus [m/s]
     vA = sqrt(2.*T./mR);     
   % Frequency of alpha particle inside nucleus striking potential barrier
     f = vA ./ (2.*R0);    % [1/s]   
   % transmission pobability
     P = exp(K1 - K2.* T.^(-0.5));
   % decay constant  [1/s]   
     gamma = f.* P;
% Half-life [s]   
   tH = log(2) ./ gamma;
end

