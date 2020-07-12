% qpAlphaDecayEAN.m

% ALPHA DECAY HALF-LIFE COMPUTATION 
%   using the Finite Difference Method to solve the [1D] time independent
%   Schrodinger Equation for a nuclear potential. 

% VARIABLES:  potential energy U, kinetic energy T
%             transition probability P, radial coordinates x
%             half life h
%             polonium P (1) / uranium U (2)
%             experimental E / analytical A / numerical N

% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% Badr-Eddine Elfatehy Chfira: badreddine.el@um.es
% University of Murcia (UMU, Spain) 

% DOING PHYSICS ONLINE: 
% ../mphome.htm
% Documentation
% http://www.physics.usyd.edu.au/teach_res/mp/doc/qpAlphaDecay.pdf



close all
clear 
clc


% CONSTANTS ===============================================================
   h_N = 6.62607004e-34;      % Planck constant [J.s]
   hbar = h_N/(2*pi);         % [J.s]
   e = 1.602e-19;             % Fundamental charge  [C]
   mA = 6.64465675e-27;       % mass alpha particle nucleus  [kg]
   eps0 = 8.854187e-12;       % permittivity of free space
   Ese = e*1e6;               % Energy scaling factor [J  MeV] 
   Lse = 1e-15;               % Length scaling factor [m fm]
   k = 1/(4*pi*eps0);         % Coulomb constant   
   num = 501;                 % Grid point for geiger-Nuttall Law

% Numerical simulations
  % Depth of potenital well [-40 MeV]
    U0 = -40;  
  % Max radial coordinate [250 fm]
    xMax = 250;    
  % Number of grid points   
    N = 1e6;                   
 
   
% POLONIUM (1)
  Z1 = 84;
  A1 = [218 217 216 215 214 213 212 210 209 208]';
  Rf1 = 1.07;   % nuclear radius factor
% Po Experimental kinetic energies T [J] 
     T1_E = [6.002 6.537 6.778 7.37 7.687 8.376 8.785 5.304 4.883 5.11]'.*Ese;
% Po Experimental half-lives [s]     
     h1_E  = [186 1.53 0.15 1.80e-3 1.62e-4 3.70e-6 3.0e-7 1.20e7 3.22e9 9.15e7]'; 

% URANIUM (2)
  Z2 = 92;
  A2 = [238 236 235 234 233 232 230 228 227 226 225 223]';
  Rf2 = 1.15;    % nuclear radius factor
% U  Experimental kinetic energies T [J] 
     T2_E = [4.151 4.445 4.215 4.722 4.729 5.236 5.818 6.410 6.860 7.402 7.875 8.780]'.*Ese;  
% U  Experimental half-lives [s]
     h2_E = [1.41e17 7.39e14 2.22e16 7.75e12 5.02e12 2.23e9 1.80e6 5.46e2 66 0.26 84e-3 18e-6]';
     

% INPUTS =================================================================
 
for flagE = 1:4

% Polonium    
if flagE == 1
  Z = Z1;
  A = A1;
  T =  T1_E./Ese;    % KE alpha particle [J]  
  Rf = Rf1;          % Nuclear radius factor   
end

% Uranium
if flagE == 2
  Z = Z2;
  A = A2;
  T =  T2_E./Ese;   % KE alpha particle  [J] 
  Rf = Rf2;         % Nuclear radius factor    
end

% Polonium: Geiger-Nutall Law
if flagE == 3 
   num = 501;
   Z = Z1;
   A = 214.*ones(num,1);
   T = linspace(4,9,num)';
   Rf = Rf1;
end

% Uranium: Geiger-Nutall Law
if flagE == 4 
   num = 501;
   Z = Z2;
   A = 230.*ones(num,1);
   T = linspace(4,9,num)';
   Rf = Rf2;
end

% SETUP ===============================================================
% Radial coordinate x
   xMin = 0;
   x = linspace(xMin,xMax,N);
   dx = x(2)-x(1);
  
% Wavefunction: numerical simulations   psi = zeros(1,N);  
   psi = zeros(N,1);
   psi(N) = 1; psi(N-1) = 2; 
% Potential energy function: numerical simulations  [MeV]   
   U = U0 .* ones(1,N); 
% Half-lives: humerical h_N / h_A analytical      
   h_N = zeros(length(A),1);
   h_A = zeros(length(A),1);
  

  
% FDM =================================================================   
 for c = 1: length(A) 
 % reduced mass  [kg] 
  mR = mA*(A(c)-4)/A(c);
 % Potential well width (radius of alpha + daughter)
   x0 = Rf * (4^(1/3) + (A(c)-4)^(1/3));
 % Position at U = T; 
  x1 = 2*(Z-2)*e^2/(4*pi*eps0*T(c)*Ese*Lse);
 % Constants 
  Cse = hbar^2/(2*mR) / (Lse^2*Ese);     
  K1 = (4*e/hbar).*(mR./(pi*eps0)).^0.5 .* (Lse.*x0).^0.5 .* (Z-2).^0.5; 
  K2 = (e^2/(hbar * eps0)) .* (mR./2).^0.5 .* (Z-2) ; 
    
% Nuclear potential
   for n = 2:N
     if x(n) >  x0
       U(n) = k * 2 * (Z-2) * e^2 / (Ese * Lse * x(n));
     end
   end

% Wavefunction
  for n = N-1:-1:2    
   SEconst = Lse^2 .* Ese .* (2*mR/hbar^2).*(T(c) - U(n)).* dx^2;
   psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
  end
  
  
if flagE < 3  
% CALCULATIONS ========================================================
% NUMERICAL
% Position in arrays
% index for x position when E = U outside well  
   index = find(x>x1, 1 );  
% index for position outside well   
   index1 = find(x>x0, 1 );
% index for postion to calc Aout   
   index2 = round(0.9*N);       
% velocity of alpha particle inside well
   v_in = sqrt(2*(T(c)- U0)*Ese/mR);
% frequency at which particles strike barrier on the RHS of well   
   f =  v_in / (2*x0*Lse);               
% amplitude of wavefunction inside and outsie nucleus   
   Ain  = max(abs(psi(1:index1)));     
   Aout = max(abs(psi(index2:N)));
% probability of alpha particle escaping
   P = (Aout / Ain)^2;  
% decay constant
   gamma = f * P; 
% half-life   
   h_N(c) = log(2) / gamma;   
end   

% ANALYTICAL
% Theoretical values   
   P = exp(K1 - K2.* (T(c).*Ese).^(-0.5)); 
% decay constant  [1/s]   
     gamma = f.* P;
% Half-life [s]   
   h_A(c) = log(2) ./ gamma;
 end

% STORE VALUES ======================================================== 
if flagE == 1
    A1 = num2str(A,'%3.0f \n');
    T1 = T;
    h1_N = h_N;
    h1_A = h_A;
    
end
if flagE == 2
    A2 = num2str(A,'%3.0f \n');
    T2 = T;
    h2_N = h_N;
    h2_A = h_A;
end

if flagE == 3
    T1_GN = T;
    h1_GN = h_A;
end

if flagE == 4
    T2_GN = T;
    h2_GN = h_A;
end

end  
% DISPLAY OUTPUTS =======================================================
  format SHORTE
% Polonium  
  TP   = num2str(T1_E./Ese,'%2.3f \n');
  hP_E = num2str(h1_E,'%2.2e \n');
  hP_A = num2str(h1_A,'%2.1e \n');
  hP_N = num2str(h1_N,'%2.1e \n');
  RP_EA = num2str(h1_E./h1_A,'%2.1f \n');
  RP_EN = num2str(h1_E./h1_N,'%2.1f \n');
  table(A1,TP,hP_E,hP_A,RP_EA,hP_N,RP_EN)

% Uranium  
  TU   = num2str(T2_E./Ese,'%2.3f \n');
  hU_E = num2str(h2_E,'%2.2e \n');
  hU_A = num2str(h2_A,'%2.1e \n');
  hU_N = num2str(h2_N,'%2.1e \n');
  RU_EA = num2str(h2_E./h2_A,'%2.1f \n');
  RU_EN = num2str(h2_E./h2_N,'%2.1f \n');
    table(A2,TU,hU_E,hU_A,RU_EA,hU_N,RU_EN)
 

%  GRAPHICS ============================================================
figure(1)
  pos = [0.05 0.1 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  hold on
  box on
 
% Polonium: GeigerNuttall Law
   xP = (Ese.*T1_GN).^(-0.5); yP = log(h1_GN);
      plot(xP,yP,'b','lineWidth',2)
 % Uranium:  GeigerNuttall Law   
    xP = (Ese.*T2_GN).^(-0.5); yP = log(h2_GN);
      plot(xP,yP,'r','lineWidth',2)     
    
% Polonium: Experimental values
  xP = T1_E.^(-0.5); yP = log(h1_E);
    %plot(xP,yP,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],'MarkerSize',8)
    plot(xP,yP,'o','linewidth',2,'MarkerEdgeColor',[0 0 1],'MarkerSize',8)    
% Polonium: Analytical values
  xP = T1_E.^(-0.5); yP = log(h1_A);
    plot(xP,yP,'o','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1],'MarkerSize',8)
% Polonium: Numerical values
   xP = T1_E.^(-0.5); yP = log(h1_N);
    plot(xP,yP,'s','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',8)
            
% Uranium: Experimental values
  xP = T2_E.^(-0.5); yP = log(h2_E);
    %plot(xP,yP,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],'MarkerSize',8)
    plot(xP,yP,'o','linewidth',2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8)    
% Uranium: Analytical values
  xP = T2_E.^(-0.5); yP = log(h2_A);
    plot(xP,yP,'o','MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1],'MarkerSize',8)
% Polonium: Numerical values
   xP = T2_E.^(-0.5); yP = log(h2_N);
    plot(xP,yP,'s','MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[0 0 0],'MarkerSize',8) 

    
  xlabel('1 / \surdT    (T [J]) ')
  ylabel('log_{e}( t_{1/2} )','fontsize',14) 
  legend('Po','U','location','northwest');
  grid on
  set(gca,'fontsize',14)  

figure(2)
  pos = [0.40 0.1 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  hold on
  box on

% Polonium: GeigerNuttall Law
   xP = T1_GN; yP = log10(h1_GN);
      plot(xP,yP,'b','lineWidth',2)
 % Uranium:  GeigerNuttall Law   
    xP = T2_GN; yP = log10(h2_GN);
      plot(xP,yP,'r','lineWidth',2)     
  
  
% Polonium experimental values
 xP = T1_E./Ese; yP = log10(h1_E);
  plot(xP,yP,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',8);
  
% Uranium experimental values
 xP = T2_E./Ese; yP = log10(h2_E);
  plot(xP,yP,'o','MarkerFaceColor',[1 0 0],'MarkerSize',8,'MarkerEdgeColor',[1 0 0])

 grid on
 xlabel('T   [MeV ]) ')
 ylabel('log_{10}( t_{1/2} )','fontsize',14) 
 legend('Po','U','location','northeast');
 grid on
 set(gca,'fontsize',14)  


