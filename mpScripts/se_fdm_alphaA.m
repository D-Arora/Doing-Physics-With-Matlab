% se_fdm_alphaA.m

% ALPHA DECAY HALF-LIFE COMPUTATION 
%   using the Finite Difference Method to solve the [1D] time independent
%   Schrodinger Equation for a nuclear potential. 

% VARIABLES:  potential energy U, kinetic energy T
%             transition probability P, radial coordinates x.

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% Documentation
% http://www.physics.usyd.edu.au/teach_res/mp/doc/qp_alpha_decay_BE.pdf


close all
clear 
clc

% INPUTS =================================================================

% Atomic number of element
Z = 92;
% Mass number of element
A = 236;
% Kinetic energy of escaped alpha particle  [MeV]
T = 6.17;
% Nuclear radius factor  polonium (1.07) uranium (1.15) 
Rf = 1.15;

% Angular momentum quantum number LA [0, 1, 2, ... ]
LA = 5;

% Woods Saxon Nuclear Potential  flagWS = 1 / square well flagWS = 0;
%   diffuseness parameter a   [0.25 to 1]
flagWS = 1;
a = 1;

% Depth of potential well [-40 MeV]
   U0 = -40;  
% Max radial coordinate [250 fm]
   xMax = 250;    
% Number of grid points   
   N = 1e6;                   

   
% CONSTANTS ===============================================================
   h = 6.62607004e-34;        % Planck constant [J.s]
   hbar = h/(2*pi);           % [J.s]
   e = 1.602e-19;             % Fundamental charge  [C]
   mA = 6.64465675e-27;       % mass alpha particle nucleus  [kg]
   eps0 = 8.854187e-12;       % permittivity of free space

   mR = mA*(A-4)/A;                     % reduced mass  [kg]
   Ese = 1.6e-13;                       % Energy scaling factor [J  MeV] 
   Lse = 1e-15;                         % Length scaling factor [m fm]
   Cse = hbar^2/(2*mR) / (Lse^2*Ese);   % Schrodinger Eq constant       
   k = 1/(4*pi*eps0);                   % Coulomb constant   
 
   
% SETUP CALCULATIONS ==================================================
% Radial coordinate x
   xMin = 0;
   x = linspace(xMin,xMax,N);
   dx = x(2)-x(1);
% Potential well width (radius of alpha + daughter)
   x0 = Rf * (4^(1/3) + (A-4)^(1/3));
% Potentail energy = KE alpha particle   
 %  x1 = 2*(Z-2)*e^2/(4*pi*eps0*T*Ese*Lse);

   
% FDM =================================================================
   psi = zeros(1,N);  
   psi(N) = 1; psi(N-1) = 2; 
   U = U0 .* ones(1,N);
   
% Woods Saxon Nuclear Potential
   if flagWS == 1
     U = ( U0./(1+exp(-x0./a))).* ones(1,N); 
   end
   
% Nuclear potential
   for n = 2:N
     if x(n) >  x0
       U(n) = k * 2 * (Z-2) * e^2 / (Ese * Lse * x(n));
       U(n) = U(n) + (LA.*(LA+1).*Cse./((x(n)).^2));
         if flagWS == 1
           U(n) = U(n)+ U0./(1+exp((x(n)-x0)./a));
         end
     end
   end

% Wavefunction
  for n = N-1:-1:2    
    SEconst = (T - U(n)).* dx^2./Cse;
    psi(n-1) = (2 - SEconst) * psi(n) - psi(n+1);
  end

  psi = psi ./ max(psi);
  
% Second turning point for each value of the orbital angular momentum, LA= 0, 1 ,2, 3
%   T = U
x1 = max(double(solve(@(z) LA.*(LA+1).*Cse./(z.^2) + 2.*(Z-2).*e.^2./ (4.*pi.*eps0.*Ese.* Lse.* z) == T))); % U(R1) = T   
  
  
% CALCULATIONS ========================================================
% Position in arrays
% index for x position when E = U outside well  
   index = find(x>x1, 1 );  
% index for position outside well   
   index1 = find(x>x0, 1 );
% index for postion to calc Aout   
   index2 = round(0.9*N);       
% velocity of alpha particle inside well
   v_in = sqrt(2*(T-U0)*Ese/mR);
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
   half_life = log(2) / gamma;    

%  Wavelength --> kinetic energy
   [~, z2] = findpeaks(psi(1:index1));
   wLin = (x(z2(end))-x(z2(1)))/(length(z2)-1);
   [z1, z2] = findpeaks(psi(index2:N));
   wLout = (x(z2(end)+index2)-x(z2(1)+index2))/(length(z2)-1);

   p_in = h/(wLin*1e-15);
   T_in = p_in^2/(2*mR)/(e*1e6);
  
   p_out = h/(wLout*1e-15);
   T_out = p_out^2/(2*mR)/(e*1e6);
   
 % kinetic energy T - U 
   K_in = T - U0;   
   K_out = T - U(end);
   
   Umax = max(U);
   
   
% DISPLAY OUTPUTS =======================================================
  table(A,Z,T,Rf,N,LA,a,half_life,gamma,P)
  table(x0,Umax,x1,v_in,f,T_in,T_out,K_in,K_out)
  
  
% GRAPHICS ==============================================================
figure(1)
  set(gcf,'color',[1 1 1]);
  set(gcf,'Units','Normalized') 
  set(gcf,'Position',[0.1 0.1 0.4 0.8]) 

subplot(4,1,1)
  plot(x,U,'r','lineWidth',3)
  hold on
  plot([x(1) x(end)],[T T],'b','linewidth',2);
  ylabel('T, U (MeV)','fontsize',14);
  plot([x(index) x(index)],[-abs(U0) abs(U0)],'k','lineWidth',1);
  plot([0 0],[-40 40],'r','lineWidth',3);
  xlim([0 100])
  grid on
  tm1 = 'T = ';
  tm2 = num2str(T, '%3.3f  MeV');
  tm3 = '    t_{1/2} = ';
  tm4 = num2str(half_life, '%3.1e  s');
  tm = [tm1 tm2 tm3 tm4];
  title(tm)
  
  set(gca,'fontsize',14);

subplot(4,1,2)
  set(gcf,'color',[1 1 1]);
  plot(x,psi,'lineWidth',3)
  hold on
  set(gca,'Ylim',[-1.2*max(psi) 1.2*max(psi)]);
  plot([x(index) x(index)],[-1.2*max(psi) 1.2*max(psi)],'k','lineWidth',1);
  set(gca,'fontsize',14);
  ylabel('\psi ','fontsize',16);
  grid on

subplot(4,1,3)
  set(gcf,'color',[1 1 1]);
  plot(x,psi,'lineWidth',3)
  hold on
  set(gca,'Ylim',[-1.2*max(psi) 1.2*max(psi)]);
  plot([x(index) x(index)],[-1.2*max(psi) 1.2*max(psi)],'k','lineWidth',1);
  set(gca,'fontsize',14);
  ylabel('\psi ','fontsize',16);
  xlim([0 20])
 grid on

subplot(4,1,4)
 set(gcf,'color',[1 1 1]);
 plot(x,psi,'lineWidth',3)
  hold on
  plot([xMin xMax],[0 0],'k');
  xlim([200 250])

  set(gca,'Ylim',[-1.2*max(psi(index:N)) 1.2*max(psi(index:N))]);
  set(gca,'fontsize',14);
  xlabel('radial coordinate   \it {x}  [ fm ]','fontsize',14)
  ylabel('\psi ','fontsize',16);
  grid on



