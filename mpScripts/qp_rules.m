% qp_rules.m

% HYDROGEN ATOM: Selection Rules and Transition Rates
% Calls qp_fh.m to find solution to radial equation numerically, giving it
% n and l to construct l dependent effective potential, returning the nth
% of many solutions found.
% Calls simpson1d.m to find area under curve:  num must be an ODD number
% The plot for the prob density can be changed by commenting /
%   uncommenting out the lines of code for the plot.
% The code for the animated gif can be commented if you do not want it
%   saved.
% HydrogenElectricDipoleTransitionRates
% Calculates lifetime = 1/transitionrate for low lying states of hydrogen
%    in the electric dipole approximation by numerical integration of
%    transition dipole moments

% Supporting documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/qp_rules.pdf

% Ian Cooper
% School of Physics, University of Sydney
% Email:  ian.cooper@sydney.edu.au
% Doing Physics with Matlab
%      https://d-arora.github.io/Doing-Physics-With-Matlab/

% Corrections:
%    Duncan Carlsmith
%    Professor of Physics
%    University of Wisconsin-Madison
%    duncan@hep.wisc.edu
% 180315

clear all
close all
clc


% Inputs ===============================================================
% State #2  State #1   n L Mm_L
    state = [4 0 0 2 1 1];
% r_max: Adjust to give accurate results of radial wavefunction
    r_max = 10e-10;      
% Show animation of probability cloud in Fig. Window (yes 1) (no 0)
    flagS = 1;
% Probability cloud: (1 surf plot)  (0 contour plot)
    flagPC = 1;    
% Save animation as gif  (yes 1)   (no 0)
    f_gif = 0;

% must be an ODD number
    num = 1901;               

% Setup for saving images (im) ==========================================
   ag_name = 'ag_rules2.gif';
% Delay in seconds before displaying the next image  
   delay = 1;  
% Frame counter start
   nt = 1; 
   
   
% CALCULATIONS ==========================================================
   n(1) = state(1);   n(2) = state(4);
   L(1) = state(2);   L(2) = state(5);
   mL(1) = state(3);  mL(2) = state(6);

% CONSTANTS -------------------------------------------------------------
h = 6.626e-34;         % J.s
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
eps0 = 8.854e-12;      % F/m
c = 3e8;               % m./s

% Initialize arryas ------------------------------------------------------
N = zeros(1,2);
IP = zeros(1,3);
IT = zeros(1,3);
IR = zeros(1,3);
Ixyz = zeros(1,3);


% Azimuthal ==============================================================
phi = linspace(0,2*pi,num)';
PHI = zeros(num,2);

f_PHI = zeros(num,3);

PHI(:,1) = exp((1i*mL(1)).* phi)./sqrt(2*pi);
PHI(:,2) = exp((-1i*mL(2)).* phi)./sqrt(2*pi);

PHI_PHI = PHI(:,2) .* PHI(:,1);

f_PHI(:,1) = cos(phi) .* PHI_PHI;
f_PHI(:,2) = sin(phi) .* PHI_PHI;
f_PHI(:,3) = PHI_PHI;

% Dipole integral - azimuthal
for cc = 1 : 3
    IP(cc) = (simpson1d((f_PHI(:,cc))',0,2*pi));
end

% Angular ================================================================
theta = linspace(0,pi,num)';
THETA = zeros(num,2);

THETA1 = legendre(L(1),cos(theta))';  % returns P_L^{|m|}
THETA2 = legendre(L(2),cos(theta))';

THETA(:,1) = THETA1(:,abs(mL(1))+1);
THETA(:,2) = THETA2(:,abs(mL(2))+1);

sinT = sin(theta);

% Normalize angular wavefunctions ---------------------------------------
for cc = 1 : 2
   N(cc) = simpson1d((sinT .* THETA(:,cc) .* THETA(:,cc))',0,pi);
   THETA(:,cc) = THETA(:,cc) / sqrt(N(cc));
end

% Dipole integral - angular
sin2T = sinT.^2;
cosT = cos(theta);
THETA_THETA = THETA(:,1) .* THETA(:,2);

f_THETA(:,1) = sinT .* THETA_THETA .*sinT;
f_THETA(:,2) = f_THETA(:,1);
f_THETA(:,3) = cosT .* THETA_THETA .* sinT;

% Dipole integral - angular
for cc = 1 : 3
    IT(cc) = (simpson1d((f_THETA(:,cc))',0,pi));
end

% RADIAL WAVEFUNCTION ===================================================
[ EB(1), R1, r] = qp_fh(n(1), L(1),num, r_max);
[ EB(2), R2, r] = qp_fh(n(2), L(2),num, r_max);

% Normalize radial wavefunctions
N(1) = simpson1d((R1 .* R1),0,max(r));
N(2) = simpson1d((R2 .* R2),0,max(r));
R1 = R1 ./ sqrt(N(1));
R2 = R2 ./ sqrt(N(2));

f_R = r .* R1 .* R2;

% Dipole integral - radial
for cc = 1 : 3
   IR(cc) = simpson1d(f_R,0,max(r));
end

% DIPOLE INTEGRAL / ELECTRIC DIPOLE MOMENT / TRANSITION RATE

for cc = 1: 3
   Ixyz(cc) = IR(cc) * IP(cc) * IT(cc);    
end

Idipole = sqrt(sum(Ixyz*Ixyz'));

p21 = e .* Idipole;

f = abs((EB(1) - EB(2)) * e /h);
wL = 1e9 * c /f;
TR = 16 * pi^3 * p21^2 * f^3 / (3 * eps0 * h * c^3);
lifeTime = 1e9/TR;



% GRAPHICS ==============================================================

if flagS == 1
figure(1);   % probability cloud
   set(gcf,'Name','Cloud2','NumberTitle','off')
   set(gcf,'Units','normalized','Position',[0.34 0.1 0.25 0.4]);
   set(gcf,'color',[1 1 1]);

for k = 1 : 20
    
if k < 11; a1 = 1-k/10;  a2 = 1-a1;
else
    a2 = 1 - (k-10)/10; a1 = 1-a2;
end
Rmax = max(r);
numR = 201;
R = linspace(0,Rmax,numR);
nump = length(R);
maxp = max(R);
xxx = zeros(1,numR);
probcloud = zeros(nump,nump);

dp = 2*maxp/(nump-1);
xp = -maxp : dp : maxp; 
yp = -maxp : dp : maxp;

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

N = round(1 + (num-1)*(rp/max(r)));
leg01 = legendre(L(1),cos(theta));
leg01 = leg01(mL(1)+1,:);
psi1 = R1(N) .* leg01;

leg01 = legendre(L(2),cos(theta));
leg01 = leg01(mL(2)+1,:);
psi2 = R2(N) .* leg01;

psi = a1 .* psi1 + a2 .* psi2;

probcloud(c1,c2) = psi .* psi;  

else
    probcloud(c1,c2) = 0;

end   %if  rp < maxp 

end   %for c1
end   %for c2

%probcloud = probcloud./max(max(probcloud));
 probS = probcloud .^0.5;

   set(gca,'ZLim',[0 60000]);
   pcolor(xp,yp,probS);
   if flagPC == 1
      surf(xp,yp,probS);
      view(-45,72)
   end
   colormap copper
   colorbar off;
   shading flat;
   axis square
   axis off
  
   tm1 = num2str(k);
   tm2 = '       ';
   tm3 = num2str(state(1:3));
   tm4 = num2str(state(4:6));
   tm = [tm1 tm2 tm3 tm2 tm4];
   title(tm);
   set(gca,'fontsize',14)
   pause(0.001)

    if f_gif > 0 
       frame = getframe(1);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
     %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
    end
   
end  % for k
end

figure(99)  % numerical results
  pos = [0.05 0.05 0.28 0.52];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  axis([0 100 0 100])
  
  tm1 = 'Initial State (n  l m_l)   ';
  tm2 = num2str(state(1),'%2.0f  ');
  tm3 = num2str(state(2),'%2.0f  ');
  tm4 = num2str(state(3),'%2.0f  ');
  tm = [tm1 tm2 tm3 tm4];
  text(0,100,tm,'fontsize',14)
  tm1 = '  Binding Energy  E_B = ';
  tm2 = num2str(EB(1),'%4.3f eV ');
  tm = [tm1 tm2];
  text(0,90,tm,'fontsize',14)
  
  tm1 = 'Final State (n  l m_l)   ';
  tm2 = num2str(state(4),'%2.0f  ');
  tm3 = num2str(state(5),'%2.0f  ');
  tm4 = num2str(state(6),'%2.0f  ');
  tm = [tm1 tm2 tm3 tm4];
  text(0,80,tm,'fontsize',14)
  tm1 = '  Binding Energy  E_B = ';
  tm2 = num2str(EB(2),'%4.3f eV ');
  tm = [tm1 tm2];
  text(0,70,tm,'fontsize',14)
  
  tm1 = 'Photon wavelength  \lambda = ';
  tm2 = num2str(wL,'%4.1f nm ');
  tm = [tm1 tm2];
  text(0,60,tm,'fontsize',14)
  
  tm1 = 'Integral  I_P = ';
  tm2 = num2str(real(IP),'%3.4f ');
  tm = [tm1 tm2];
  text(0,50,tm,'fontsize',14)
  
  tm1 = 'Integral  I_T = ';
  tm2 = num2str(IT,'%3.4f ');
  tm = [tm1 tm2];
  text(0,40,tm,'fontsize',14)
  
  tm1 = 'Integral  I_R = ';
  tm2 = num2str(IR,'%3.3e ');
  tm = [tm1 tm2];
  text(0,30,tm,'fontsize',14)
  
  tm1 = 'Integral  I_{xyz} = ';
  tm2 = num2str(real(Ixyz),'%3.3e    ');
  tm = [tm1 tm2];
  text(0,20,tm,'fontsize',14)
  
  tm1 = 'Transition Rate  T_R = ';
  tm2 = num2str(TR,'%3.3e  s^{-1} ');
  tm = [tm1 tm2];
  text(0,10,tm,'fontsize',14)
  
  tm1 = 'Lifetime  t_L = ';
  tm2 = num2str(lifeTime,'%3.3f  ns ');
  tm = [tm1 tm2];
  text(0,0,tm,'fontsize',14)
  
  axis off
 
 figure(88)  % numerical results
  pos = [0.60 0.05 0.28 0.52];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  subplot(2,1,1)
  plot(r,R1,'b','linewidth',2)
  grid on
  tmS = '  ';
  tm1 = num2str(state(1),'%2.0f');
  tm2 = num2str(state(2),'%2.0f');
  tm3 = num2str(state(3),'%2.0f');
  tm4 =  'may need to adjust r_{max}';
  tm = [tmS tm1 tmS tm2 tmS tm3 tmS tmS tm4];
  title(tm)
  ylabel('\psi_R')
  set(gca,'fontsize',14)
  
  subplot(2,1,2)
  plot(r,R2,'r','linewidth',2)
  grid on
  ylabel('\psi_R')
  xlabel('radial distance  r  [m]')
  tmS = '  ';
  tm1 = num2str(state(4),'%2.0f');
  tm2 = num2str(state(5),'%2.0f');
  tm3 = num2str(state(6),'%2.0f');
  tm = [tm1 tmS tm2 tmS tm3];
  title(tm)
  set(gca,'fontsize',14)
  
  
% FUNCTIONS =========================================================
  
function integral = simpson1d(f,a,b)

% [1D] integration - Simpson's 1/3 rule
%      f function    a = lower bound    b = upper bound
%      Must have odd number of data points
%      Simpson's coefficients   1 4 2 4 ... 2 4 1

num = length(f);               % number of data points

sc = 2*ones(num,1);
sc(2:2:num-1) = 4;
sc(1) = 1; sc(num) = 1;

h = (b-a)/(num-1);

integral = (h/3) * f * sc;
end


function [ EB, Psi, r] = qp_fh( pqn, L, num, r_max)
% Solution of Schrodinger Equation for hydrogen atom

% INPUTS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Z = 1;

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
E = zeros(n,1);
while flag == 0
    E(n) = e_values(n,n);
    if E(n) > 0, flag = 1; end % if
    n = n + 1;
end  % while
%E(n-1) = [];
%n  = n-2;


% Corresponding Eigenfunctions 1, 2, ... ,n: Normalizing the wavefunction
 psi = zeros(num,n);
for cn = 1 : n
   psi(:,cn) = [0; e_funct(:,cn); 0];
   area = simpson1d((psi(:,cn).* psi(:,cn))',r_min,r_max);
   psi(:,cn) = psi(:,cn)/sqrt(area);
end % for

qn = pqn-L;

Psi = psi(:,qn)';
EB = -E(pqn-L);

end


