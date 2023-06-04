% qmS004.m
% 230320

clear
close all
clc

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Spinor |Z>  = z1|zU> + z2|zD>    |ZU> = [1;0] |ZD> = [0;1]
% Enter coefficients z1 and z2
  z1 = 1+1i;
  z2 = 1i;

% Direction angles for orthonormal basis vectors |U> and |D>
% polar (theta) and azimuthal (phi) angles [deg]
  thetaDeg = 45;
  phiDeg   = 30;

% SETUP ========================================================== 
% Spinor |Z>  = z1|zU> + z2|zD>    |ZU> = [1;0] |ZD> = [0;1]
% Normalized
  A = sqrt(z1'*z1 + z2'*z2);
  z1 = z1/A; z2 = z2/A;
  Z = [z1;z2];

% Direction matrix U an D
% Angles deg to rad
  theta = deg2rad(thetaDeg); phi = deg2rad(phiDeg);
% U and D matrices  
  U(1) = cos(theta/2)*exp(1i*phi);
  U(2) = sin(theta/2);
  D(1) = sin(theta/2);
  D(2) = -cos(theta/2)*exp(-1i*phi);
  U = [U(1);U(2)];
  D = [D(1);D(2)];

% Spinor |N>  = b1|U> + b2|D>
  n2 = (z1*U(2) - z2*U(1))/(U(2)*D(1) - U(1)*D(2));
  n1 = (z1-n2*D(1))/U(1);
  N  = n1*U + n2*D;

% CHECK theta and phi   [3d] orientation of spin vector
  [thetaU, phiU, xU, yU, zU] = Direction(U(1),U(2));
  [thetaD, phiD, xD, yD, zD] = Direction(D(1),D(2));
  [thetaN, phiN, xN, yN, zN] = Direction(N(1),N(2));

% SPIN STATE MATRICES
  P(1,1) = cos(theta);
  P(1,2) = exp(-1i*phi)*sin(theta);
  P(2,1) = exp(1i*phi)*sin(theta);
  P(2,2) = -cos(theta);

% PAULI SPIN OPERATOR MATRICES
  Px = [0 1; 1 0];
  Py = [0 -1i; 1i 0];
  Pz = [1 0; 0 -1];

% Normalized  and Orthogonal
  UU = U'*U;
  DD = D'*D;
  UD = U'*D;

% Two ways of calculating probabiities of up state and down state  
   probzU = z1'*z1;
   probzD = z2'*z2;
   probnU = n1'*n1;
   probnD = n2'*n2;

   probNU = (U'*N)'*(U'*N);
   probND = (D'*N)'*(D'*N);

% (x,y,z) components for direction vector
  x = sin(theta)*cos(phi); y = sin(theta)*sin(phi); z = cos(theta);

% Operator acting on spim vectors
%  P  Px  Py  Pz
  PU = P*U;
  PD = P*D;
  PA = P*N;
% Enter operator to act on |N> to give spinor |B>  >>>>>>>>>>>>>>>>>
  txt ='Operator  B = Pz*|N>';
  B = Pz*N;
  [thetaB, phiB, xB, yB, zB] = Direction(N(1),N(2));

%  Expectation values" probabilities / sandwich rule
   Pavg = (probnU - probnD)/2;
%   Pmean = A'*P*A;


% OUTPUTS
  fprintf('|Z>  z1 = %2.3f%+.3fi    z2  = %2.3f%+.3fi \n',real(z1),imag(z1),real(z2),imag(z2));
  disp(Z)
  fprintf('|N>  n1 = %2.3f%+.3fi    n2  = %2.3f%+.3fi \n',real(n1),imag(n1),real(n2),imag(n2));
  disp(N)
  disp('Spin in direction (theta, phi)')
  fprintf('   theta = %2.3f deg    phi = %2.3f deg \n',thetaDeg, phiDeg)
  fprintf('   U(1) = %2.3f%+.3fi    U(2) = %2.3f%+.3fi \n',real(U(1)),imag(U(1)),real(U(2)),imag(U(2)))
  fprintf('   D(1) = %2.3f%+.3fi    D(2) = %2.3f%+.3fi \n',real(D(1)),imag(D(1)),real(D(2)),imag(D(2)))
  disp('Probabilities')
  fprintf('   |Z> state:  prob(ZUp) = %2.3f   prob(ZDown) = %2.3f  \n',probzU,probzD)
  fprintf('   |N> state:  prob(NUp) = %2.3f   prob(NDown) = %2.3f  \n',probnU,probnD)
%    fprintf('thetaPU = %2.3f deg      phiPU = %2.3f deg \n',thetaU, phiU)
%    fprintf('thetaPD = %2.3f deg      phiPD = %2.3f deg \n',thetaD, phiD)
%    disp('  ')
%)
%    disp('  ')
%    disp(' Pauli Matrix  '); disp(P)
%    disp('  ')
%    disp('Normalized and Orthogonal')
%    fprintf('U*U = %2.2f   D*D = %2.2f   U*D = %2.2f   \n',UU,DD,UD)
%    disp('  ')
%    disp('Operator acting on spinor')
%    fprintf('P*U = %2.2f   P*D = %2.2f   \n',PU,PD)
%    disp('  ')

%    fprintf('xU = %2.2f   yU = %2.2f   zU = %2.2f   \n',xU,yU,zU)
%    fprintf('xD = %2.2f   yD = %2.2f   zD = %2.2f   \n',xD,yD,zD)
%    disp('  ')
%    disp('Calculating Probabilies')
%    fprintf('Coeff: probU = %2.2f    probD = %2.2f   \n',probU,probD)
%    fprintf('Inner products: probAU = %2.2f    probAD = %2.2f   \n',probAU,probAD)
%    disp('  ')
     disp('Expectation values')
     fprintf('   <N> = %2.2f   \n',Pavg)
     disp(txt)
     fprintf('   |N>   xN  = %2.2f    yN = %2.2f    zN = %2.2f   \n',xN,yN,zN)
     fprintf('   |B>   xB  = %2.2f    yB = %2.2f    zB = %2.2f   \n',xB,yB,zB)
% GRAPHICS ========================================================   
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.06 0.05 0.35 0.65]);
   set(gcf,'color','w')
  
subplot(2,1,1)  
   plot3([0,0],[0,0],[0,1],'b','LineWidth',2)
   hold on
   plot3([0,0],[0,0],[0,-1],'r','LineWidth',2)
   plot3([0,xN],[0,yN],[0,zN],'k','LineWidth',2)
   
   ax = gca;
   ax.Box = 'on';
   ax.BoxStyle = 'full';
  
   Hplot = plot(0,0,'ko');
   set(Hplot,'markersize',8,'markerfacecolor','k')
   grid on
   xlim([-1 1])
   ylim([-1 1])
   zlim([-1 1])
   view(-27,28)

   xlabel('x')
   ylabel('y')
   zlabel('z')
   txt = sprintf('|Z> = z_1|Z_U> + z_2|Z_D>   \\theta = %2.1f^o   \\phi = %2.1f^o deg \n', 0,0);
   title(txt)

   legend('up','down','|Z> = |N>','Orientation','horizontal','Location','southoutside','box','off')

   ax = gca;
   ax.BoxStyle = 'full';
   set(gca,'fontsize',12)

subplot(2,1,2)  
   plot3([0,xU],[0,yU],[0,zU],'b','LineWidth',2)
   hold on
   plot3([0,xD],[0,yD],[0,zD],'r','LineWidth',2)
   plot3([0,xN],[0,yN],[0,zN],'k','LineWidth',2)
   
   ax = gca;
   ax.Box = 'on';
   ax.BoxStyle = 'full';
  
   Hplot = plot(0,0,'ko');
   set(Hplot,'markersize',8,'markerfacecolor','k')
   grid on
   xlim([-1 1])
   ylim([-1 1])
   zlim([-1 1])
   view(-27,28)

   xlabel('x')
   ylabel('y')
   zlabel('z')
   txt = sprintf('|N> = n_1|N_U> + n_2|N_D>   \\theta = %2.1f^o   \\phi = %2.1f^o deg \n', thetaDeg,phiDeg);
   title(txt)

   legend('up','down','|N> = |Z>','Orientation','horizontal','Location','southoutside','box','off')

   ax = gca;
   ax.BoxStyle = 'full';
   set(gca,'fontsize',12)

% FUNCTIONS  ===========================================================

function [theta, phi, x, y, z] = Direction(s1,s2)
     s1Mag = abs(s1);
     theta = 2*acosd(s1Mag);
     phi1 = rad2deg(angle(s1));
     phi2 = rad2deg(angle(s2));
     phi = phi1 - phi2;
     x = sind(theta)*cosd(phi);
     y = sind(theta)*sind(phi);
     z = cosd(theta);
end
