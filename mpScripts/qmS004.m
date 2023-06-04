% qmS004.m
% 230320

clear
close all
clc

% Enter ploar and arimuthal angles  [deg]
  theta = 0;
  phi = 0;
  
  thetaR = deg2rad(theta); phiR = deg2rad(phi);

% Spin operator matrices - theta / phi direction
  P(1,1) = cos(thetaR);
  P(1,2) = exp(-1i*phiR)*sin(thetaR);
  P(2,1) = exp(1i*phiR)*sin(thetaR);
  P(2,2) = -cos(thetaR);

% PAULI SPIN OPERATOR MATRICES ==========================================
  Px = [0 1; 1 0];
  Py = [0 -1i; 1i 0];
  Pz = [1 0; 0 -1]; 
  

% Spinors: up U and down D
 % U(1) = cos(thetaR/2); U(2) = sin(thetaR/2)*exp(-1i*phiR);
  U(1) = cos(thetaR/2)*exp(1i*phiR);
  U(2) = sin(thetaR/2);
  
  D(1) = sin(thetaR/2);
  D(2) = -cos(thetaR/2)*exp(-1i*phiR);

  U = [U(1);U(2)];
  D = [D(1);D(2)];

% Normalized  and Orthogonal
  UU = U'*U;
  DD = D'*D;
  UD = U'*D;

% CHECK theta and phi   [3d] orientation of spin vector
  [thetaU, phiU, xU, yU, zU] = Direction(U(1),U(2));
  [thetaD, phiD, xD, yD, zD] = Direction(D(1),D(2));

  x = sin(thetaR)*cos(phiR); y = sin(thetaR)*sin(phiR); z = cos(thetaR);

% Operator acting on spinor
  PU = P*U;
  PD = P*D;

% General spin state and predicting probabilites
  a(1) = 3;  
  a(2) = 5;
  NA = sqrt( abs(a(1))^2 + abs(a(2)^2) );
  a = a./NA;
  A = a(1)*U + a(2)*D;

% Two ways of calculating probabiities of up state and down state  
  probU = a(1)'*a(1);
  probD = a(2)'*a(2);
  probAU = (U'*A)'*(U'*A);
  probAD = (D'*A)'*(D'*A);


% Expectation values" probabilities / sandwich rule
  Pavg = (probU - probD)/2;
  Pmean = A'*P*A;
% OUTPUTS
   fprintf('  theta = %2.3f deg    phi = %2.3f deg \n',theta, phi)
   fprintf('thetaPU = %2.3f deg      phiPU = %2.3f deg \n',thetaU, phiU)
   fprintf('thetaPD = %2.3f deg      phiPD = %2.3f deg \n',thetaD, phiD)
   disp('  ')
   fprintf('U(1) = %2.3f%+.3fi    U(2) = %2.3f%+.3fi \n',real(U(1)),imag(U(1)),real(U(2)),imag(U(2)))
   fprintf('D(1) = %2.3f%+.3fi    D(2) = %2.3f%+.3fi \n',real(D(1)),imag(D(1)),real(D(2)),imag(D(2)))
   disp('  ')
   disp(' Pauli Matrix  '); disp(P)
   disp('  ')
   disp('Normalized and Orthogonal')
   fprintf('U*U = %2.2f   D*D = %2.2f   U*D = %2.2f   \n',UU,DD,UD)
   disp('  ')
   disp('Operator acting on spinor')
   fprintf('P*U = %2.2f   P*D = %2.2f   \n',PU,PD)
   disp('  ')
   disp('[3D] orientation spin vector')
   fprintf('x  = %2.2f    y = %2.2f    z = %2.2f   \n',x,y,z)
   fprintf('xU = %2.2f   yU = %2.2f   zU = %2.2f   \n',xU,yU,zU)
   fprintf('xD = %2.2f   yD = %2.2f   zD = %2.2f   \n',xD,yD,zD)
   disp('  ')
   disp('Calculating Probabilies')
   fprintf('Coeff: probU = %2.2f    probD = %2.2f   \n',probU,probD)
   fprintf('Inner products: probAU = %2.2f    probAD = %2.2f   \n',probAU,probAD)
   disp('  ')
   disp('Expectation values')
   fprintf('prob: <P> = %2.2f   \n',Pavg)
   fprintf('matrices: <P> = %2.2f   \n',Pmean)


   

% GRAPHICS ========================================================   
   figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.06 0.05 0.25 0.35]);
   set(gcf,'color','w')
  
   plot3([0,xU],[0,yU],[0,zU],'b','LineWidth',2)
   hold on
   plot3([0,xD],[0,yD],[0,zD],'r','LineWidth',2)
   
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
   txt = sprintf('\\theta = %2.1f^o   \\phi = %2.1f^o deg \n', theta,phi);
   title(txt)

   legend('up','down','Orientation','horizontal','Location','southoutside','box','off')

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
