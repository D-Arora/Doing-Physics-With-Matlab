% qmS001.m

% SPIN STATES for 1/2 spin particles

% 230314


clear 
close all
clc


% Enter spinsor components >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  A(1) = 3;
  A(2) = 1;

% Check - normalize spin state
  N = sqrt( abs(A(1))^2 + abs(A(2)^2) );

% SPINOR:  spin vector
  A = (A/N)';
%   s1 = (1/5)*(1+sqrt(8)*1i)*(1);
%   s2 = (1/5)*(2*sqrt(2)*(1+1i)*(1));

% s1 = (1/sqrt(2))*(1) 
% s2 = (1/sqrt(2))*(1)

% s1 = 0
% s2 = 1

% Select Pauli matrix operator  1 --> PX / 2 --> PY / 3 --> PZ
%  PM = 1;

% Vector spinors
%  S = [s1;s2];

% [3D] orientation in space
[theta, phi, x, y, z] = Direction(A(1),A(2));


% PAULI SPIN OPERATOR MATRICES ==========================================
% PX = [0 1; 1 0];
% PY = [0 -1i; 1i 0];
% PZ = [1 0; 0 -1];
% 
% % Rotated spinsor vectors
% SX = PX * S;
% SY = PY * S;
% SZ = PZ * S;
% 
% if PM == 1, S1 = SX(1); S2 = SX(2); end
% if PM == 2, S1 = SY(1); S2 = SY(2); end
% if PM == 3, S1 = SZ(1); S2 = SZ(2); end
% 
% % [3D] orientation in space of new spinor vector
% [thetaPS, phiPS, xPS, yPS, zPS] = pauli(S1,S2);
% 
% % OUTPUT
% fprintf('Operator PM = %2.0f  \n',PM)
% fprintf('theta = %2.3f deg   phi = %2.2f deg \n', theta,phi)
% fprintf('thetaPS = %2.3f deg   phiPS = %2.2f deg \n', thetaPS,phiPS)
% fprintf('x = %2.3f   y = %2.2f   z = %2.3f \n', x,y,z)
% fprintf('xPS = %2.3f   yPS = %2.2f   zPS = %2.3f \n', xPS,yPS,zPS)

% GRAPHICS =========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.06 0.05 0.20 0.20]);
   set(gcf,'color','w')
  
   plot3([0,x],[0,y],[0,z],'b','LineWidth',2)
   hold on
   ax = gca;
   ax.Box = 'on';
   ax.BoxStyle = 'full';
  % plot3([0,xPS],[0,yPS],[0,zPS],'r','LineWidth',2)

   Hplot = plot(0,0,'ko');
   set(Hplot,'markersize',8,'markerfacecolor','k')
   grid on
   xlim([-1 1])
   ylim([-1 1])
   zlim([-1 1])

   xlabel('x')
   ylabel('y')
   zlabel('z')

   ax = gca;
   ax.BoxStyle = 'full';
   set(gca,'fontsize',12)

%%
% clc
% theta = pi
% phi   = 0
% 
% [nP, nM] = spinN(theta,phi);
% nP



%%
% PAULI SPIN OPERATOR MATRICES ==========================================
% PX = [0 1; 1 0];
% PY = [0 -1i; 1i 0];
% PZ = [1 0; 0 -1];
% 
% zP = [1; 0];
% zM = [0; 1];
% 
% PP = (PX + 1i*PY)/2;
% PM = (PX - 1i*PY)/2;
% 
% PP*zP
% PP*zM
% PM*zP
% PM*zM
% 
% PX*PX
% PY*PY
% PZ*PZ
% 
% PX*PY
% 1i*PZ
% 
% PY*PZ
% 1i*PX
% 
% PZ*PX
% 1i*PY
% 
% % state s
% s = [3; 5]
% N = sqrt(s(1)^2+s(2)^2)
% s = s/N
% 
% zP'*s
% zM'*s

%%

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

function [nP, nM] = spinN(theta,phi)
 %  nP = [ cos(theta/2)*exp(-1i*phi/2);  sin(theta/2)*exp(-1i*phi/2)];
   nP = [ cos(theta/2)*exp(1i*phi);  sin(theta/2)];

   nM = [ sin(theta/2)*exp(-1i*phi/2); -cos(theta/2)*exp(1i*phi/2)];
   
end
