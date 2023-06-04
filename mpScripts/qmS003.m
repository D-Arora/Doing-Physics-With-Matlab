% qmS001.m

% SPIN STATES for 1/2 spin particles

% 230314


clear 
close all
clc


% INPUTS Enter spinsor components >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  A(1) = 5+2i;
  A(2) = 3-6i;
 
% SPIN STATES: normalized spinors
  NA = sqrt( abs(A(1))^2 + abs(A(2)^2) );
  A = [A(1); A(2)];
  A = A./NA;
 
% [3D] orientation in space spinor A
  [theta, phi, x, y, z] = Direction(A(1),A(2));

% PAULI SPIN OPERATOR MATRICES ==========================================
  Px = [0 1; 1 0];
  Py = [0 -1i; 1i 0];
  Pz = [1 0; 0 -1]; 

% New Spinors  
  Bx = Px*A;
  By = Py*A;
  Bz = Pz*A;

% Commutation relationships
  PX = (Py*Pz - Pz*Py)/2i;
  PY = (Pz*Px - Px*Pz)/2i;
  PZ = (Px*Py - Py*Px)/2i;

% [3D] spin orientation

 [thetaBx, phiBx, xBx, yBx, zBx] = Direction(Bx(1),Bx(2));
 [thetaBy, phiBy, xBy, yBy, zBy] = Direction(By(1),By(2));
 [thetaBz, phiBz, xBz, yBz, zBz] = Direction(Bz(1),Bz(2));


% OUTPUT =============================================================
 fprintf('A(1) = %2.3f%+.3fi    A(2) = %2.3f%+.3fi \n',real(A(1)),imag(A(1)),real(A(2)),imag(A(2)))
% fprintf('A(1) = %2.3f %2.3fi   A(2) = %2.3f  %2.3fi   \n',real(A(1)),imag(A(1)), real(A(2)), imag(A(2)))
 fprintf('theta = %2.3f deg   phi = %2.3f deg \n', theta,phi)
 fprintf('x = %2.3f   y = %2.3f   z = %2.3f \n', x,y,z)
 disp('  ')
 fprintf('Bx(1) = %2.3f%+.3fi    Bx(2) = %2.3f%+.3fi \n',real(Bx(1)),imag(Bx(1)),real(Bx(2)),imag(Bx(2)))
 fprintf('thetaBx = %2.3f deg   phiBx = %2.3f deg \n', thetaBx,phiBx)
 fprintf('xBx = %2.3f   yBx = %2.3f   zBx = %2.3f \n', xBx,yBx,zBx)
 disp('  ')
 fprintf('By(1) = %2.3f%+.3fi    By(2) = %2.3f%+.3fi \n',real(By(1)),imag(By(1)),real(By(2)),imag(By(2)))
 fprintf('thetaBy = %2.3f deg   phiBy = %2.3f deg \n', thetaBy,phiBy)
 fprintf('xBy = %2.3f   yBy = %2.3f   zBy = %2.3f \n', xBy,yBy,zBy)
 disp(' ')
 fprintf('Bz(1) = %2.3f%+.3fi    Bz(2) = %2.3f%+.3fi \n',real(Bz(1)),imag(Bz(1)),real(Bz(2)),imag(Bz(2)))
 fprintf('thetaBz = %2.3f deg   phiBz = %2.3f deg \n', thetaBz,phiBz)
 fprintf('xBz = %2.3f   yBz = %2.3f   zBz = %2.3f \n', xBz,yBz,zBz)
 disp('  ')
 disp('COMMUTATION RELATIONSHIPS  ')
 disp('PX'); disp(PX)
 disp('PY'); disp(PY)
 disp('PZ'); disp(PZ)
 


% GRAPHICS =========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.06 0.05 0.15 0.5]);
   set(gcf,'color','w')

subplot(3,1,1)
   plot3([0,x],[0,y],[0,z],'b','LineWidth',2)
   hold on
   plot3([0,xBx],[0,yBx],[0,zBx],'r','LineWidth',2)
   
   ax = gca;
   ax.Box = 'on';
   ax.BoxStyle = 'full';
  
   Hplot = plot(0,0,'ko');
   set(Hplot,'markersize',8,'markerfacecolor','k')
   grid on
   xlim([-1 1])
   ylim([-1 1])
   zlim([-1 1])
   view(-41,48)

   xlabel('x')
   ylabel('y')
   zlabel('z')
   title('Bx = Px A')

   ax = gca;
   ax.BoxStyle = 'full';
   set(gca,'fontsize',12)

subplot(3,1,2)
   plot3([0,x],[0,y],[0,z],'b','LineWidth',2)
   hold on
   plot3([0,xBy],[0,yBy],[0,zBy],'r','LineWidth',2)
   
   ax = gca;
   ax.Box = 'on';
   ax.BoxStyle = 'full';
  
   Hplot = plot(0,0,'ko');
   set(Hplot,'markersize',8,'markerfacecolor','k')
   grid on
   xlim([-1 1])
   ylim([-1 1])
   zlim([-1 1])
   view(-41,48)

   xlabel('x')
   ylabel('y')
   zlabel('z')
   title('By = Py A')

   ax = gca;
   ax.BoxStyle = 'full';
   set(gca,'fontsize',12)

subplot(3,1,3)
   plot3([0,x],[0,y],[0,z],'b','LineWidth',2)
   hold on
   plot3([0,xBz],[0,yBz],[0,zBz],'r','LineWidth',2)
   
   ax = gca;
   ax.Box = 'on';
   ax.BoxStyle = 'full';
  
   Hplot = plot(0,0,'ko');
   set(Hplot,'markersize',8,'markerfacecolor','k')
   grid on
   xlim([-1 1])
   ylim([-1 1])
   zlim([-1 1])
   view(-41,48)

   xlabel('x')
   ylabel('y')
   zlabel('z')
   title('Bz = Pz A')

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
