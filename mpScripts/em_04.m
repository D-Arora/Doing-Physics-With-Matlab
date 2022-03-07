% em_04.m

% DOT AND CROSS PRODUCT

% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

% DOWNLOAD Scripts from
%   https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%   https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% 220226  Matlab Version R2021b

clear
close all
clc

% INPUTS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  A = [-2;3;2];
  B = [3;2;1];
  C = [-2;2;-3];
 
% CALCULATIONS   =======================================================
% Magnitudes
  Ar = norm(A);
  Br = norm(B);
  Cr = norm(C);
% dot and cross products
  AdotB = dot(A,B);
  BdotA = dot(B,A);
  AcrossB = cross(A,B);
  BcrossA = cross(B,A);
% Angle between A and B
  angleAB = acos(AdotB/(Ar*Br));
  angleABdeg = rad2deg(angleAB);
% Unit vector n
  n = cross(A,B)/(Ar*Br*sin(angleAB));
% Area of parallelogram
  areaAB = AdotB;
% Triple products
  AdBxC = dot(A,cross(B,C));  
  CdAxB = dot(C,cross(A,B));
  BdCxA = dot(B,cross(C,A));
  BdAxC = dot(B,cross(A,C));
  volABC =  AdBxC;  
% Triple cross products
  AxBC  =  cross(A, cross(B,C));   % Ax(BxC)
  ABxC  =  cross(cross(A,B),C );   % (AxB)xC
  BAC_CAB = B.*dot(A,C) - C.*dot(A,B);


% OUTPUT   ==============================================================
disp('  ')
fprintf('  Vector A   %2.2f   %2.2f   %2.2f  \n  ',A)
fprintf('Vector B   %2.2f   %2.2f   %2.2f  \n  ',B)
fprintf('Vector C   %2.2f   %2.2f   %2.2f  \n  ',C)
fprintf('Ar = %2.4f    Br = %2.4f   Cr = %2.4f  \n  ',Ar,Br,Cr)
fprintf('AdotB = %2.4f    BdotA = %2.4f  \n  ',AdotB, BdotA)
fprintf('AcrossB   %2.4f   %2.4f   %2.4f  \n  ',AcrossB')
fprintf('BcrossA   %2.4f   %2.4f   %2.4f  \n  ',BcrossA')
fprintf('AB angle   %2.4f rad   %2.4f deg  \n  ',angleAB,angleABdeg)
fprintf('Unit vector n   %2.4f   %2.4f   %2.4f  \n  ',n')
fprintf('Area parallelogram AB = %2.4f  \n  ',areaAB')
fprintf('AdBxC = %2.4f  \n  ',AdBxC)
fprintf('CdAxC = %2.4f  \n  ',CdAxB)
fprintf('BdCxA = %2.4f  \n  ',BdCxA)
fprintf('BdAxC = %2.4f  \n  ',BdAxC)
fprintf('volume parallelepiped = %2.4f  \n  ',volABC)
fprintf('Ax(BxC)  %2.4f   %2.4f   %2.4f  \n  ',AxBC')
fprintf('(AxB)xC  %2.4f   %2.4f   %2.4f  \n  ',ABxC')
fprintf('BAC_CAB  %2.4f   %2.4f   %2.4f  \n  ',BAC_CAB')


% GRAPHICS   ============================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.15 0.2 0.30]);

  xP = [0 A(1)]; yP = [0 A(2)]; zP = [0 A(3)];
    plot3(xP,yP,zP,'b','linewidth',2);
  hold on
  xP = [0 B(1)]; yP = [0 B(2)]; zP = [0 B(3)];
    plot3(xP,yP,zP,'r','linewidth',2);

  xP = [0 n(1)]; yP = [0 n(2)]; zP = [0 n(3)];
    plot3(xP,yP,zP,'k','linewidth',2);

  xlabel('X'); ylabel('Y');zlabel('Z')
   title('A(blue)    B(red)   n(black)')
  grid on; box on
  axis equal
  view(-46,36)
  set(gca,'fontsize',12)


figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.3 0.15 0.2 0.3]);

  xP = [0 A(1)]; yP = [0 A(2)]; zP = [0 A(3)];
    plot3(xP,yP,zP,'b','linewidth',2)
  hold on
  xP = [0 B(1)]; yP = [0 B(2)]; zP = [0 B(3)];
    plot3(xP,yP,zP,'r','linewidth',2)

  xP = [0 C(1)]; yP = [0 C(2)]; zP = [0 C(3)];
    plot3(xP,yP,zP,'k','linewidth',2)

  xlabel('X'); ylabel('Y');zlabel('Z')
  title('A(blue)    B(red)   C(black)')
  grid on; box on
  axis equal
  view(-46,36)
  set(gca,'fontsize',12)
