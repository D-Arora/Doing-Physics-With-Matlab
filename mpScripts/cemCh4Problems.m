% cemCh4problems.m

% VECTOR ANALYSIS: PROBLEMS AND SOLUTIONS
% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220806 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%  https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemCh4.pdf


%% PROBLEM 1
close all
clear all
clc

A = [1 1]
B = [-1 1]
Amag = norm(A)
Bmag = norm(B)
AdotB = dot(A,B)
theta = acosd(AdotB/(Amag*Bmag))

C = A-B
Cmag = norm(C)
CdotC = dot(C,C)
LHS = Cmag^2
RHS = Amag^2 + Bmag^2 - 2*Amag*Bmag*cosd(theta)


%%   PROBLEM 2
close all
clear all
clc
A = [1 1 1]
B = [2 -1 1]
C = [-2 1 1]
D = [2 2 0]
AB_C = cross(cross(A,B),C)
A_BC = cross(A,cross(B,C))
AB_D = cross(cross(A,B),D)
A_BD = cross(A,cross(B,D))

%% PROBLEM 3
close
clear all
clc
A = [4 3 0]      % base diagonal
B = [0 4 1]      % side 1 diagonal
C = [3 0 1]      % side 2 diagonal

Amag = norm(A)
Bmag = norm(B)
Cmag = norm(C)

AdotB = dot(A,B)
AdotC = dot(A,C)

theta1 = acosd( AdotB/(Amag*Bmag) )
theta2 = acosd( AdotC/(Amag*Cmag) )

%% PROBLEM 4
clear
close
clc
S = [3 9 8]
F = [5 7 9]
R = F - S
Rmag = norm(R)
Rhat = R./Rmag

%% PROBLEM 5
clear ll
clear
clc

syms r x y z drdx
r = sqrt(x^2 + y^2 + z^2)
drdx = diff(r,x)
drdy = diff(r,y)
drdz = diff(r,z)

%% PROBLEM 6
clear 
clear
clc

syms f x y z drdx e
%f = x^3 + y^4 + z^5
%f = x*y^2*z^3
f = e^x*sin(y)*log(z)
drdx = diff(f,x)
drdy = diff(f,y)
drdz = diff(f,z)

%% PROBLEM 7
clear 
clear
clc

syms r1 r2 x y z 
r2 = x^2 + y^2 + z^2
grad2_x = diff(r2,x)
grad2_y = diff(r2,y)
grad2_z = diff(r2,z)

r1 = 1/sqrt(r2)
grad1_x = diff(r1,x)
grad1_y = diff(r1,y)
grad1_z = diff(r1,z)

%% PROBLEM 8
clear
close
clc

syms x y z 
R = [x^2*y*z^5 3*x*y^4*z^2 -2*x*z]
vars = [x y z];
divR = divergence(R,vars)
crossR = cross(R,vars)
grad2_x = diff(R(1),x,2)
grad2_y = diff(R(2),y,2)
grad2_z = diff(R(3),z,2)

A = -2*sin(x^2)*sin(4*y)*sin(3*z^3)
lapA = laplacian(A)

%% PROBLEM 9
close all
clc
clear

syms x y z
V = [-y x z]
vars = [x y z];
divV = divergence(V,vars)
crossV = cross(V,vars)

N = 201;
X = linspace(-10,10, N); Y = X; Z = X;
[xx, yy, zz] = meshgrid(X,Y,Z);
Vxx = -yy;   Vyy = xx;  Vzz = zz;

divV = divergence(xx, yy, zz, Vxx, Vyy, Vzz);
[curlVxx, curlVyy, curlVzz] = curl(xx, yy, zz, Vxx, Vyy, Vzz);


% GRAPHICS =============================================================
    minX = -10; minY = -10; maxX = 10; maxY = 10;
    dx = 1:20:N; dy = dx;   dz = 1;
 
    figure(1)
    set(gcf,'units','normalized','position',[0.05 0.2 0.3 0.4]);
    p1 = xx(dx,dy,dz); p2 = yy(dx,dy,dz); p3 = zz(dx,dy,dz); 
    p4 = Vxx(dx,dy,dz); p5 = Vyy(dx,dy,dz); p6 = Vzz(dx,dy,dz); 
    h = quiver3(p1, p2, p3, p4, p5, p6);
    set(h,'color',[0 0 1],'linewidth', 2);
    axis tight
    set(gca,'xLim',[minX, maxX]);
    set(gca,'yLim',[minY, maxY]);
    % set(gca,'zLim',[minZ, maxZ]);
    title('Vector Field V');
    xlabel('x'); ylabel('y'); zlabel('z');
    set(gca,'fontsize',14)
    view(-40,90)
    box on
    axis tight


%% PROBLEM 10
close all
clc
clear

syms x y z
V = [-x*y*z  (x+y)*z  x^3*y^5*z^6]
vars = [x y z];
divV = divergence(V,vars)

%% PROBLEM 11
close all
clc
clear

syms x y z k

V = [cos(k*x)/k (sin(k*y))/k 0]
vars = [x y z];
divV = divergence(V,vars)

lambda = 25;
k = 2*pi/lambda; num = 101;
X = linspace(0,100,num);   dx = X(2) - X(1);
Y = X; Z = X;
[xx, yy, zz] = meshgrid(X,Y,Z);
Vxx = cos(k.*xx)/k;   Vyy = sin(k.*yy)/k;  Vzz = zeros(num,num,num);
divV = divergence(xx, yy, zz, Vxx, Vyy, Vzz);
xP = xx(:,:,1); yP = yy(:,:,1); P = divV(:,:,1);

figure(1)
 set(gcf,'units','normalized','position',[0.05 0.2 0.3 0.4]);
 pcolor(xP,yP,P)
 shading("interp")
 colorbar
 hold on

 xQ = 2:17:100; yQ = 2:21:100;

[xxQ, yyQ] = meshgrid(xQ,yQ);
Vx = cos(k*xxQ)./k; Vy = sin(k*yyQ)./k;

h = quiver(xQ,yQ,Vx,Vy,'r','linewidth',2);
set(h,'AutoScale','on', 'AutoScaleFactor',0.6)


%% PROBLEM 12
   clear
   clc
   close all

% Limits >>>>>
  xMin = 1;
  xMax = 0;

% X range
  num = 999;
  x = linspace(xMin,xMax,num);

% Function >>>>>
  F = x.*(-2.*x + 2).^2 -x -4;

% Value of integral
  S = simpson1d(F,xMin,xMax);
  fprintf('Integral S = %2.4f  \n',S)



