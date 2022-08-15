% cemStokes.m

% VECTOR ANALYSIS
%    Stokes TheoremFundamental Theorm for Curls
% Ian Cooper
%   email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemCh3.pdf

% COMPUTE SURFACE INTEGRAL SA -------------------------------------------
   clear; close all; clc
    tic
 % [3D] Grid
   N  = 99;
   x = zeros(1,N); y = linspace(0,1,N); z = y;
   dx = 0; dy = y(2) - y(1); dz = z(2) - z(1);
   [yy, zz] = meshgrid(y,z);
   xx = zeros(N,N);

% Compute curl components: enter vector function V(Vx Vy Vz)
  syms x y z 
   Vx = 0;
   Vy = 2*x*z + 3*y^2;
   Vz = 4*y*z^2;
   disp('Curl grad(V x n)')
   Dx = diff(Vz,y) - diff(Vy,z)     % x component of curl
   Dy = diff(Vx,z) - diff(Vz,x)     % y component of curl
   Dz = diff(Vy,x) - diff(Vx,y)     % z component of curl
   
   curlx = double(subs(Dx,{x,y,z},{0,yy,zz}));
   curly = double(subs(Dy,{x,y,z},{0,yy,zz}));
   curlz = double(subs(Dz,{x,y,z},{0,yy,zz}));

  nA = [1 0 0]; % surface unit vector
  curlVdotA = curlx.*nA(1) + curly.*nA(2) + curlz.*nA(3);
  SA = simpson2d(curlVdotA,0,1,0,1);   % surface integral
  fprintf('Surface integral SA = %2.6f \n ',SA)
  disp('  ')

% COMPUTE LINE INTEGRAL SL ----------------------------------------
% Enter function and parameters for each path
  clear VL x y  z
  x = zeros(1,N); y = linspace(0,1,N); z = y;
  SL = zeros(1,4);
for c = 1:4
  if c == 1   % y 0 --> 1
     xL = 0; yL = y; zL = 0;
     VL = 2.*xL.*zL + 3.*yL.^2;
     SL(c) = simpson1d(VL,0,1);
  end
  if c == 2  % z 0 --> 1
     xL = 0; yL = 1; zL = z;
     VL = 4*yL*zL.^2;
     SL(c) = simpson1d(VL,0,1);
  end
  if c == 3  % y 1 --> 0
     xL = 0; yL = y; zL = 1;
     VL = 2.*xL.*zL + 3.*yL.^2;
     SL(c) = -simpson1d(VL,0,1);
  end
  if c == 4  % z 1 --> 0
     xL = 0; yL = 0; zL = z;
     VL = 4*yL*zL.^2;
     SL(c) = -simpson1d(VL,0,1);
  end
end
% Output resultsto Command Window
     SLtot = sum(SL);
     disp('   ')
     disp('Line integral: individual paths')
     fprintf('   %2.6f  ',SL)
     disp('  ')
     fprintf('Line integral  SLtot = %2.6f \n ',SLtot)
     disp('  ')


%% CELL 2
% Compute line integral  SL  -------------------------------------
  clear; close all; clc
  tic
% Enter function and parameters for each path
  N = 99;
  x = zeros(1,N); y = linspace(0,2,N); z = linspace(0,3,N);
  SL = zeros(1,3);
for c = 1:3
  if c == 1   % y 0 --> 2
     xL = 0; yL = y; zL = 0;
     SL(c) = 0;
  end
  if c == 2  % z 3 --> 0
     xL = 0; yL = 1; zL = z;
     SL(c) = 0;
  end
  if c == 3  % y 2 --> 0
     xL = 0; yL = y; zL = -(3/2).*yL + 3;
     VL = 2.*yL.*zL;
     SL(c) = -simpson1d(VL,0,max(yL));
  end
end
SLtot = sum(SL);
     disp('   ')
     disp('Line integral: individual paths')
     fprintf('   %2.6f  ',SL)
     disp('  ')
     fprintf('Line integral  SLtot = %2.6f \n ',SLtot)
     disp('  ')

 % COMPUE SURFACE INTEGRAL SA -------------------------------------------    
  % [3D] Grid
   clear; close all
   N  = 99;
   x = zeros(1,N); y = linspace(0,2,N); z = linspace(0,3,N);
   dx = 0; dy = y(2) - y(1); dz = z(2) - z(1);
   [yy, zz] = meshgrid(y,z);
   xx = zeros(N,N);

% Compute curl components: enter vector function V(Vx Vy Vz)
  syms x y z 
   Vx = x*y;
   Vy = 2*y*z;
   Vz = 3*x*z;
   disp('Curl grad(V x n)')
   Dx = diff(Vz,y) - diff(Vy,z)     % x component of curl
   Dy = diff(Vx,z) - diff(Vz,x)     % y component of curl
   Dz = diff(Vy,x) - diff(Vx,y)     % z component of curl
   
   curlx = double(subs(Dx,{x,y,z},{0,yy,zz}));
   curly = double(subs(Dy,{x,y,z},{0,yy,zz}));
   curlz = double(subs(Dz,{x,y,z},{0,yy,zz}));

  nA = [1 0 0]; % surface unit vector
  curlVdotA = curlx.*nA(1) + curly.*nA(2) + curlz.*nA(3);

  AM = ones(N,N);
  for cy = 1:N
  for cz = 1:N
      if zz(cy,cz) > -(3/2)*yy(cy,cz)+3; AM(cy,cz) = 0; end
  end
  end

  M = AM.*curlVdotA;
  
  SA = simpson2d(M,0,2,0,3);   % surface integral
  fprintf('Surface integral SA = %2.6f \n ',SA)
  disp('  ')    

    toc
