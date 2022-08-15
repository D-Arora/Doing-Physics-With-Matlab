% cemCh3.m

% VECTOR ANALYSIS: INTEGRATION
% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemCh3.pdf



%% CELL 1:   Symbolic Integration
close all
clc
clear

% Indefintite integral [1D]
syms x a
% Input functiuon >>>>>
  f = 1/(a^2 + x^2)
% Integral  
  F = int(f,x)

% Definite integral [1D]
  f = exp(-x)*sin(x)/x
% Integral  
  F = int(f,x,0,inf) 
% Definite integral
  f = x*exp(-x)
% Integral
  f = int(f,x,0, inf)

%% CELL 2:   [1D] integration using Simpson's rule
   clear; close all; clc
% Input: number of grid point N (odd number),  lower limit a, upper limit b
  N = 999;
  a = 0;  b = 1;

  x = linspace(a,b,N);
% Input function >>>>>
  q = 0.3; r = 0.9; s = 6;
  f = 1./((x-q).^2+0.01) + 1./((x-r).^2 + 0.04) - s;
   
% Evaluate integral
  F = simpson1d(f,a,b);
% Output to Commmand Window
  fprintf('integral  F = %2.10e  \n',F)


%% CELL 3:   [2D] integration:  Matlab function  integral2
  clear; close all; clc
% Limits   >>>>>
  N = 999;
  xMin = -3; xMax= 3;
  yMin= -4; yMax = 4;
% Grid
 % x = linspace(xMin,xMax, N);
 % y = linspace(yMin,yMax, N);

% Compute integral
  funct = @(x,y) 1 + 2.*x.^2 + y.^2;
  F = integral2(funct,xMin,xMax,yMin,yMax);

  % Output to Commmand Window
  fprintf('integral  F = %2.10e  \n',F)

%% CELL 4:   [2D] integration:  simpson2d.m 
 
  clear; close all; clc
% Limits   >>>>>
  N = 999;
  xMin = -3; xMax= 3;
  yMin= -4; yMax = 4;

 % Grid
 x = linspace(xMin,xMax, N);
 y = linspace(yMin,yMax, N);
 [xx, yy] = meshgrid(x,y);

 f = 1 + 2.*xx.^2 + yy.^2;

 F = simpson2d(f,xMin,xMax,yMin,yMax);
 % Output to Commmand Window
  fprintf('integral  F = %2.10e  \n',F)


 
 %% CELL 5:   [2D] integration:  simpson2d.m   circle
 
  clear; close all; clc
% Limits   >>>>>
  N = 999;
  xMin = -1; xMax = 1;
  yMin = -1; yMax = 1;

 % Grid
 x = linspace(xMin,xMax, N);
 y = linspace(yMin,yMax, N);
 [xx, yy] = meshgrid(x,y);

 f = ones(N,N);
 f((xx.^2 + yy.^2) > 1) = 0;

 F = simpson2d(f,xMin,xMax,yMin,yMax);
 % Output to Commmand Window
  fprintf('integral  F = %2.10e  \n',F)


%% CELL 6:   [3D] integration

 clear; close all; clc
 
 xmin = 0; xmax = 2;
 ymin = 0; ymax = 4;
 zmin = 0; zmax = 8;

 fun = @(x,y,z) x.^2.*y.^3.*z.^4;
 F = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax);

 fprintf('integral  F = %2.10e  \n',F)

 %% CELL 7:     [3D] integration

 clear; close all; clc
 
 xmin = -1;
 xmax = 1;
 ymin = @(x)     -sqrt(1 - x.^2);
 ymax = @(x)      sqrt(1 - x.^2);
 zmin = @(x,y)   -sqrt(1 - x.^2 - y.^2);
 zmax = @(x,y)    sqrt(1 - x.^2 - y.^2);

 fun = @(x,y,z) x.*cos(y) + x.^2.*cos(z);
 F = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax);

 fprintf('integral  F = %2.10e  \n',F)

%% CELL 8
clear; close all; clc

N = 9999;

% PATH 1
  a = 0; b = 1+4i;
  s = linspace(a,b,N);
  f = ( (s + 1)./(s-1-2*1i) );
  F1 = simpson1d(f,a,b);

  a = b; b = 2 + 4i;
  s = linspace(a,b,N);
  f = ( (s + 1)./(s-1-2*1i) );
  F2 = simpson1d(f,a,b);
  F = F1 + F2;

% Output to Commmand Window
 disp('  ')
  disp('')
  fprintf('integral F1   FR = %2.5e  FI =%2.5e        \n',real(F1), imag(F1))
  fprintf('integral F2   FR = %2.5e  FI = %2.5e       \n',real(F2), imag(F2))
  fprintf('integral F = F1 + F2 = %2.10e  FI = %2.5e  \n',real(F),imag(F))


%%   CELL 9
   clear; close all; clc
   x1 = -1; x2 = 1; N = 9999;
   x = linspace(x1,x2,N);
   dx = x(2) - x(1);
   
   yN = sqrt(1-x.^2);

   dyN = gradient(yN,dx);
   FN = sqrt(1 + dyN.^2);
   LN = simpson1d(FN,x1,x2);
   fprintf('LN = %2.5f  \n',LN)

% Graphics
   figure(1)
   yyaxis left
   xP = x; yP = yN;
   plot(xP, yP,'b','linewidth',2)
   xlim([-1.2 1.2]);
   ylim([-1.2 1.2])
   grid on; box on
   axis square
   xlabel('x'); ylabel('y')
   set(gca,'YColor','b')
   set(gca,'fontsize',14)

   yyaxis right
   yP = dyN;
   plot(xP, yP,'r','linewidth',2)
   ylim([-10 10])
   grid on
   ylabel('dy/dx)')
   set(gca,'YColor','r')

 %% CELL 10
    clear; close all; clc
    tMin = 0; tMax = 1; N = 999;
    t = linspace(tMin,tMax,N);
    dt = t(2) - t(1);
    x = sin(5*t);
    y = t.^2 - t;
    z = t;
% Vector function
    V = [x.^2.*y.*z;  x.*y;  2.*y.*z];
    r = [x; y; z];
% Smooth [3D] curve
    rDash = gradient(r,dt);
% Curve gradient
    VdotDash = dot(V,rDash);
% Line Integral
    L = abs(simpson1d(VdotDash, tMin, tMax));
% Output
    fprintf('line integral  L = %2.6e  \n',L)
% GRAPHICS 
    figure(1)
    plot3(x,y,z,'b','LineWidth',2)
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on; box on
    set(gca,'fontsize',14)


 %% CELL 11
 % Fundamental Theorm for divergence 

    clear; close all; clc;

    tic
% [3D] grid >>>>>     
    N = 299;
    x = linspace(0,1,N); y = x; z = x;
    dx = x(2) - x(1); dy = y(2) - y(1); dz = z(2) - z(1);
    [xx,yy,zz] = meshgrid(x,y,z);
% Vector function   >>>>>    
    Vxx = yy.^2;  Vyy = 2.*xx.*yy + zz.^2; Vzz = 2.*yy.*zz;
% Divergence of vector function
    divV = divergence(xx, yy, zz, Vxx, Vyy, Vzz);
% Approximate integral 
    FN  = (sum(divV,'all'))*dx*dy*dz;

% Symbolic computation of integral    
    syms x y z
    V = [y^2  2*x*y + z^2  2*y*z];
    vars = [x y z];
% Display divergence of vector function in Command Window    
    divV = divergence(V,vars)
% >>>>> Enter divergence function from output in Command Window    
    fun = @(x,y,z) 2*x + 2*y;
% Compute volume integral   
    FS = integral3(fun,0,1,0,1,0,1);
% Output
   fprintf('volume integral (numeric)   FN = %2.6f  \n',FN)
   fprintf('volume integral (symbolic)  FS = %2.6f  \n',FS)
   disp('  ')

% Compute area integrals S and Stot
    x = linspace(0,1,N); y = linspace(0,1,N); z = linspace(0,1,N);
    S = zeros(6,1); Stot = 0; Amag = 1;
for c = 1:6
 if c == 1
    xx = 0; [yy, zz] = meshgrid(y,z);
    nA = Amag.*[-1 0 0];
 end
 if c == 2
    xx = 1; [yy, zz] = meshgrid(y,z);
    nA = Amag.*[1 0 0];
 end
 if c == 3
    yy = 0; [xx, zz] = meshgrid(x,z);
    nA = Amag.*[0 -1 0];
 end
 if c == 4
    yy = 1; [xx, zz] = meshgrid(x,z);
    nA = Amag.*[0 1 0];
 end
 if c == 5
    zz = 0; [xx, yy] = meshgrid(x,y);
    nA = Amag.*[0 0 -1];
 end
 if c == 6
    zz = 1; [xx, yy] = meshgrid(x,y);
    nA = Amag.*[0 0 1];
 end
    Vx = yy.^2; Vy = 2.*xx.*yy + zz.^2; Vz = 2*yy*zz;
    VdotA = Vx.*nA(1) + Vy.*nA(2) + Vz.*nA(3);
    S(c) = simpson2d(VdotA,0,1,0,1);
    Stot = Stot + S(c);
end
   disp('surface integrals')
   fprintf('  %2.3f ',S)
   disp('  ')
   fprintf('Surface Integral Stot = %2.6f  \n ',Stot)
   disp('  ')


 %% CELL 12
 % Fundamental Theorm for divergence 

    clear; close all; clc;

    tic
% [3D] grid >>>>>     
    N = 299;
    x = linspace(0,2,N); y = x; z = x;
    dx = x(2) - x(1); dy = y(2) - y(1); dz = z(2) - z(1);
    [xx,yy,zz] = meshgrid(x,y,z);
% Vector function   >>>>>    
    Vxx = xx.*yy;  Vyy = 2.*yy.*zz; Vzz = 3.*xx.*zz;
% Divergence of vector function
    divV = divergence(xx, yy, zz, Vxx, Vyy, Vzz);
% Approximate integral 
    FN  = (sum(divV,'all'))*dx*dy*dz;

% Symbolic computation of integral    
    syms x y z
    V = [x*y 2*y*z 3*x*z];
    vars = [x y z];
% Display divergence of vector function in Command Window    
    divV = divergence(V,vars)
% Enter divergence function from output in Command Window    
    fun = @(x,y,z) 3*x + y + 2*z;
% Compute volume integral   
    FS = integral3(fun,0,2,0,2,0,2);
% Output
   fprintf('volume integral (numeric)   FN = %2.6f  \n',FN)
   fprintf('volume integral (symbolic)  FS = %2.6f  \n',FS)
   disp('  ')

% Compute area integrals S and Stot
    x = linspace(0,2,N); y = linspace(0,2,N); z = linspace(0,2,N);
    S = zeros(6,1); Stot = 0; Amag = 4;
for c = 1:6
 if c == 1
    xx = 0; [yy, zz] = meshgrid(y,z);
    nA = Amag.*[-1 0 0];
 end
 if c == 2
    xx = 2; [yy, zz] = meshgrid(y,z);
    nA = Amag.*[1 0 0];
 end
 if c == 3
    yy = 0; [xx, zz] = meshgrid(x,z);
    nA = Amag.*[0 -1 0];
 end
 if c == 4
    yy = 2; [xx, zz] = meshgrid(x,z);
    nA = Amag.*[0 1 0];
 end
 if c == 5
    zz = 0; [xx, yy] = meshgrid(x,y);
    nA = Amag.*[0 0 -1];
 end
 if c == 6
    zz = 2; [xx, yy] = meshgrid(x,y);
    nA = Amag.*[0 0 1];
 end
    Vx = xx.*yy; Vy = 2.*yy.*zz; Vz = 3.*xx.*zz;
    VdotA = Vx.*nA(1) + Vy.*nA(2) + Vz.*nA(3);
    S(c) = simpson2d(VdotA,0,1,0,1);
    Stot = Stot + S(c);
end
   disp('surface integrals')
   fprintf('  %2.3f ',S)
   disp('  ')
   fprintf('Surface Integral Stot = %2.6f  \n ',Stot)
   disp('  ')


 %% CELL 13  Fundamental theorem of curls
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

% Compute line integral  SL
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

%% CELL 14
% Compute line integral  SL
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

 % ----------------------------------------------------------    
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

