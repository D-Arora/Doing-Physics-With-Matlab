% ns_izh2D.m

clear
close all
clc


global   a b c d vpeak I uF vF dt


% CONSTANTS
  a = -0.02; b = -1; c = -60; d = 8; vpeak = 30;   I = 80;
% a = 0.2; b = 2; c = -56; d = -16; vpeak = 30;     I= -108; 
 
 D11 = 0.1; D22 = 0.1;
 
% SETUP 
% Grid
  xMax = 1;
  N = 100;
  xG = linspace(0,xMax,N);
  dx = xG(2) - xG(1);
  [x, y] = meshgrid(xG,xG);

% Time
  dt = 0.5*dx^2/D11;
  NT = 40;
  
% Constants for 2nd derivates   constants < 1 for stability
  uF = (dt/(dx)^2)*D11;
  vF = (dt/(dx)^2)*D22;

% Potentials
% Current values / Next values / Initial values / boundary values
  u1 = -65 + rand(N,N);
  v1 = -65*b + rand(N,N);
  u2 = zeros(N,N); v2 = zeros(N,N);

% Boundary Conditions  
  u1(1,1) = u1(2,2);
  u1(1,N) = u1(2,N-1);
  u1(N,1) = u1(N-1,2);
  u1(N,N) = u1(N-1,N-1);
    
  v1(1,1) = v1(2,2);
  v1(1,N) = v1(2,N-1);
  v1(N,1) = v1(N-1,2);
  v1(N,N) = v1(N-1,N-1);
  

% Update u and v variables

figure(1)

for k = 1:NT
  u2 = U(N, u1, v1, dt, uF, I);
  u1 = u2;
  
  v2 = V(N, u1, v1, dt, vF, a, b);  
  v1 = v2;
  
  for c1 = 2:N-1
      u1(1,c1) = u1(2,c1);
      u1(c1,1) = u1(c1,2);
      u1(N,c1) = u1(N-1,c1);
      u1(c1,N) = u1(c1,N-1);
      v1(1,c1) = v1(2,c1);
      v1(c1,1) = v1(c1,2);
      v1(N,c1) = v1(N-1,c1);
      v1(c1,N) = v1(c1,N-1);
  end
  
  for c1 = 1:N
  for c2 = 1: N
  if u1(c1,c2) > 30
      u1(c1,c2) = c;
      v1(c1,c2) = c + d;
  end
  end
  end
    
 % surf(x,y,u1)
 % zlim([-100 31])
  pcolor(x,y,u1)
  shading interp
  colorbar

  txt = sprintf('t = %2.0f   max %2.1f   min  %2.1e', k, max(max(u1)), min(min(u1)));
  title(txt)
  box on
  pause(0.5)
end  


  % FUNCTIONS  ========================================================
  function u2 = U(N, u1, v1, dt, uF, I)
   u2 = zeros(N,N);
   x = 2:N-1;
   y = 2:N-1;

   u2(x,y) = u1(x,y) + dt*(0.04*u1(x,y)^2 + 5*u1(x,y) + 140 - v1(x,y) + I);
   u2(x,y) = u2(x,y) + uF*(u1(x+1,y) - 2*u1(x,y) + u1(x-1,y));
   u2(x,y) = u2(x,y) + uF*(u1(x,y+1) - 2*u1(x,y) + u1(x,y-1));
  end

function v2 = V(N, u1, v1,dt, vF, a, b)
  v2 = zeros(N,N);
  x = 2:N-1;
  y = 2:N-1;

  v2(x,y) = v1(x,y) + dt*a*(b*u1(x,y) - v1(x,y));
  v2(x,y) = v2(x,y) + vF*(v1(x+1,y) - 2*v1(x,y) + v1(x-1,y));
  v2(x,y) = v2(x,y) + vF*(v1(x,y+1) - 2*v1(x,y) + v1(x,y-1));
end