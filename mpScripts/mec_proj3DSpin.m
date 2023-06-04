% mec_proj3D1.m

% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6044256/

% https://www.chemeurope.com/en/encyclopedia/Density_of_air.html

% https://sethna.lassp.cornell.edu/SimScience/fluid/red/golf_exp.html

% https://folk.ntnu.no/leifh/teaching/tkt4140/._main021.html



close all
clc
clear

global yLand

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Club head speed [km/h]
  vC = 204.4/1.45;
% Smash factor
  sF = 1.45;
% Launch angle [deg]
  theta = 14.1;
% Face angle [deg]
  phi = 0;
% Initial angular velocity components  [rpm]
  rpmX = 000; rpmY = 000; rpmZ = 5999;
% max simulation time  [s]
  tMax = 10;
% Initial displacement components  [m]  
  x0 = 0; y0 = 0; z0 = 0;
% Landing spot: Terminate simulation at y = yLand
  yLand = 0;
% air density  [1.204  kg/m^3] 
  rho = 1.204;              
% color for some plots
  col = [0 0 1];   

% SETUP  ======================================================  
  m = 45.93e-3;         % golf ball mass [kg];
  r = 42.67e-3/2;       % golf ball radius  [m]
  g = 9.8;              % acceleration due to gravity  [m/s^2]
  
  
  A = pi*r^2;           % golf ball cross-sectional area [m^2]
  k = 0.5*rho*A/m;      % force coeff

  % Initial ba;; velocity and components  [m/s]
  v0 = (vC/3.6)*sF;
  vx0 = v0*cosd(theta)*cosd(phi);
  vy0 = v0*sind(theta);
  vz0 = v0*cosd(theta)*sind(phi);

% ******************************************************************
% Angular velocity compoents   [rad/s]  
  wx = rpmX*pi/30;  wy = rpmY*pi/30; wz = rpmZ*pi/30;
  NS = 299;
  wz = linspace(2001,5999,NS)*pi/30;
  w = zeros(NS,1); X = w; H = w;
  for s = 1:NS
   % w = sqrt(wx^2 + wy^2 + wz(s)^2);
 
% Model parameters passed to solver for ODE 
  K = [g k wx wy wz(s)];

% ODE solver
  tSpan = linspace(0,tMax,1201);
% tSpan = [0 tMax];
% Events: check final vertical displacement to stop ode
% options = odeset('RelTol',1e-8,'AbsTol',1e-6,'Events',@landing);
 options = odeset('Events',@landing);
 u0 = [x0 vx0 y0 vy0 z0 vz0];

% Solve ODE
[t, u] = ode45(@(t,u) EqM(t,u,K), tSpan, u0, options);

% Extract displacements, velocities
  x =  u(:,1);
  vx = u(:,2);
  y =  u(:,3);
  vy = u(:,4);
  z  = u(:,5);
  vz = u(:,6);

%  Drag force  D and Lift force L 
%    N = length(x);
%    Dx = zeros(N,1); Dy = Dx; Dz = Dx; Lx = Dx; Ly = Dx; Lz = Dx;
%    w = sqrt(wx^2 + wy^2 + wz(s)^2);
% 
 %   for c = 1:N
%        v = sqrt( vx(c)^2 + vy(c)^2 + vz(c)^2 );  
%        [CD,CL] = cDcL(v,w);
%        KD = m*k*CD;  KL = m*k*CL;
%        Dx(c) = -KD*v*vx(c);
%        Dy(c) = -KD*v*vy(c);
%        Dz(c) = -KD*v*vz(c);
% 
%        Lx(c) = KL*(v/w)*( wy*vz(c) - wz(s)*vy(c) );
%        Ly(c) = KL*(v/w)*( wz(s)*vx(c) - wx*vz(c) );
%        Lz(c) = KL*(v/w)*( wx*vy(c) - wy*vx(c) );
%    end
%       W = m*g;
%       Dx = Dx./W; Dy = Dy./W; Dz = Dz./W;
%       Lx = Lx./W; Ly = Ly./W; Lz = Dz./W;
 
%  RESULTS  ==================================================
   tF = t(end);       % simulation time
   vF = sqrt(vx(end)^2 + vy(end)^2);   % final velocity
   xF = x(end);       % final X displacement
   yF = y(end);       % final y displacement
   zF = z(end);       % final z displacement
%    vxF = vx(end);     % final X velocity
%    vyF = vy(end);     % final y velocity
%    vzF = vz(end);     % final z velocity
    yMax = max(y); xMax = max(x); zMax = max(z);
%    Av = abs(atand(vyF/vxF));

   X(s) = xMax;
   H(s) = yMax;
  end

%   disp('INPUTS')
%   fprintf('vC = %3.1f kph = %3.1f m/s   sF = %3.2f   theta = %3.1f deg   phi = %3.1f deg  \n',vC,vC/3.6,sF,theta,phi)
%   fprintf('Initial velocity of ball: v0 = %3.1f kph = %3.1f m/s   vx0 = %3.1f m/s   vy0 = %3.1f m/s   vz0 = %3.1f m/s  \n',v0*3.6,v0,vx0,vy0,vz0)
%   fprintf('Spin of ball [rpm]: w0 = %3.1f   wx = %3.1f   wy = %3.1f   wz = %3.1f  \n',w*30/pi,rpmX,rpmY,rpmZ)
%   disp('RESULTS')
%   fprintf('tF = %3.1f s   xMax = %3.1f m yMax = %3.1f m  zMax = %3.1f m  \n',tF,xMax,yMax,zMax)
%   fprintf('xF = %3.1f m   yF = %3.1f m   zF = %3.1f m \n',xF,yF,zF) 
%   fprintf('vF = %3.1f m/s   vxF = %3.1f m/s   vyF = %3.1f m/s   vzF = %3.1f m/s   Av = %3.1f deg \n',vF, vxF,vyF,vzF, Av) 

% GRAPHICS  ==================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.06 0.05 0.35 0.50]);
  set(gcf,'color','w')  
 
  subplot(2,1,1)
  plot(wz*30/pi,X,'color',col,'linewidth',2)
  hold on; grid on
  txtX = 'spin rate w_z  [ rpm ]';
  txtY = 'Range  [ m ]';
  txtT = '  ';
  plotSetup(txtX,txtY,txtT)
  hold on
  subplot(2,1,2)
  plot(wz*30/pi,H,'color',col,'linewidth',2)
  hold on; grid on
  txtX = 'spin rate w_z  [ rpm ]';
  txtY = 'Max height  [ m ]';
  txtT = '  ';
  plotSetup(txtX,txtY,txtT)
  hold on

% FUNCTIONS  ==================================================

function uDot = EqM(~,u,K)
  g = K(1); k = K(2); 
  wx = K(3); wy = K(4); wz = K(5);
  uDot = zeros(6,1);
  vx = u(2); vy = u(4); vz = u(6);
  v = sqrt(vx^2 + vy^2+ vz^2);
  w = sqrt(wx^2 + wy^2+ wz^2);

  r = 42.67e-3/2;       % golf ball radius  [m]
  rpm = w*30/pi;
  R = rpm*r/v;
%  CD = 0.1403 - 0.3406*R*log(R) + 0.3747*R^1.5;
%  CL = 0.3996 + 0.1583*log(R) + 0.03790*R^(-0.5);
  

  S = r*w/v;
 
  [CD,CL] = cDcL(v,w);
 % CL = -3.25*S^2 + 1.99*S;
 % CL = 0; CD = 0;
 %CL = 0;

  uDot(1) = vx;
  uDot(2) = -k*( CD*v*vx - CL*(v/w)*(wy*vz - wz*vy) ) ;
  uDot(3) = vy;
  uDot(4) = -k*( CD*v*vy - CL*(v/w)*(wz*vx - wx*vz) ) - g ; 
  uDot(5) = vz;
  uDot(6) = -k*( CD*v*vz - CL*(v/w)*(wx*vy - wy*vz) ) ;
end

function [check,stop,direction] = landing(~,u)      % EVENTS
    global yLand
    check = (u(3) < yLand);    % Input final vetical displacement
    stop = 1;              % Stop execution
    direction = 1;         % Projectile falling +1 / rising -1
end

function plotSetup(txtX,txtY,txtT)
grid on
xlabel(txtX)
ylabel(txtY)
title(txtT,'FontWeight','normal')
set(gca,'fontsize',14)
end

 
function [cD, cL] = cDcL(v,w)

nrpm = w*30/pi;
v0 = [13.7, 21.6, 29.9, 38.4, 46.9, 55.2, 63.1, 71.9, 80.2, 88.1];

nrpm0 = linspace(2000,6000,21);

CD = zeros(21,10);

CD(1,:)  = [0.3624, 0.2885, 0.2765, 0.2529, 0.2472, 0.2481, 0.2467, 0.2470, 0.2470, 0.2470];
CD(2,:)  = [0.3806, 0.3102, 0.2853, 0.259 , 0.2507, 0.2498, 0.2485, 0.2486, 0.2484, 0.2484];
CD(3,:)  = [0.3954, 0.3288, 0.2937, 0.2649, 0.2543, 0.2516, 0.2504, 0.2502, 0.2497, 0.2497];
CD(4,:)  = [0.407 , 0.3443, 0.3018, 0.2708, 0.258 , 0.2535, 0.2522, 0.2518, 0.2511, 0.2511];
CD(5,:)  = [0.4153, 0.3566, 0.3095, 0.2765, 0.2617, 0.2556, 0.2541, 0.2534, 0.2524, 0.2524];
CD(6,:)  = [0.4203, 0.3658, 0.3169, 0.2822, 0.2655, 0.2578, 0.256 , 0.255 , 0.2538, 0.2538];
CD(7,:)  = [0.412 , 0.3719, 0.324 , 0.2876, 0.2693, 0.2602, 0.2579, 0.2566, 0.2551, 0.2551];
CD(8,:)  = [0.396 , 0.3749, 0.3308, 0.293 , 0.2732, 0.2627, 0.2599, 0.2582, 0.2565, 0.2565];
CD(9,:)  = [0.3876, 0.3766, 0.3372, 0.2983, 0.2772, 0.2653, 0.2619, 0.2598, 0.2578, 0.2578];
CD(10,:) = [0.41  , 0.3854, 0.3433, 0.3034, 0.2811, 0.2681, 0.2639, 0.2614, 0.2592, 0.2592];
CD(11,:) = [0.4288, 0.4003, 0.349 , 0.3084, 0.2852, 0.271 , 0.2659, 0.263 , 0.2605, 0.2605];
CD(12,:) = [0.4445, 0.4082, 0.3544, 0.3133, 0.2893, 0.2741, 0.268 , 0.2646, 0.2619, 0.2619];
CD(13,:) = [0.4575, 0.4153, 0.3595, 0.318 , 0.2934, 0.2772, 0.2701, 0.2662, 0.2632, 0.2632];
CD(14,:) = [0.4682, 0.4215, 0.3643, 0.3227, 0.2976, 0.2806, 0.2722, 0.2678, 0.2646, 0.2646];
CD(15,:) = [0.4772, 0.4269, 0.3687, 0.3272, 0.3019, 0.284 , 0.2743, 0.2694, 0.2659, 0.2659];
CD(16,:) = [0.4848, 0.4314, 0.3728, 0.3316, 0.3062, 0.2876, 0.2765, 0.271 , 0.2673, 0.2673];
CD(17,:) = [0.4914, 0.435 , 0.3765, 0.3358, 0.3105, 0.2913, 0.2787, 0.2726, 0.2686, 0.2686];
CD(18,:) = [0.4976, 0.4377, 0.3799, 0.34  , 0.3149, 0.2952, 0.2809, 0.2742, 0.27  , 0.27  ];
CD(19,:) = [0.5039, 0.4395, 0.383 , 0.344 , 0.3194, 0.2992, 0.2831, 0.2758, 0.2713, 0.2713];
CD(20,:) = [0.5105, 0.4405, 0.3858, 0.3479, 0.3239, 0.3034, 0.2854, 0.2774, 0.2727, 0.2727];
CD(21,:) = [0.518 , 0.4406, 0.3882, 0.3517, 0.3285, 0.3076, 0.2877, 0.279 , 0.274 , 0.274 ];

CL(1,:)  = [0.104 , 0.1846, 0.246 , 0.1984, 0.1762, 0.1538, 0.1418, 0.136 , 0.128 , 0.1276];
CL(2,:)  = [0.1936, 0.2318, 0.259 , 0.2089, 0.1824, 0.1603, 0.1476, 0.1405, 0.1321, 0.1324];
CL(3,:)  = [0.2608, 0.2694, 0.2715, 0.2191, 0.1886, 0.1668, 0.1533, 0.145 , 0.1362, 0.1367];
CL(4,:)  = [0.309 , 0.2986, 0.2835, 0.2289, 0.1949, 0.1731, 0.1589, 0.1494, 0.1403, 0.1404];
CL(5,:)  = [0.3418, 0.3205, 0.2947, 0.2384, 0.2011, 0.1794, 0.1646, 0.1539, 0.1444, 0.1436];
CL(6,:)  = [0.3624, 0.3362, 0.3049, 0.2475, 0.2073, 0.1856, 0.1702, 0.1584, 0.1485, 0.1464];
CL(7,:)  = [0.3743, 0.347 , 0.314 , 0.2562, 0.2135, 0.1916, 0.1757, 0.1629, 0.1526, 0.1488];
CL(8,:)  = [0.3808, 0.3541, 0.3217, 0.2644, 0.2197, 0.1976, 0.1813, 0.1673, 0.1567, 0.1508];
CL(9,:)  = [0.3854, 0.3584, 0.328 , 0.2722, 0.2259, 0.2035, 0.1868, 0.1718, 0.1608, 0.1524];
CL(10,:) = [0.3915, 0.3614, 0.3325, 0.2795, 0.2322, 0.2092, 0.1923, 0.1763, 0.1649, 0.1539];
CL(11,:) = [0.4005, 0.364 , 0.3347, 0.2862, 0.2384, 0.2149, 0.1977, 0.1807, 0.169 , 0.1582];
CL(12,:) = [0.401 , 0.3696, 0.338 , 0.2925, 0.2446, 0.2205, 0.2032, 0.1852, 0.1731, 0.1624];
CL(13,:) = [0.4026, 0.3748, 0.3412, 0.2982, 0.2508, 0.2259, 0.2086, 0.1897, 0.1772, 0.1666];
CL(14,:) = [0.405 , 0.3797, 0.344 , 0.3033, 0.257 , 0.2313, 0.2139, 0.1941, 0.1813, 0.1707];
CL(15,:) = [0.4084, 0.3845, 0.3468, 0.3078, 0.2632, 0.2366, 0.2193, 0.1986, 0.1854, 0.1749];
CL(16,:) = [0.413 , 0.3894, 0.3496, 0.3119, 0.2695, 0.2418, 0.2246, 0.2031, 0.1895, 0.1791];
CL(17,:) = [0.4192, 0.3947, 0.3527, 0.3149, 0.2757, 0.2469, 0.2298, 0.2076, 0.1936, 0.1833];
CL(18,:) = [0.4271, 0.4004, 0.3563, 0.3184, 0.2819, 0.2518, 0.2351, 0.212 , 0.1977, 0.1875];
CL(19,:) = [0.4371, 0.4069, 0.3607, 0.3233, 0.2881, 0.2567, 0.2403, 0.2165, 0.2018, 0.1916];
CL(20,:) = [0.4494, 0.4143, 0.366 , 0.33  , 0.2943, 0.2615, 0.2455, 0.221 , 0.2059, 0.1958];
CL(21,:) = [0.4644, 0.4227, 0.3724, 0.3393, 0.3005, 0.2662, 0.2507, 0.2254, 0.21  , 0.2   ];

%CD = CD./max(max(CD));
%CL = CL./max(max(CL));


cD = interp2(v0,nrpm0,CD,v,nrpm);
cL = interp2(v0,nrpm0,CL,v,nrpm);
end
