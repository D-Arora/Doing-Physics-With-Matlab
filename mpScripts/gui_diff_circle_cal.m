% gui_diff_circle_cal

%function gui_diff_circle_cal()

wL = boxA * 1e-9;
aQ = boxB * 1e-3;
zP = boxC * 1e-3;
xPmax = boxD * 1e-3;


% =======================================================================
% INPUTS 
% =======================================================================
n1 = 100;             % Aperture: grid points for inner ring
n2 = 200;            % Aperture: grid points for outer ring
nR = 201;             % Aperture: number of rings    must be ODD
nP = 109;             % Screen (Observation plane XY): must be odd
%wL = 632.8e-9;        % wavelength [m]
%aQ = 1*1e-4;             % radius of circular aperture  [m]

uQmax = 1e-3;         % Aperture: incident energy density  [W.m^-2]

%zP = 1;               % Aperture to Screen distance  [m]
%zP = 6.5*a;
yP = 0;               % Observation point P
%xPmax = 2e-2;         % Max radial distance from Z axis



% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

EQmax = sqrt(2*uQmax/(cL*nRI*eps0));   % Aperture: Electric Field {V/m]

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

d_RL = 4*aQ^2/wL;          % Rayleigh distance  

% Aperture Space -------------------------------------------------------
zQ = 0;
A = zeros(nR,1);     % intgeral for each ring in aperture space
n = zeros(nR,1);     % number of points for ring in aperture space

UQ = uQmax * pi * aQ^2;   % Theoretical energy emitted from aperture  [J/s]


% Ring structure
%   radius of ring r  [m]     no. of data points around a ring n
%   Greater the circumference of a ring -->  more grid points
%   Width of each ring  dr
%   Total no. grid points for Aperture  nQ
rMax = aQ;
rMin = eps;
r = linspace(rMin, rMax, nR);
dr = r(2)-r(1);
m = (n2-n1) / (nR-1);
b = n2 - m * nR;

for c = 1 : nR
   n(c) = 2*round(0.5*(m * c + b))+1;
end
nQ = sum(n);

% energy emitted from aperture [J/s]
  
UQring = zeros(1,nR);
  for c = 1 : nR
     f = uQmax .* ones(1,n(c)); 
     UQring(c) = r(c) * dr * simpson1d(f,0,2*pi);
  end
  UQtheory = sum(UQring);
  
  
% Observation Space -----------------------------------------------------
xP = linspace(0,xPmax,nP);
dxP = xP(2)-xP(1);
% optical coordinates
vP = (2*pi*aQ/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
EP = zeros(1,nP);


% =======================================================================
% Computing the Rayleigh-Sommerfeld 1 Diffraction Integral
%    EP along the X axis   [V/m]
%    uP energy density (irradiance) along the X axis  [W/m^2]
%    UP energy from aperture to screen  [J/s]

for cP = 1 : nP
   for cQ = 1 : nR
      t = linspace(0,2*pi,n(cQ)); 
      xQ = r(cQ) .* cos(t);
      yQ = r(cQ) .* sin(t);
      unit = ones(1,n(cQ));
      rPQ = fn_distancePQ(xP(cP),yP,zP,xQ,yQ,zQ);
      rPQ3 = rPQ .* rPQ .* rPQ;
      kk = ik .* rPQ;
      MP1 = exp(kk);
      MP1 = MP1 ./ rPQ3;
      MP2 = zP .* (ik .* rPQ - unit);
      
      f = EQmax .* MP1 .* MP2;
      
      A(cQ) = r(cQ) * simpson1d(f,0,2*pi);
   end
     EP(cP) = dr * sum(A)/ (2*pi);
end

uP = (cL*nRI*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));
uPdB = 10 .* log10(uP./uPmax);

% Energy received - XY Observation plane  -------------------------------
%   AP      energy within a ring  [J/s]
%   UPsum   energy within XY Observation plane  [J/s]
%   UPsum   energy enclosed within rings of increasing radius  
AP    = zeros(1,nP); 
UPsum = zeros(1,nP);
for cP = 1 : nP
    AP(cP) = 2 * pi * xP(cP) * uP(cP) * dxP ;
end
UP = sum(AP);

for cP = 1 : nP-1
 UPsum(cP+1) = UPsum(cP) + AP(cP);
end

% ========================================================================
% GRAPHICS 
% ========================================================================
graphColor = ColorCode(wL);
% plot 1

x = xP .* 1000;
y = uP;
pos = [0.05 0.45 0.4 0.28];
plot1 = subplot('Position',pos);
%plot(x,y,'Color',graphColor,'lineWidth',2);
h_area = area(x,y,'EdgeColor',graphColor,'FaceColor',graphColor);
%plot(x,y,'k','lineWidth',2);
tx = 'radial position  r  (mm)';
ty = 'irradiance  I  (W/m^2)';
xlabel(tx);
ylabel(ty);

% plot 2
y = uPdB;
pos = [0.55 0.45 0.4 0.28];
plot2 = subplot('Position',pos);
plot(x,y,'Color',graphColor,'lineWidth',2);
ty = 'irradiance  I  (dB)';
xlabel(tx);
ylabel(ty);

% plot 3
y = 100 * UPsum ./ UQ;
pos = [0.1 0.08 0.25 0.25];
plot3 = subplot('Position',pos);
plot(x,y,'k','lineWidth',2);
ty = '% Energy enclosed with in circle';
xlabel(tx);
ylabel(ty);
grid on

% [3D] plots --------------------------------------------------------
  r = xP;
  t = linspace(0,2*pi,nP);

  xx = zeros(nP,nP);
  yy = zeros(nP,nP);

  for c = 1: nP
    xx(c,:) = r .* cos(t(c));
    yy(c,:) = r .* sin(t(c));
  end

  uPxy = meshgrid(uP,uP);


% plot 4
pos = [0.4 0.08 0.30 0.30];
plot4 = subplot('Position',pos);
pcolor(xx,yy,10.*(uPxy).^0.3);
shading interp
axis equal
shadingMap = gray(64);
colormap(shadingMap(:,1)*graphColor);
%colormap(gray)
axis off

% plot 5
pos = [0.7 0.08 0.30 0.30];
plot5 = subplot('Position',pos);
surf(xx,yy,10.*(uPxy).^0.3);
shading interp
%colormap(jet)
axis off
view(-38,64);
