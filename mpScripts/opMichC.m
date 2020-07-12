% opMichC.m

% Michelson Interferometer
%   Point source illumination of mirrors
%   Mirrors precisely aligned at right angles to the beam


% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mo/doc/op_michelson.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/


% Ian Cooper   email: matlabvisualphysics@gmail.com
% 191018   Matlab 2019b

clear
clc
close all

tic

% INPUTS  ============================================================= 

% Wavelength [m]    380 < wL < 750  nm 
  wL = 550e-9;       
 
% Detector Screen in XY plane: observation point  P(xp,yp,zp) 
% Number of grid points np x np
  np = 1559;           

% Max angle [deg]
  Amax = 5;

% Source points S1 and S2  
% Distance between sources points 1 and 2  
  S21 = 1e-3; 
 

% SETUP  ==============================================================

% angular wave number [rad/m]
  k = 2*pi/wL;      
% Number of source points
  ns = 2; 
% Axial distance from Source 1 to Detector Screen
  zSP = 0.2;
  
% Detector Screen: x and y coordinates given by variable p
  pmax = zSP*tand(Amax);
  pmin = -pmax;
  dp = (pmax-pmin)/(np-1);
  p = pmin : dp :pmax;
  
  xp = zeros(np,np);
  yp = zeros(np,np);
  zp = zSP*ones(np,np);  
  
  for c = 1 : np
   xp(:,c) = p(c);
   yp(c,:) = -p(c);
  end  %for c   

% Electric field from point sources
   E = zeros(np,np);   

% Point Sources S1 and S2: (x,y,z) coordinates / phases for reflection
  S = zeros(ns,3);     %  coordinates for source points
  S(1,:) = [0 0 0];
  S(2,:) = [0 0 S21];
  phi = zeros(ns,1);   
  phi(1) = 0;          % zero phase change mirror M1
  phi(2) = pi;         % pi   phase change mirror M2

% Loop for source points
 for cs = 1 : ns  
     
   xs = S(cs,1);   
   ys = S(cs,2);
   zs = S(cs,3);

% Distance matrix d: distance from source points S to detector points P 
  dx = xp - xs;
  dy = yp - ys;
  dz = zp - zs;

  dx = dx .*dx;
  dy = dy .*dy;
  dz = dz .*dz;

  d = dx + dy + dz;
  d = sqrt(d);

  % Electric field from Source points on Detector Screen
    E = E + exp(1i*k.*d + 1i*phi(cs))./d;

end %for cs

% Intensity on Detector Screen
  intensity = E .*conj(E);
  if max(max(intensity)) > 0
    intensity = intensity ./max(max(intensity));
  end
  
% GRAPHICS  ===========================================================  
figure(1)
 pos = [0.1 0.1 0.31 0.65];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   
subplot(2,1,1)
  xP = atand(xp./zSP); yP = atand(yp./zSP);
  pcolor(xP,yP,intensity)
  col = ColorCode(wL);
  colorMap = [linspace(0,col(1),256)',linspace(0,col(2),256)',linspace(0,col(3),256)' ];
  colormap(colorMap); 
  axis square
  shading interp
  set(gca,'xtick',-Amax:1:Amax);
  set(gca,'ytick',-Amax:1:Amax);
  xlabel('\theta_x  [ deg ]')
  ylabel('\theta_y  [ deg ]')
  set(gca,'fontsize',12)

 
subplot(2,1,2)
  n0 = ceil(np/2);
  xP = atand(p./zSP);
  yP = intensity(n0,:);
  plot(xP,yP,'color',col,'linewidth',2)

  set(gca,'xtick',-Amax:1:Amax);
  grid on
  xlabel('\theta  [ deg ]')
  ylabel('intensity [ a.u. ]')
  set(gca,'fontsize',12)
% [Imax, pPeaks] = findpeaks(yP);
% AnglePeaks = xP(pPeaks);
% 
% S1P = sqrt(p(pPeaks).^2 + p(pPeaks).^2 + zSP^2);
% S2P = sqrt(p(pPeaks).^2 + p(pPeaks).^2 + (zSP-S21)^2);
% 
% orderN = (S1P-S2P)./wL;




toc
