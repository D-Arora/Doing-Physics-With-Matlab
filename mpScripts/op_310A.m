% op_301.m


clear 
close all
clc

tic

% INPUT SECTION =======================================================

% Wavelegnth  [m]
wL = 550e-9;

% Source S1 and S2  [m]
xS(1) = 2.0e-4;   yS(1) = 0;
xS(2) = -2.0e-4;  yS(2) = 0;
zS(1) = -1;    zS(2) = -1;
ES(1) = 1;  ES(2) = 1;

% Aperture
% Number of rings (must be an ODD number)
nQR = 55;  
% Additional partitons for each ring
nQA = 13 * 4 + 1;
% Radius  [m] / Position  [m]
a = 1e-3;
zQ = 0;

% Observation space  P
 nX = 55;
 nY = 55;
 nXY = 55;
 xP_max = 6e-4; 
 yP_max = 6e-4;
 
 zP = 1;
%     if a > 2e-3, rPmax = 1e-5; end   
%     if a > 4e-3, rPmax = 0.5e-5; end
    

% SETUP  ============================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

NF = a^2 / (wL * zS(1));       % FRESNEL NUMBER
NA = a / sqrt(a^2 + zS(1)^2);  % NUMERICAL APERTURE 


% Aperture Space / Ring Structure =====================================
    %   radius of ring r  [m]     no. of data points around a ring n
    %   Greater the circumference of a ring -->  more grid points
    %   Width of each ring  dr / grid points for a ring nQ
    %   Total no. grid points for Aperture  nQ_sum

    rQ = linspace(eps,a,nQR);  % radii of each ring
    drQ = rQ(2) - rQ(1);
    tQ = linspace(0,2*pi,nQA);  % angles for grid points
    xQ = zeros(nQR,nQA);  % X values for grid points
    yQ = zeros(nQR,nQA);  % Y values for grid points
    rSQ1 = zeros(nQR,nQA); % distance Q to S
    rSQ2 = zeros(nQR,nQA); % distance Q to S
    EQ1 = zeros(nQR,nQA);
    EQ2 = zeros(nQR,nQA);
    SQring = zeros(1,nQR);  % flux for a ring


  for c = 1 : nQR
    xQ(c,:) = rQ(c) .* cos(tQ);
    yQ(c,:) = rQ(c) .* sin(tQ);
  end

% distance between aperture grid points and source point  [m]
% electric field at each aperture grid point [V/m]
  for c = 1 : nQR
    rSQ1(c,:) = distancePQ(xS(1),yS(1),zS(1),xQ(c,:),yQ(c,:),zQ);
    EQ1(c,:) = ES(1) * exp(-ik .* rSQ1(c,:)) ./ rSQ1(c,:);
    rSQ2(c,:) = distancePQ(xS(2),yS(2),zS(2),xQ(c,:),yQ(c,:),zQ);
    EQ2(c,:) = ES(2) * exp(-ik .* rSQ2(c,:)) ./ rSQ2(c,:);
  end
    EQ = EQ1 + EQ2;
% % irradiance distribution
%     S_Q  = (cL*nRI*eps0/2) .* abs(EQ).^2;
%     S1_Q = (cL*nRI*eps0/2) .* abs(EQ1).^2;
%     S2_Q = (cL*nRI*eps0/2) .* abs(EQ2).^2;
% % radiant flux
%     for c = 1 : nQR
%         SQring(c) = rQ(c) * drQ * simpson1d(S_Q(c,:),0,2*pi);
%         S1Qring(c) = rQ(c) * drQ * simpson1d(S1_Q(c,:),0,2*pi);
%         S2Qring(c) = rQ(c) * drQ * simpson1d(S2_Q(c,:),0,2*pi);
%     end
%         SQ = sum(SQring);
%         S1Q = sum(S1Qring);
%         S2Q = sum(S2Qring);
  
% Observation Space / Ring Structure  ================================    
   xP = linspace(-xP_max,xP_max,nX);
   yP = linspace(-yP_max,yP_max,nY);
   
   [xxP, yyP] = meshgrid(xP,yP);
   A = zeros(nQR,1); A1 = zeros(nQR,1); A2 = zeros(nQR,1);
   EP = zeros(nX,nY); EP1 = zeros(nX,nY); EP2 = zeros(nX,nY);
  for c1 = 1 : nX
  for c2 = 1 : nY 
  for c3 = 1 : nQR 
    rPQ = distancePQ(xP(c1),yP(c2),zP,xQ(c3,:),yQ(c3,:),zQ);
    rPQ3 = rPQ .* rPQ .* rPQ;
    kk = ik .* rPQ;
    MP1 = exp(kk);
    MP1 = MP1 ./ rPQ3;
    %unit = ones(1,nQ(c3));
    MP2 = zP .* (ik .* rPQ - 1);
    f = EQ(c3,:) .* MP1 .* MP2;
    f1 = EQ1(c3,:) .* MP1 .* MP2;
    f2 = EQ2(c3,:) .* MP1 .* MP2;
    A(c3)  = rQ(c3) * simpson1d(f,0,2*pi); 
    A1(c3) = rQ(c3) * simpson1d(f1,0,2*pi); 
    A2(c3) = rQ(c3) * simpson1d(f2,0,2*pi); 
    
  end
    EP(c1,c2) = drQ * sum(A)/ (2*pi);
    EP1(c1,c2) = drQ * sum(A1)/ (2*pi);
    EP2(c1,c2) = drQ * sum(A2)/ (2*pi);
  end
  end

SP = (cL*nRI*eps0/2) .* abs(EP).^2;
SP_max = max(max(SP));
SP_dB = 10 .* log10(SP./SP_max);

SP1 = (cL*nRI*eps0/2) .* abs(EP1).^2;
SP1_max = max(max(SP1));
SP1_dB = 10 .* log10(SP1./SP1_max);

SP2 = (cL*nRI*eps0/2) .* abs(EP2).^2;
SP2_max = max(max(SP2));
SP2_dB = 10 .* log10(SP2./SP2_max);

% Angular spread
theta1 = atand(abs(xS(1)/zS(1)));
theta2 = atand(abs(xS(2)/zS(2)));
theta = abs(theta1) + abs(theta2);


% Saddle point
if abs(xS(1)) == abs(xS(2))
   saddle = SP(1,nXY(1))/SP_max;
end

if abs(xS(1)) > abs(xS(2))
    c = 1;
    while SP(c,nXY(1)) > SP(c+1,nXY(1))
     c = c+1 ;  
    end
     saddle = SP(c,nXY(1))/SP_max;
end

if abs(xS(1)) < abs(xS(2))
    c = 1;
    while SP(c,nXY(3)) > SP(c+1,nXY(3))
     c = c+1 ;  
    end
     saddle = SP(c,nXY(3))/SP_max;
end


% GRAPHICS 
% ========================================================================
col = ColorCode(wL);
fs = 14;
figure(1)
  pos = [0.1 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

%    tx = 'x_P (blue) &  y_P (red)  [\mum]';
%    ty = 'irradiance w_e  [W.m^{-2}]';
%    t1 = 'saddle point ratio =  ';
%    t2 = num2str(saddle,3);
%    tm = [t1 t2];
   x = xP;   % +X axis  
  % if xS(1) == xS(2), x =  (xP(:,nXY(1))+xS(1)).*1e6; end;
  y = SP(:, round(nY/2))./max(SP(:, round(nY/2)));
   plot(x,y,'color',col,'linewidth',2)
   
   hold on
   
  y = SP1(:, round(nY/2))./(eps+max(SP1(:, round(nY/2))));
    plot(x,y,'color',col,'linewidth',1)
   
  y = SP2(:, round(nY/2))./(eps+max(SP2(:, round(nY/2))));
    plot(x,y,'color',col,'linewidth',1)
     
   
%    x = xP(:,nXY(3)).*1e6;   % -X axis
%   % if xS(1) == xS(2), x =  (xP(:,nXY(1))+xS(1)).*1e6; end;
%    y = we_P(:, nXY(3));
%    plot(x,y,'b','linewidth',2)
% 
%    x = yP(:,nXY(2)).*1e6;   % +Y axis
%    y = we_P(:, nXY(2));
%    plot(x,y,'r','linewidth',0.5)
% 
%    x = yP(:,nXY(4)).*1e6;   % -Y axis
%    y = we_P(:, nXY(4));
%    plot(x,y,'r','linewidth',0.5)

%   xlabel(tx);
%   ylabel(ty);
%   title(tm)
  grid on
  set(gca,'fontsize',fs);

  
figure(2)
  pos = [0.4 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w')
  
  x = xP; y = yP;
  z = (SP.^0.9)';
  pcolor(x,y,z)
  myMap = zeros(64,3);
     for c1 = 1:64
      myMap(c1,:) = col;
     end
   
     shading interp
     axis equal
     shadingMap = gray(64);
     colormap(shadingMap(:,1)*col);
  
   %colorbar 
%    xlabel('x_P - x_S  [m])');
%    ylabel('y_P  [m]')
   axis equal
   grid on
   axis off
%    title('short time exposure');

toc

% FUNCTIONS =========================================================== 
function d = distancePQ(xP,yP,zP,xQ,yQ,zQ)
 % Function to calculate the distance between two points P and Q
   d = sqrt((xP-xQ).^2 + (yP-yQ).^2 + (zP-zQ)^2);
end        