% op_301.m


clear 
close all
clc

tic

% INPUT SECTION =======================================================

% Wavelegnth  [m]
wL = 550e-9;

% Source S1     position  [m]   electric feild [a.u.]
xS(1) = 5e-5;
yS(1) = 0;
zS(1) = 1; 
ES(1) = 1;

% Source S2;
xS(2) = -xS(1);
yS(2) = 0;
zS(2) = zS(1);
ES(2) = ES(1);

% Aperture
% Radius a  [m] / Position zQ  [m]
a = 4e-3;
zQ = 0;
% Number of rings (must be an ODD number)
nQR = 105;  
% Additional partitons for each ring   
nQA = 36 * 4 + 1;


% Observation space  P
 zP = 22e-3;
 nX = 105;
 nY = 105;
% xP_max = 15e-6; 
% yP_max = 15e-6;
 

    

% SETUP  ============================================================
k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

% Angular spread:  angles  [rad]
  theta1 = atan(abs(xS(1)/zS(1)));
  theta2 = atan(abs(xS(2)/zS(2)));
  theta = abs(theta1) + abs(theta2);
  
  thetaMin = 0.61*wL/a;
  xMin = zP * tan(thetaMin);
  
  x1Max = zP * tan(theta1);
  x2Max = zP * tan(theta2);
  if xS(1) > 0; x1Max = -x1Max; end
  if xS(2) > 0; x2Max = -x2Max; end
  
  
  x1Min(1) = x1Max + xMin;
  x1Min(2) = x1Max - xMin;
  x2Min(1) = x2Max + xMin;
  x2Min(2) = x2Max - xMin;
  
  xP_max = 1 * (abs(max(x1Min)) + abs(max(x2Min)));
  yP_max = xP_max;

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

  
% Observation (Screen) Space   ================================    
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

% Irradiance  S  [a.u.]  
  SP = abs(EP).^2;
  SP_max = max(max(SP));
  
  SP1 = abs(EP1).^2;
  SP2 = abs(EP2).^2;
  
%   Normalized irradiance  
    SP = SP ./ SP_max;
    SP1 = SP1 ./ SP_max;
    SP2 = SP2 ./ SP_max;

%   Irradiance  [dB]
    SP_dB  = 10 .* log10(SP);
    SP1_dB = 10 .* log10(SP1);
    SP2_dB = 10 .* log10(SP2);



% Saddle point (assume saddel point is located at screen postion xP = 0)
   nX0 = find(xP == 0);      % index for xP = 0
   nY0 = find(yP == 0);      % index for yP = 0
   saddle = SP(nY0,nX0);
   

% OUTPUT DISPLAY  =====================================================

fprintf('x1_Max = %2.3e m \n',x1Max)
fprintf('x1_Min = %2.3e m   %2.3e \n',x1Min(2),x1Min(1))
fprintf('x2_Max = %2.3e m \n',x2Max)
fprintf('x2_Min = %2.3e m   %2.3e \n',x2Min(2),x2Min(1))


% GRAPHICS  ===========================================================
col = ColorCode(wL);
fs = 14;
figure(1)
  pos = [0.1 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');



   x = xP;   % +X axis  
  % if xS(1) == xS(2), x =  (xP(:,nXY(1))+xS(1)).*1e6; end;
  y = SP(:, round(nY/2))./max(SP(:, round(nY/2)));
   plot(x,y,'color',col,'linewidth',2)
   
   hold on
   
  y = SP1(:, round(nY/2))./(eps+max(SP1(:, round(nY/2))));
    plot(x,y,'color',col,'linewidth',1)
   
  y = SP2(:, round(nY/2))./(eps+max(SP2(:, round(nY/2))));
    plot(x,y,'color',col,'linewidth',1)
     
  ylim([-0.1 1.1])
  
  ylabel('irradiance  S  [ a.u.]');
  xlabel('x_{screen}  [ m ]');
  t1 = 'saddle point =  ';
  t2 = num2str(saddle,3);
  tm = [t1 t2];
  title(tm)  
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