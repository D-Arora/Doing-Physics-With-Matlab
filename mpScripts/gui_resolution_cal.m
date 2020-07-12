% gui_resolution_cal.m

wL = boxA * 1e-9;
a = boxB * 1e-3;
xS(1) = boxC * 1e-6;
xS(2) = boxD * 1e-6;
zP = boxE * 1e-3;



% =======================================================================
% INPUTS 
% =======================================================================

% Aperture space
    nQR = 55;               % number of rings - must be odd
    nQA = 13 * 4 + 1;               % additional partitons for each ring
    zQ = 0;
    
% Source  S
    zS = 22e-3; 
    yS = 0;            % always set yS = 0
    ES = 1e-3;
    
% Observation space  P
    %nP = 14 * 4 + 1;             % Screen (Observation plane XY): must be odd  N * 4 + 1
    nPR = 55;
    nPA = 13 * 4 + 1;
    rPmax = 2e-5;         
    if a > 2e-3, rPmax = 1e-5; end;   
    if a > 4e-3, rPmax = 0.5e-5; end;  
    
% ========================================================================
% SETUP 
% ========================================================================
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nRI = 1;                % refractive index

k = 2*pi/wL;            % propagation constant  [rad/s]
ik = 1i*k;              % j k

NF = a^2 / (wL * zS);       % FRESNEL NUMBER
NA = a / sqrt(a^2 + zS^2);  % NUMERICAL APERTURE

% ========================================================================
% Aperture Space / Ring Structure 
% ========================================================================
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
    %EQ = zeros(nQR,nQA);  % electric field at Q
    EQ1 = zeros(nQR,nQA);
    EQ2 = zeros(nQR,nQA);
    WE_Qring = zeros(1,nQR);  % flux for a ring


for c = 1 : nQR
    xQ(c,:) = rQ(c) .* cos(tQ);
    yQ(c,:) = rQ(c) .* sin(tQ);
   % tQ = tQ + 2*pi*c/nQR;    
end




% distance between aperture grid points and source point
% electric field at each aperture grid point [V/m]
for c = 1 : nQR
    rSQ1(c,:) = fn_distancePQ(xS(1),yS,zS,xQ(c,:),yQ(c,:),zQ);
    EQ1(c,:) = ES * exp(-ik .* rSQ1(c,:)) ./ rSQ1(c,:);
    rSQ2(c,:) = fn_distancePQ(xS(2),yS,zS,xQ(c,:),yQ(c,:),zQ);
    EQ2(c,:) = ES * exp(-ik .* rSQ2(c,:)) ./ rSQ2(c,:);
    
end
    EQ = EQ1 + EQ2;
% irradiance distribution
    we_Q = (cL*nRI*eps0/2) .* abs(EQ).^2;
% radiant flux
    for c = 1 : nQR
        WE_Qring(c) = rQ(c) * drQ * simpson1d(we_Q(c,:),0,2*pi);
    end
        WE_Q = sum(WE_Qring);

        
% ========================================================================
% Observation / Ring Structure 
% ========================================================================    
    
    rP = linspace(eps,rPmax,nPR);  % radii of each ring
    drP = rP(2) - rP(1);

    tP = linspace(0,2*pi,nPA);  % angles for grid points
    xP = zeros(nPR,nPA);  % X values for grid points
    yP = zeros(nPR,nPA);  % Y values for grid points
       
for c = 1 : nPR
    xP(c,:) = rP(c) .* cos(tP);
    yP(c,:) = rP(c) .* sin(tP);
end

% X & Y axes / Optical coordinates - radial
%  nXY:   +X (1)   -X (3)   +Y (2)   -Y (4)

    nXY = [1 (nPA+1)/2 ; 1+(nPA-1)/4 nPA - (nPA-1)/4];
      
%   
A = zeros(1,nQR);
EP = zeros(nPR,nPA);
for c1 = 1 : nPR            % P rings
for c2 = 1 : nPA            % P angles
for c3 = 1 : nQR            % Q rings    
    rPQ = fn_distancePQ(xP(c1,c2),yP(c1,c2),zP,xQ(c3,:),yQ(c3,:),zQ);
    rPQ3 = rPQ .* rPQ .* rPQ;
    kk = ik .* rPQ;
    MP1 = exp(kk);
    MP1 = MP1 ./ rPQ3;
    %unit = ones(1,nQ(c3));
    MP2 = zP .* (ik .* rPQ - 1);
    f = EQ(c3,:) .* MP1 .* MP2;
    A(c3) = rQ(c3) * simpson1d(f,0,2*pi); 

end
    EP(c1,c2) = drQ * sum(A)/ (2*pi);
end
end

we_P = (cL*nRI*eps0/2) .* abs(EP).^2;
we_P_max = max(max(we_P));
we_P_dB = 10 .* log10(we_P./we_P_max);

if abs(xS(1)) == abs(xS(2))
   saddle = we_P(1,nXY(1))/we_P_max;
end

if abs(xS(1)) > abs(xS(2))
    c = 1;
    while we_P(c,nXY(1)) > we_P(c+1,nXY(1));
     c = c+1 ;  
    end
     saddle = we_P(c,nXY(1))/we_P_max;
end

if abs(xS(1)) < abs(xS(2))
    c = 1;
    while we_P(c,nXY(3)) > we_P(c+1,nXY(3));
     c = c+1 ;  
    end
     saddle = we_P(c,nXY(3))/we_P_max;
end

% ========================================================================
% GRAPHICS 
% ========================================================================
graphColor = ColorCode(wL);
% plot 1
hold off
pos = [0.05 0.45 0.4 0.28];
plot1 = subplot('Position',pos);
 tx = 'x_P (blue) &  y_P (red)  [\mum]';
   ty = 'irradiance w_e  [W.m^{-2}]';
   t1 = 'saddle point ratio =  ';
   t2 = num2str(saddle,3);
   tm = [t1 t2];
   x = xP(:,nXY(1)).*1e6;   % +X axis  
  % if xS(1) == xS(2), x =  (xP(:,nXY(1))+xS(1)).*1e6; end;
   y = we_P(:, nXY(1));
   plot(x,y,'b','linewidth',2)
   hold on
   x = xP(:,nXY(3)).*1e6;   % -X axis
  % if xS(1) == xS(2), x =  (xP(:,nXY(1))+xS(1)).*1e6; end;
   y = we_P(:, nXY(3));
   plot(x,y,'b','linewidth',2)

   x = yP(:,nXY(2)).*1e6;   % +Y axis
   y = we_P(:, nXY(2));
   plot(x,y,'r','linewidth',0.5)

   x = yP(:,nXY(4)).*1e6;   % -Y axis
   y = we_P(:, nXY(4));
   plot(x,y,'r','linewidth',0.5)

  xlabel(tx);
  ylabel(ty);
  title(tm)
  grid on
  set(gca,'fontsize',fs);
  
% plot 2
alpha12 = atan(xS./zS);
alpha = sum(abs(alpha12));
alpha_deg = alpha * 180 / pi;
t1 = 'ang. separation  =  ';
t2 = num2str(alpha,'%3.2e');
t3 = ' rad =  ';
t4 = num2str(alpha_deg,2);
t5 = '  deg  ';
tm = [t1 t2 t3 t4 t5];
hold off
pos = [0.55 0.45 0.4 0.28];
plot2 = subplot('Position',pos);
tx = 'x_P (blue) &  y_P (red)  [\mum]';
   ty = 'irradiance w_e  [dB]';
   x = xP(:,nXY(1)).*1e6;   % +X axis                
   y = we_P_dB(:, nXY(1));
   plot(x,y,'b','linewidth',2)
   hold on
   x = xP(:,nXY(3)).*1e6;   % -X axis
   y = we_P_dB(:, nXY(3));
   plot(x,y,'b','linewidth',2)

   x = yP(:,nXY(2)).*1e6;   % +Y axis
   y = we_P_dB(:, nXY(2));
   plot(x,y,'r','linewidth',0.5)

   x = yP(:,nXY(4)).*1e6;   % -Y axis
   y = we_P_dB(:, nXY(4));
   plot(x,y,'r','linewidth',0.5)

  xlabel(tx);
  ylabel(ty);
  title(tm);
  grid on
  set(gca,'fontsize',fs);
  
  % plot 3
  hold off
  pos = [0 0.08 0.25 0.25];
  plot3 = subplot('Position',pos);
   x = xP; y = yP;
   %z = we_P_dB;
   z = we_P.^0.9;
   pcolor(x,y,z)
   shading interp
   colormap(hot)
   %colorbar 
   xlabel('x_P - x_S  [m])');
   ylabel('y_P  [m]')
   axis equal
   axis off
   title('short time exposure');
   
 % plot4
 hold off
   pos = [0.25 0.08 0.25 0.25];
   plot4 = subplot('Position',pos);
   % x = xP; y = yP;
   % z = we_P_dB;
   surf(x,y,z)
   shading interp
   colormap(hot)
   axis off
   view(20,70);
   
% plot 5
hold off
   pos = [0.50 0.08 0.25 0.25];
   plot3 = subplot('Position',pos);
   z = we_P.^0.2;
   pcolor(x,y,z)
   shading interp
   colormap(hot)
   %colorbar 
   xlabel('x_P - x_S  [m])');
   ylabel('y_P  [m]')
   axis equal
   axis off
    title('longer time exposure');
   
 % plot 6
 hold off
   pos = [0.75 0.08 0.25 0.25];
   plot4 = subplot('Position',pos);
   % x = xP; y = yP;
   % z = we_P_dB;
   surf(x,y,z)
   shading interp
   colormap(hot)
   axis off
   view(20,70);