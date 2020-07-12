% cemVectors.m
%24 feb 16



%%  CELL 1: Clear all variables & Command Window, close all Figure Windows 
close all;
clear all;
clc;

%% CELL 2: Input the Cartesian components for the vectors A, B & C
A = [1 2 3];
B = [-1 -3 -5];
C = [2 4 -3];

%% CELL 3: Cartesian to spherical / cylindrical
% Magnitude (R) and directions of a vector V
% magtiiude R; polar angle phi; azimuthal angle theta; XY plane radius rho
% Input the Cartesian components or the vector that has already been defined for V;
 % Vx = sqrt(3/2); Vy = 1/sqrt(2); Vz = sqrt(2);
  % V = C;
  V = [sqrt(3/2) sqrt(1/2) sqrt(2)];
  %V = [3 5 -6];
% Calculations 
  %R = norm(V);
   %rho = sqrt(V(1)^2+V(2)^2);
   %phi = atan2(V(2),V(1));
   %if phi < 0, phi = phi + 2*pi; end

  %theta = acos(V(3)/R);

  [phi rho z] = cart2pol(V(1),V(2),V(3));
  [theta phi radius] = cart2sph(V(1),V(2),V(3));
  
  phiD = rad2deg(phi);
  thetaD = rad2deg(theta);
  
  close all
  
figure(1)
   set(gcf,'units','normalized','position',[0.2 0.1 0.5 0.5]);
   fs = 16;
   pos = [0.5 0.3 0.5 0.6];
   subplot('Position',pos);
% vector V
   xP = [0 V(1)]; yP = [0 V(2)]; zP = [0, V(3)];
   h = plot3(xP,yP,zP);
   set(h,'color',[0 0 1],'linewidth',3);
  
 %  set(gca,'Position',pos);
   hold on
% Vx
   xP = [0 V(1)]; yP = [0 0];     % 
   h = plot(xP,yP);
   set(h,'color',[1 0 1],'linewidth',2);
% Vy
   xP = [0 0]; yP = [0 V(2)];
   h = plot(xP,yP);
   set(h,'color',[1 0 1],'linewidth',2);
% Vz
   xP = [0 0]; yP = [0 0]; zP =[0 V(3)];
   h = plot3(xP,yP,zP);
   set(h,'color',[1 0 1],'linewidth',2);

   grid on
   xlabel('x'); ylabel('y'); zlabel('z');
   axis equal
   box on
   set(gca,'fontsize',14);
   
   pos2 = [0.05 0.3 0.35 0.6];
   subplot('Position',pos2);
   axis([0 100 0 100]);
   
   h = text(5,100,'[3D] VECTOR');
   set(h,'fontsize',fs,'color',[0 0 1]);
   
   d = 15; tx = 2; ty = 80;
     
   tm1 = 'V_x  =  ';
   tm2 = num2str(V(1),'%2.3e\n');
   tm = [tm1 tm2];
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'V_y  =  ';
   tm2 = num2str(V(2),'%2.3e\n');
   tm = [tm1 tm2];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'V_z  =  ';
   tm2 = num2str(V(3),'%2.3e\n');
   tm = [tm1 tm2];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'magnitude V,  R   =  ';
   tm2 = num2str(radius,'%2.3e\n');
   tm = [tm1 tm2];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'XY magnitude,  \rho   =  ';
   tm2 = num2str(rho,'%2.3e\n');
   tm = [tm1 tm2];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'azimuthal angle,  \phi   =  ';
   tm2 = num2str(theta,'%2.2f\n');
   tm3 = '   rad  =  ';
   tm4 = num2str(thetaD,'%2.2f\n');
    tm5 = '^o'
   tm = [tm1 tm2 tm3 tm4 tm5];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'polar angle,  \theta   =  ';
   tm2 = num2str(phi,'%2.2f\n');
   tm3 = '   rad  =  ';
   tm4 = num2str(phiD,'%2.2f\n');
  
   tm = [tm1 tm2 tm3 tm4 tm5];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   
   box on
   axis off
   
%% CELL 4: Cylindrical to Cartesian coordinates
  % Inputs - Cylindical coordinates: rho, phi [radians], z
      rho = sqrt(2)
      phi = pi/4
      z = 2
  
  % Calculations Cartesian coordinates vector V(Vx, Vy, Vz)
      Vx = rho * cos(phi)
      Vy = rho * sin(phi)
      Vz = z

%% CELL 5: Spherical to Cartesian coordinates
  % Inputs - spherical coordinates: magnitude R, theta [rad] phi [rad]
      R = sqrt(3)
      theta = pi/4
      phi = pi/4
      
  % Calculations Cartesian coordinates vector V(Vx, Vy, Vz)
      Vx = R * sin(theta) * cos(phi)
      Vy = R * sin(theta) * sin(phi)
      Vz = R * cos(theta)    
      
%% CELL 6: DOT (SCALAR) PRODUCT
  % input for Cartesian components for vectors A and B: dot product equals C
     A = [1 0 1]
     B = [0 1 1]

  % magnitude of vectors A and B, and dot product C
      Amag = norm(A)
      Bmag = norm(B)
      C = dot(A,B)
  
 % angle between vectors [radians] and [degrees]    
     theta_rad = acos(C/(Amag * Bmag))
     theta_deg = acosd(C/(Amag * Bmag)) 


%% CELL 7: CROSS (VECTOR) PRODUCT
  % input for Cartesian components for vectors A and B: dot product equals C
     A = [3 5 7]
     B = [2 4 6]

  % magnitude of vectors A and B, and dot product C
      Amag = norm(A)
      Bmag = norm(B)
      C = cross(A,B)
      Cmag = norm(C)
      
 % angle between vectors [radians] and [degrees]    
     theta_rad = asin(Cmag/(Amag * Bmag))
     theta_deg = asind(Cmag/(Amag * Bmag)) 
     
%% CELL 8: Triple Products
    A = [12 3 4]
    B = [4 5 6]
    C = [-1 -3 -5]

    V1 = dot(A,cross(B,C))
    V2 = dot(B,cross(C,A))
    V3 = dot(C,cross(A,B))

    V4 = dot(A,cross(C,B))
    V5 = dot(B,cross(A,C))
    V6 = dot(C,cross(B,A))

    V7 = cross(A,cross(B,C))
    V8 = B .* dot(A,C) - C .* dot(A,B)
    V9 = cross(cross(A,B),C)

%% CELL 9:    [3D] Rotation - vector transformation of coordinates
   theta = zeros(3,3)
   V = [2 3 0];
   
% rotation angle in XY plane  (deg)
   thetaR = 30
% Rotation matrix   
   theta(1,1) = thetaR;
   theta(1,2) = 90-thetaR;
   theta(1,3) = 90;
   
   theta(2,1) = 90+thetaR;
   theta(2,2) = thetaR;
   theta(2,3) = 90;
   
   theta(3,1) = 90;
   theta(3,2) = 90;
   theta(3,3) = 0;

   R = cosd(theta)

% Transformation of vector
   Vdash = R * V'
   
   