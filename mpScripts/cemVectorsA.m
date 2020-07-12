% cemVectorsA.m
% 10 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Orthogonal components of a [3D] vector
% Inputs: V(1)   V(2)   V(3)   Cartesian component of vector V
% Outputs:   Cylindical components Vrho   Vphi   Vz
% Outputs:   Spherical components  Vtheta   Vphi   VR

clear all
close all
clc

% INPUT vector in Cartesian components ------------------------------------ 
    V = [sqrt(3/2) sqrt(1/2) sqrt(2)];

% CALCULATIONS ------------------------------------------------------------
% angles calculated in radians
% angles from radians to degrees  VthetaD and VphiD
   [Vphi, Vrho, Vz]   = cart2pol(V(1),V(2),V(3));
   [Vtheta, Vphi, VR] = cart2sph(V(1),V(2),V(3));
  
   VphiD = rad2deg(Vphi);
   VthetaD = rad2deg(Vtheta);
  
% GRAPHICS ---------------------------------------------------------------  
figure(1)
   set(gcf,'units','normalized','position',[0.2 0.1 0.5 0.5]);
   fs = 16;
   pos = [0.5 0.3 0.5 0.6];
   subplot('Position',pos);

   % vector V
   xP = [0 V(1)]; yP = [0 V(2)]; zP = [0, V(3)];
   h = plot3(xP,yP,zP);
   set(h,'color',[0 0 1],'linewidth',3);
  
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

   
   xlabel('x'); ylabel('y'); zlabel('z');
   grid on
   axis equal
   box on
   set(gca,'fontsize',14);
 
 %   Text output 
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
   
   tm1 = 'magnitude V_R  =  ';
   tm2 = num2str(VR,'%2.3e\n');
   tm = [tm1 tm2];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'XY magnitude,  V_\rho   =  ';
   tm2 = num2str(Vrho,'%2.3e\n');
   tm = [tm1 tm2];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'azimuthal angle,  V_\phi   =  ';
   tm2 = num2str(Vtheta,'%2.2f\n');
   tm3 = '   rad  =  ';
   tm4 = num2str(VthetaD,'%2.2f\n');
   tm5 = '^o';
   tm = [tm1 tm2 tm3 tm4 tm5];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'polar angle,  V_\theta   =  ';
   tm2 = num2str(Vphi,'%2.2f\n');
   tm3 = '   rad  =  ';
   tm4 = num2str(VphiD,'%2.2f\n');
   tm = [tm1 tm2 tm3 tm4 tm5];
   ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   
   box on
   axis off