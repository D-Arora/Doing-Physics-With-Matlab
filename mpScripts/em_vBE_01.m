% em_vBE_01.m
% 15 11 08
% Motion of a charged partilce in uniform cross B and E fields
% S.I. units used unless otherwise stated
% Default value and unit are given in ()
% Ian Cooper   School of Physics   University of Sydney
% ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

clear all
close all
clc

% =========================================================================
% INPUT PARAMETERS 
% =========================================================================

% number of time steps
   N = 2840;

% Dimensions of volume element and field region
   xMax = 0.4;   xMin = -xMax;
   yMax = xMax; yMax = 0.2;   yMin = -yMax;
   zMax = xMax;   zMin = -zMax;
   
   
   xFMax = 0.35;    xFMin = -xFMax;
   yFMax = xFMax;   yFMin = -yFMax;
   zFMax = xFMax;   zFMin = -zMax;   
   
% Particle  e = 1.602e-19   proton m = 1.67e-27
   q = 1.602e-19;     
   m = 1.67e-27;      

% Fields
   B = 0.8;
   E = 25e5;

% Initial velocites and displacements
   ux = 8e5;   uy = 8e5;   uz = 1e4;
   x0 = -0.4;   y0 = 0;   z0 = 0;

   
% =========================================================================
% Setup
% =========================================================================

% time step
   if B == 0; h = 1e-9;
   else
     h = abs(0.01 * m / (q * B));
   end
   %h = 1e-10;
   
% constants
   k1 = q * B * h / (2 * m);   k2 = q * E *h^2 / m;   k3 = 1/(1+k1^2);
% Initialize arrays
   t = [1:N] .* h;
   x  = zeros(N,1);    y  = zeros(N,1);    z  = uz .* t';
   vx = zeros(N,1);   vy = zeros(N,1);   vz = uz .* ones(N,1);
   v = zeros(N,1);
   ax = zeros(N,1);   ay = zeros(N,1);   az = zeros(N,1);
   
% Time step 1: n = 1   
  x(1)  = x0;    y(1)  = y0;
  vx(1) = ux;   vy(1) = uy;
  ax(1) = q * B * uy / m;   ay(1) = (q/m) * (E- ux * B);

% Time step 2: n = 2
x(2) = x(1) + vx(1)*h;   y(2) = y(1) + vy(1)*h;

% =========================================================================
% TIME LOOPS
% =========================================================================

% Displacement
   for n = 2 : N-1
      if abs(x(n)) < xFMax && abs(y(n)) < yFMax
         x(n+1) = k3 * (2*x(n) + (k1^2-1)*x(n-1) + ...
               2*k1*(y(n) - y(n-1)) + k1*k2);
         y(n+1) = 2*y(n) - y(n-1) - k1*x(n+1) + k1*x(n-1) + k2;  
      else
      x(n+1) = 1 * (2*x(n) + (0-1)*x(n-1) + ...
               2*0*(y(n) - y(n-1)) + 0);
         y(n+1) = 2*y(n) - y(n-1) - 0*x(n+1) + 0*x(n-1) + 0;   
      end
   end
 % Velocity  -------------------------------------------------------------
 for n = 2 : N-1
     vx(n) = (x(n+1) - x(n-1))/(2*h);
     vy(n) = (y(n+1) - y(n-1))/(2*h);
 end
     vx(N) = (x(N)-x(N-1))/h; vy(N) = (y(N)-y(N-1))/h;  
     v = sqrt(vx.^2 + vy.^2 + vz.^2);
     
 % Acceleration
 for n = 1 : N
  if abs(x(n)) < xFMax && abs(y(n)) < yFMax
      ax(n) = (q*B/m) .* vy(n);   ay(n) = (q/m) .* (E - vx(n) .* B);
  else
      ax(n) = 0;   ay(n) = 0;    
  end
 end
  
 
% =========================================================================
% Graphics
% =========================================================================
  
figure(1) % --------------------------------------------------------------
   %set(gcf,'PaperType','A4');
   hold on
   fs = 14;
   set(gcf,'units','normalized','position',[0.05,0.6,0.3,0.3]); 
  
    h_rect = rectangle('Position',[xFMin, yFMin 2*xFMax 2*yFMax]);
    col = [0.8 0.9 0.9];
    set(h_rect,'FaceColor',col,'EdgeColor',col);
    
   xP = x; yP = y;
   plot(xP,yP,'b','LineWidth',2)

   grid on
   box on
   axis equal

   xlabel('x  [m]');
   ylabel('y  [m]');
   axis([xMin xMax yMin yMax]);
   set(gca,'fontsize',fs);
   
  
figure (2) % -------------------------------------------------------------
   set(gcf,'units','normalized','position',[0.36,0.6,0.3,0.3]); 
   xP = t; yP = vx;
   plot(xP,yP,'b','LineWidth',2)
   hold on
   xP = t; yP = vy;
   plot(xP,yP,'r','LineWidth',2)
   yP = v;
   plot(xP,yP,'m','LineWidth',2)
   xlabel('time  t  [s]');
   ylabel('v  [m/s]');
   legend('v_x','v_y', '|v|');
   grid on
   set(gca,'fontsize',14);
   
figure (3) % -------------------------------------------------------------
   set(gcf,'units','normalized','position',[0.670,0.6,0.3,0.3]); 
   plot3(x,y,z,'b','LineWidth',2);
   xlabel('x  [m]'); ylabel('y  [m]'); zlabel('z  [m]');
   grid on
   set(gca,'fontsize',14);

figure (4) % -------------------------------------------------------------
   set(gcf,'units','normalized','position',[0.05,0.1,0.3,0.3]); 
   xP = t; yP = x;
   plot(xP,yP,'b','LineWidth',2)
   hold on
   xP = t; yP = y;
   plot(xP,yP,'r','LineWidth',2)
   xlabel('time  t  (s)');
   ylabel('s  [m]');
   legend('x','y');
   grid on
   set(gca,'fontsize',14);   

figure (5) % -------------------------------------------------------------
   set(gcf,'units','normalized','position',[0.36,0.1,0.3,0.3]); 
   xP = t; yP = ax;
   plot(xP,yP,'b','LineWidth',2)
   hold on
   xP = t; yP = ay;
   plot(xP,yP,'r','LineWidth',2)
   xlabel('time  t  [s]');
   ylabel('a  [m/s^2]');
   legend('a_x','a_y');
   grid on
   set(gca,'fontsize',14);      
  
 figure (6) % -------------------------------------------------------------
   set(gcf,'units','normalized','position',[0.67,0.1,0.3,0.4]); 
   xP = 0; yP = 0;
   plot(xP,yP,'b','LineWidth',2) 
   axis([0 100 0 100]);
   fs = 12;
     px1 = 10; py1 = 98; dpx = 5; dpy = 7; px2 = 50;

% Number of elements  N
   tx1 = 'Number of time steps  N = ';
   tx2 = num2str(N,'%4.0f\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);

% charge q
   py1 = py1 - dpy;
   tx1 = 'Charge  [C]  q = ';
   tx2 = num2str(q,'%2.3e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);

   % charge m
   py1 = py1 - dpy;
   tx1 = 'Mass  [kg]  m = ';
   tx2 = num2str(m,'%2.3e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);
   
% magnetic field B
   py1 = py1 - dpy;
   tx1 = 'Magnetic field [T]  B = ';
   tx2 = num2str(B,'%2.2f\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);
   
% electric field E
   py1 = py1 - dpy;
   tx1 = 'Electric field [V/m]  E = ';
   tx2 = num2str(E,'%2.2e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);
   axis off
   
% initial positions x y z
   py1 = py1 - 1*dpy;
   tx1 = 'Initial values (t = 0 s) for displacement [m]';
   tx2 = ' ';
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);  
   
   py1 = py1 - 1*dpy;
   tx1 = '   x_0 = ';
   tx2 = num2str(x(1),'%2.2f\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);   
   
   py1 = py1 - dpy;
   tx1 = '   y_0 = ';
   tx2 = num2str(y(1),'%2.2f\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs); 
   
   py1 = py1 - dpy;
   tx1 = '   z_0 = ';
   tx2 = num2str(z(1),'%2.2f\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs); 
   
% initial velocities vx vy vz
   py1 = py1 - dpy;
   tx1 = 'Initial values (t = 0 s) for velocity [m/s]';
   tx2 = ' ';
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);  

   py1 = py1 - dpy;
   tx1 = '   u_x = ';
   tx2 = num2str(vx(1),'%2.2e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);   
   
   py1 = py1 - dpy;
   tx1 = '   u_y = ';
   tx2 = num2str(vy(1),'%2.2e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs); 
   
   py1 = py1 - dpy;
   tx1 = '   u_z = ';
   tx2 = num2str(uz(1),'%2.2e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);
   
 % time step h
   py1 = py1 - 2*dpy;
   tx1 = 'Time step [s]  h = ';
   tx2 = num2str(h,'%2.2e\n');
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);
   
   axis off
   