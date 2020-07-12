% tp_rod_11001.m
% Heat Conduction along a cylindrical bar
% S.I. units used unless otherwise stated
% Default value and unit are given in ()
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% ../mphome.htm

tic
close all
clear all
clc

% =======================================================================
% INPUTS 
% =======================================================================

% number of spatial elements (80)
    nx = 150;
% number of time steps (20000)
    nt = 60000;
% length of bar (0.2 m)
    L = 0.20;
% radius of bar (0.01 m)
    r1 = 0.01;
    r = r1 .* ones(1,nx);    
% thermal conductivity (400 W/(m.degC)
    k1 = 400;
    k2 = 34;
    k3 = 133;
    k = k1 .* ones(1,nx);
    %k(50:end) = k2;
    %k(100:end) = k3;
% specific heat capacity (380 J/(kg.degC)
    c1 = 380;  
    c2 = 1;
    c3 = 840;
    c = c1 .* ones(1,nx);
    %c(50:end) = c2;
    %c(100:end) = c3;
% density  (8900 kg/m^3)
    rho1 = 8900;
    rho2 = 1;
    rho3 = 2000;
    rho = rho1 .* ones(1,nx);
    %rho(50:end) = rho2;
    %rho(100:end) = rho3;
    
% initial temperature of bar  (0 degC)
    T = zeros(nt,nx);

% intial energy flux density
    J = zeros(nt,nx+1);

% BOUNDARY CONDITIONS  (100   20  degC)
    T(:,1) = 100;
    T(:,end) = 20;
    %T(:,1:end) = 1200;
    %T(:,end/2:end) = 50;
    
% =======================================================================
% CALCULATIONS 
% =======================================================================    

dx = L / nx;                  % width of elements
dt = 0.5 * dx^2 * min(rho) * min(c) / (2 * max(k));    % time increment
%dt = 1e-5;
x = dx/2 + dx .*(0:nx-1);
xJ = 0:dx:L;
t = 0 : (nt-1);
t = dt .* t;
K1 = (k ./ dx);                             % constant
K2 = dt ./ (rho .* c .* dx);                 % constant
area = pi * r(1)^2;   % assume uniform cross-sectional area

for ct = 1 : nt-1
for cx = 2 : nx
    J(ct+1,cx) = K1(cx) * ( T(ct,cx-1) - T(ct,cx) );
    J(ct+1,1) = J(ct+1,2);
    J(ct+1,end) = J(ct+1,end-1);
end  

for cx = 2 : nx-1
    T(ct+1,cx) = T(ct,cx) + K2(cx) * ( J(ct+1,cx) - J(ct+1,cx+1) ) ;
end
end

% energy flux dQ/dt = JA
    dQ_dt = J(end,end) .* area;

% temperature gradient dT/dx
    dT_dx = (T(end,end) - T(end,1)) / (x(end) - x(1));
    
% time constant tau
    cc = 1; indexT = round(0.5*nx+0.5);
    while T(cc,indexT) < 0.63 * T(end,indexT)
        cc = cc + 1;
    end
    tau = t(cc); 
 
 % max time for simulation
 tMax = t(end);
    
% =======================================================================
% GRAPHICS 
% ======================================================================= 

figure(11);

set(gcf,'units','normalized','position',[0.05,0.05,0.3,0.7]); 
fs = 14;

% energy flux density J vs position -------------------------------------------------
index = round(0.1 * nt) .* [1/(0.1*nt) 1 2 5 10];

% Position of axes
   ap1 = 0.15;
   ap = [ap1 0.65 0.75 0.3]; 
   axes('position',ap);

% Plot data
   xp = xJ; yp = J(index,:);
   plot(xp,yp,'linewidth',2);
   set(gca,'fontsize',fs);
   
% Titles
   xt = 'position  x  [m]';
   yt = 'flux density  J  [w/m^2]';
   xlabel(xt,'fontsize',fs);
   ylabel(yt,'fontsize',fs);
   h_title = title(' Legend: time elapsed in seconds');
   set(h_title,'FontSize',12,'fontWeight','normal');

% Legend
   tt1 = '0';
   tt2 = num2str(t(index(2)),'%4.0f\n');
   tt3 = num2str(t(index(3)),'%4.0f\n');
   tt4 = num2str(t(index(4)),'%4.0f\n');
   tt5 = num2str(t(index(5)),'%4.0f\n');
   h_legend = legend(tt1,tt2,tt3,tt4, tt5);
   set(h_legend,'location','North', 'Orientation','horizontal');

% temperature T vs position -------------------------------------------------
  index = round(0.1 * nt) .* [1/(0.1*nt) 1 2 5 10];
  xp = x; yp = T(index,:);

% Position of axes
   ap = [ap1 0.25 0.75 0.30]; 
   axes('position',ap);

% Plot data
   xp = x; yp = T(index,:);
   plot(xp,yp,'linewidth',2);

% Titles
   xt = 'position  x  [m]';
   yt = 'Temperature T  [degC]';
   xlabel(xt,'fontsize',fs);
   ylabel(yt,'fontsize',fs);
  % h_title = title(' Legend: time elapsed in seconds');
  % set(h_title,'FontSize',12,'fontWeight','normal');

set(gca,'fontsize', fs)

% temperature profile ----------------------------------------------------
   ap = [ap1 0.05 0.75 0.1]; 
   axes('position',ap);
   [xSc, ySc] = meshgrid(x,[-1,+1]);
   zSc = [T(end,:)',T(end,:)'];
   pcolor(xSc, ySc, zSc');
   %shadingMap = hot(64);
   colormap(hot)
   % colormap(shadingMap(:,1)*graphColor);
   shading interp
   axis off
   
% =========================================================================
figure(12)
set(gcf,'units','normalized','position',[0.37,0.2,0.3,0.55]); 
fs = 14;

% energy flux density J vs position -------------------------------------------------
%index = round(0.1 * nt) .* [1/(0.1*nt) 1 2 5 10];
index = round((0.25 * nx) .* [1/(0.25*nx) 1 2 3 4]);
% Position of axes
   ap1 = 0.15; ap2 = 0.60; ap3 = 0.75; ap4 = 0.35;
   ap = [ap1 ap2 ap3 ap4]; 
   axes('position',ap);

% Plot data
   %xp = xJ; yp = J(index,:);
   xp = t; yp = J(:,index);
   plot(xp,yp,'linewidth',2);
   set(gca,'Ylim',[-1e5 10e5]);
   set(gca,'fontsize',fs);
   
% Titles
   xt = 'time  t  [s]';
   yt = 'flux density  J  [w/m^2]';
   xlabel(xt,'fontsize',fs);
   ylabel(yt,'fontsize',fs);
   h_title = title(' Legend: x position along rod in meters');
   set(h_title,'FontSize',12,'fontWeight','normal');

% Legend
   tt1 = '0';
   tt2 = num2str(x(index(2)),'%1.2f\n');
   tt3 = num2str(x(index(3)),'%1.2f\n');
   tt4 = num2str(x(index(4)),'%1.2f\n');
   tt5 = num2str(x(index(5)),'%1.2f\n');
   h_legend = legend(tt1,tt2,tt3,tt4, tt5);
   set(h_legend,'location','North', 'Orientation','horizontal');

% temperature T vs time t -------------------------------------------------
  
% Position of axes
   ap2 = 0.10;
   ap = [ap1 ap2 ap3 ap4]; 
   axes('position',ap);

% Plot data
   xp = t; yp = T(:,index);
   plot(xp,yp,'linewidth',2);

% Titles
   xt = 'time  t  [s]';
   yt = 'Temperature T  [degC]';
   xlabel(xt,'fontsize',fs);
   ylabel(yt,'fontsize',fs);
  % h_title = title(' Legend: time elapsed in seconds');
  % set(h_title,'FontSize',12,'fontWeight','normal');

set(gca,'fontsize', fs)

%% =========================================================================
figure(13)    % parameters
set(gcf,'units','normalized','position',[0.68,0.1,0.28,0.75]); 
fs = 12;
 h_title = title('Parameters');
   set(h_title,'FontSize',12,'fontWeight','normal');

plot([0 100],[60 60]);

  px1 = 10; py1 = 98; dpx = 5; dpy = 5.5; px2 = 50;
% length
   tx1 = 'length  L = ';
   tx2 = num2str(L,4);
   tx3 = '  m';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1,tx); set(h_text,'fontsize',fs);

% radius
   tx1 = 'radius  R = ';
   tx2 = num2str(r(1),4);
   tx3 = '  m';
   tx = [tx1 tx2 tx3];
   h_text = text(px2,py1,tx); set(h_text,'fontsize',fs);
   
% thermal conductivity
   tx = 'Thermal conductivity  k  [W.m^{-1}.K^{-1}]';
   h_text = text(px1,py1-1*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'k_1 = ';
   tx2 = num2str(k1,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+1*dpx,py1-2*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'k_2 = ';
   tx2 = num2str(k2,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+5*dpx,py1-2*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'k_3 = ';
   tx2 = num2str(k3,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+10*dpx,py1-2*dpy,tx); set(h_text,'fontsize',fs); 
   
% specific  heat capacity
   tx = 'Specific heat capacity  c  [J.kg^{-1}.K^{-1}]';
   h_text = text(px1,py1-3*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'c_1 = ';
   tx2 = num2str(c1,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+1*dpx,py1-4*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'c_2 = ';
   tx2 = num2str(c2,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+5*dpx,py1-4*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'c_3 = ';
   tx2 = num2str(c3,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+10*dpx,py1-4*dpy,tx); set(h_text,'fontsize',fs); 
   
% density
   tx = 'Density \rho  [kg.m^{-3}]';
   h_text = text(px1,py1-5*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = '\rho_1 = ';
   tx2 = num2str(rho1,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+1*dpx,py1-6*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = '\rho_2 = ';
   tx2 = num2str(rho2,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+5*dpx,py1-6*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = '\rho_3 = ';
   tx2 = num2str(rho3,4);
   tx3 = '  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+10*dpx,py1-6*dpy,tx); set(h_text,'fontsize',fs); 
   
 % Cross-sectional area  
   tx1 = 'Cross-sectional area  A  =  ';
   tx2 = num2str(area,'%2.2e\n');
   tx3 = '   m^2  ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1-8*dpy,tx); set(h_text,'fontsize',fs);  
   
 % Temperature gradient
   tx1 = 'Temperature gradient  dT/dx  = ';
   tx2 = num2str(dT_dx,'%4.2f\n');
   tx3 = '   ^oC.m^{-1} ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1,py1-9*dpy,tx); set(h_text,'fontsize',fs); 
   
 % energy flux densities
   tx = 'At x = L/2 and t = t_{end}';
   h_text = text(px1,py1-11*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'Energy flux density  J  = ';
   tx2 = num2str(J(end,index(3)),'%2.2e\n');
   tx3 = '   W.m^{-2} ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+dpx,py1-12*dpy,tx); set(h_text,'fontsize',fs); 
   
   tx1 = 'Energy flux dQ/dt = J A  = ';
   tx2 = num2str(area*J(end,index(3)),'%2.2e\n');
   tx3 = '   W ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+dpx,py1-13*dpy,tx); set(h_text,'fontsize',fs); 
   
% time constant   
   tx = 'Time constant \tau  (tau):' ;
   h_text = text(px1,py1-15*dpy,tx); set(h_text,'fontsize',fs);
   tx = '      time for T at x = L/2 to increase to 63% T_{final}';
    h_text = text(px1,py1-16*dpy,tx); set(h_text,'fontsize',fs);
   tx = '      T(\tau, L/2) = 0.63 T(t_{end}, L/2)';
   h_text = text(px1,py1-17*dpy,tx); set(h_text,'fontsize',fs); 
   tx1 = '\tau  =  ';
   tx2 = num2str(tau,'%2.2f\n');
   tx3 = '   s ';
   tx = [tx1 tx2 tx3];
   h_text = text(px1+dpx,py1-18*dpy,tx); set(h_text,'fontsize',fs); 
   
  
   
set(gca,'Xlim',[0 100]);
set(gca,'Ylim',[0 100]);
set(gca,'fontsize',fs);

axis off

%% ========================================================================= 
toc







