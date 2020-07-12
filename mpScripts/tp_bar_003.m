% tp_bar_001.m
% Heat Conduction along a cylindrical bar
% S.I. units used unless otherwise stated
% Default value and unit are given in ()
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

tic
close all
clear all
clc

% =======================================================================
% INPUTS 
% =======================================================================

% number of spatial elements (100)
    nx = 100;
% number of time steps (100)
    nt = 40000;
% length of bar (1 m)
    L = 0.3;
% radius of bar (0.1 m)
    r1 = 0.01;
    r = r1 .* ones(1,nx);    
% thermal conductivity (400 W/(m.degC)
    k1 = 50.2;
    k2 = 385;
    k3 = 1;
    k = k1 .* ones(1,nx);
    k(33:end) = k2 ;
% specific heat capacity (380 J/(kg.degC)
    c1 = 470;  
    c2 = 390;
    c3 = 1;
    c = c1 .* ones(1,nx);
    c(33:end) = c2 ;
% density  (8900 kg/m^3)
    rho1 = 7900;
    rho2 = 8900;
    rho3 = 1;
    rho = rho1 .* ones(1,nx);
    rho(33:end) = rho2 ;
% initial temperature of bar  (0 degC)
    T = zeros(nt,nx+1);

% intial energy flux density
    J = zeros(nt,nx);

% BOUNDARY CONDITIONS  (100   0  degC)
    T(:,1) = 100;
    T(:,end) = 0;
 
    
% =======================================================================
% CALCULATIONS 
% =======================================================================    
 

dx = L / nx;                  % width of elements
dt = 0.4 * dx^2 * min(rho) * min(c) / (2 * max(k));    % time increment
%dt = 1e-5;
x = dx .*(0:nx);
xJ = dx/2 + dx .*(0:nx-1);
t = 0 : (nt-1);
t = dt .* t;
K1 = -( k ./ dx );                             % constant
K2 = dt ./ ( rho .* c .* dx );                 % constant
area = pi * r(1)^2;   % assume uniform cross-sectional area

for ct = 1 : nt-1
       
    for cx = 1 : nx
        J(ct,cx) = K1(cx) * ( T(ct,cx+1) - T(ct,cx) );
    end 

    for cx = 2 : nx-1
        T(ct+1,cx) = T(ct,cx) +  0.5*K2(cx) * ( J(ct,cx-1) - J(ct,cx+1)) ;
    end
end
   J(end,:) = J(end-1,:);
   
% energy flux dQ/dt = JA
    dQ_dt = J(end,end) .* area;

% temperature gradient dT/dx
    dT_dx = (T(end,end) - T(end,1)) / (x(end) - x(1));
    
% time constant tau
    cc = 1;
    while J(cc,end) < 0.63 * J(end,end)
        cc = cc + 1;
    end
    tau = t(cc); 
    
    
% =======================================================================
% GRAPHICS 
% ======================================================================= 

figure(1)
fs = 12;
set(gcf,'units','normalized','position',[0.05,0.2,0.4,0.6]); 

% temperature vs position
index = round(0.1 * nt) .* [1/(0.1*nt) 1 2 5 10];
tt1 = 't  =  ';
tt2 = num2str(t(index),2);
tt3 = '   s';
tt = [tt1 tt2 tt3];
xt = 'position  x  [m]';
yt = 'temperature  T  [degC]';
xp = x; yp = T(index,:);

subplot(2,1,1)
plot(xp,yp,'linewidth',2);
xlabel(xt,'fontsize',fs);
ylabel(yt,'fontsize',fs);
title(tt);
set(gca,'fontsize',fs);

% energy flux density vs position
tt1 = 'plots at times  [s]   ';
tt2 = num2str(t(index),2);
tt = [tt1 tt2];

yt = 'flux density  J  [w/m^2]';
xp = xJ; yp = J(index,:);
subplot(2,1,2)
plot(xp,yp,'linewidth',2);
xlabel(xt,'fontsize',fs);
ylabel(yt,'fontsize',fs);
set(gca,'fontsize',fs);
% 
% figure(2)
% fs = 12;
% set(gcf,'units','normalized','position',[0.5,0.2,0.4,0.6]); 
% % temperture vs time
% index = round((0.25 * nx) .* [1/(0.25*nx) 1 2 3 4]);
% tt1 = 'x   =    ';
% tt2 = num2str(x(index),2);
% tt3 = '   m';
% tt = [tt1 tt2 tt3];
% xt = 'time  t  [s]';
% yt = 'temperature  T  [degC]';
% xp = t; yp = T(:,index);
% subplot(2,1,1)
% plot(xp,yp,'linewidth',2);
% xlabel(xt,'fontsize',fs);
% ylabel(yt,'fontsize',fs);
% title(tt);
% legend(num2str(t(index(1)),2),num2str(t(index(2)),2),num2str(t(index(3)),2) ...
%     ,num2str(t(index(4)),2),num2str(t(index(5)),2) ...
%     , 'orientation','horizontal');
% set(gca,'fontsize',fs);
% 
% % energy flux density vs time
% yt = 'energy flux density  J  [W/m^2]';
% xp = t; yp = J(:,index);
% subplot(2,1,2)
% plot(xp,yp,'linewidth',2);
% xlabel(xt,'fontsize',fs);
% ylabel(yt,'fontsize',fs);
% legend(num2str(x(index(1)),2),num2str(x(index(2)),2),num2str(x(index(3)),2) ...
%     ,num2str(x(index(4)),2),num2str(x(index(5)),2) ...
%     , 'orientation','horizontal');
% f = 1;
% if T(end,1)< T(end,end); f = -1; end;
% set(gca,'Ylim',[0 f*2*max(J(end,end))]);
% set(gca,'fontsize',fs);
% 
% figure(3)
% set(gcf,'units','normalized','position',[0.1,0.05,0.8,0.1]); 
% 
% [xSc, ySc] = meshgrid(x,[-1,+1]);
% zSc = [T(end,:)',T(end,:)'];
% pcolor(xSc, ySc, zSc');
% %shadingMap = hot(64);
% colormap(hot)
% % 
% % colormap(shadingMap(:,1)*graphColor);
% shading interp
% axis off

figure(5)
set(gcf,'units','normalized','position',[0.025,0.025,0.9,0.9]); 
fs = 12;

% temperature vs position
index = round(0.1 * nt) .* [1/(0.1*nt) 1 2 5 10];
tt1 = 't  =  ';
tt2 = num2str(t(index),2);
tt3 = '   s';
tt = [tt1 tt2 tt3];
xt = 'position  x  [m]';
yt = 'temperature  T  [degC]';
xp = x; yp = T(index,:);

axes('position',[0.3,0.65,0.3,0.25]);
plot(xp,yp,'linewidth',2);
%xlabel(xt,'fontsize',fs);
ylabel(yt,'fontsize',fs);
title(tt);
set(gca,'fontsize',fs);

axes('position',[0.3,0.3,0.3,0.25]);
% energy flux density vs position
tt1 = 'plots at times  [s]   ';
tt2 = num2str(t(index),2);
tt = [tt1 tt2];

yt = 'flux density  J  [w/m^2]';
xp = xJ; yp = J(index,:);
plot(xp,yp,'linewidth',2);
xlabel(xt,'fontsize',fs);
ylabel(yt,'fontsize',fs);
f = 1;
if T(end,1)< T(end,end); f = -1; end;
set(gca,'Ylim',[0 f*2*max(J(end,:))]);
set(gca,'fontsize',fs);

axes('position',[0.65,0.65,0.3,0.25]);
index = round((0.25 * nx) .* [1/(0.25*nx) 1 2 3 4]);
tt1 = 'x   =    ';
tt2 = num2str(x(index),2);
tt3 = '   m';
tt = [tt1 tt2 tt3];
xt = 'time  t  [s]';
yt = 'temperature  T  [degC]';
xp = t; yp = T(:,index);
plot(xp,yp,'linewidth',2);
%xlabel(xt,'fontsize',fs);
%label(yt,'fontsize',fs);
title(tt);
%legend(num2str(t(index(1)),2),num2str(t(index(2)),2),num2str(t(index(3)),2) ...
%    ,num2str(t(index(4)),2),num2str(t(index(5)),2) ...
%    , 'orientation','horizontal');
set(gca,'fontsize',fs);

axes('position',[0.65,0.3,0.3,0.25]);
% energy flux density vs time
yt = 'energy flux density  J  [W/m^2]';
xp = t; yp = J(:,index);
plot(xp,yp,'linewidth',2);
xlabel(xt,'fontsize',fs);
%ylabel(yt,'fontsize',fs);
%legend(num2str(x(index(1)),2),num2str(x(index(2)),2),num2str(x(index(3)),2) ...
%    ,num2str(x(index(4)),2),num2str(x(index(5)),2) ...
%    , 'orientation','horizontal');
f = 1;
if T(end,1)< T(end,end); f = -1; end;
set(gca,'Ylim',[0 f*2*max(J(end,end))]);
set(gca,'fontsize',fs);

axes('position',[0.3,0.1,0.305,0.1]);
[xSc, ySc] = meshgrid(x,[-1,+1]);
zSc = [T(end,:)',T(end,:)'];
pcolor(xSc, ySc, zSc');
%shadingMap = hot(64);
colormap(hot)
% 
% colormap(shadingMap(:,1)*graphColor);
shading interp
axis off


axes('position',[0.02,0.1,0.2,0.8]);
title('PARAMETERS')
%axis([0,0,100,100]);
plot(50,50);
set(gca,'Xlim',[0,100]);
set(gca,'Ylim',[0,100]);
title('PARAMETERS');

tx1 = 'L = ';
tx2 = num2str(L,4);
tx3 = '  m';
tx = [tx1 tx2 tx3];
text(10,98,tx);

tx1 = 'radius = ';
tx2 = num2str(r(1),4);
tx3 = '  m';
tx = [tx1 tx2 tx3];
text(10,94,tx);

tx1 = 'k_1 = ';
tx2 = num2str(k1,4);
tx3 = '  W.m^{-1}.^oC^{-1}';
tx = [tx1 tx2 tx3];
text(10,90,tx);

tx1 = 'k_2 = ';
tx2 = num2str(k2,4);
tx3 = '  W.m^{-1}.^oC^{-1}';
tx = [tx1 tx2 tx3];
text(10,86,tx);

tx1 = 'k_3 = ';
tx2 = num2str(k3,4);
tx3 = '  W.m^{-1}.^oC^{-1}';
tx = [tx1 tx2 tx3];
text(10,82,tx);

tx1 = 'c_1 = ';
tx2 = num2str(c1,4);
tx3 = '  J.kg^{-1}.^oC^{-1}';
tx = [tx1 tx2 tx3];
text(10,78,tx);

tx1 = 'c_2 = ';
tx2 = num2str(c2,4);
tx3 = '  J.kg^{-1}.^oC^{-1}';
tx = [tx1 tx2 tx3];
text(10,74,tx);

tx1 = 'c_3 = ';
tx2 = num2str(c3,4);
tx3 = '  J.kg^{-1}.^oC^{-1}';
tx = [tx1 tx2 tx3];
text(10,70,tx);

tx1 = '\rho_1 = ';
tx2 = num2str(rho1,4);
tx3 = '  kg.m^{-3}';
tx = [tx1 tx2 tx3];
text(10,66,tx);

tx1 = '\rho_2 = ';
tx2 = num2str(rho2,4);
tx3 = '  kg.m^{-3}';
tx = [tx1 tx2 tx3];
text(10,62,tx);

tx1 = '\rho_3 = ';
tx2 = num2str(rho3,4);
tx3 = '  kg.m^{-3}';
tx = [tx1 tx2 tx3];
text(10,58,tx);

tx1 = 'cross-sectional area = ';
tx2 = num2str(area,4);
tx3 = '  m^{2}';
tx = [tx1 tx2 tx3];
text(10,54,tx);

tx = 'At x = L and time t = t_{end}';
text(10,48,tx);
tx1 = '   Energy flux dQ/dt = J A =  ';
tx2 = num2str(dQ_dt,4);
tx3 = '  W';
tx = [tx1 tx2 tx3];
text(10,40,tx);
tx1 = '   Energy flux density  J =  ';
tx2 = num2str(J(end,end),4);
tx3 = '  W';
tx = [tx1 tx2 tx3];
text(10,44,tx);

tx = 'At x = L and t = tau ';
text(10,34,tx);
tx = '    Time for J(tau,L) = 0.63 J(t_{end},L)';
text(10,30,tx);
tx1 = '    time constant  \tau  =  ';
tx2 = num2str(tau,4);
tx3 = '  s';
tx = [tx1 tx2 tx3];
text(10,26,tx);

tx1 = 'Temperature gradient  dT/dx  = ';
tx2 = num2str(dT_dx,4);
tx3 = '  ^oC.m^{-1}';
tx = [tx1 tx2 tx3];
text(10,22,tx);

axis off


toc







