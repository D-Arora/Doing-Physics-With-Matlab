% cemB01.m
% 07 May 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% Faraday's Law: mutual inductance M / self inductance L
%    time dependent magnetic fields

% SI units used unless stated otherwise

clear all
close all
clc

tic

% ========================================================================
% INPUTS  
% ========================================================================

% Grid points
   N = 5001;
% Square coil: side length sL / radius of wire a / resisitivity rho
% Position of square loop 
   sL = 1;
   a = 1e-3;
   rho = 1.68e-8;
   
   x1 = 0.1;
  % x1 = 0.05;
   x2 = x1 + sL;
   y1 = -sL/2; y2 = sL/2;
 
% current in long straight conducting wire +I in Z direction
   I10 = 2;

% time 
   NT = 5001;
   minT = 0;
   maxT = 5e-3;
   
% =======================================================================
% SETUP 

% =======================================================================

% Grid for square coil   xG   yG
   x = linspace(x1,x2,N);
   y = linspace(y1,y2,N);
   
   [xG, yG] = meshgrid(x,y);
   
% constants   permeability of free space
   mu0 = 4*pi*1e-7;
   
% time
   t =linspace(minT,maxT, NT);
   dt = t(2)-t(1);
   
% =======================================================================
% CALCULATION: ELECTRIC FIELD 
% =======================================================================

% Resistance of copper coil  R
   R = rho*4*sL / (pi*a^2);
   %R = 2*R;
% Mutual inductance for wire and square coil M:  three ways of calculating
  fn = (1./xG);
  ax = x1; bx = x2; ay = y1; by = y2;
  integral2D = simpson2d(fn,ax,bx,ay,by);   
  M1 = (mu0 / (2*pi)) * integral2D;

  fn = 1./x;
  integral1D = simpson1d(fn,ax,bx);
  M2 = (mu0 / (2*pi)) * (by-ay) * integral1D;

  M3 = (y2 - y1)* mu0/(2*pi) * log(x2/x1);

  M = M1;

% Self inductance   L
   L = 8e-7 * sL * (log(sL/a) - 0.52401);
   %L = 2*L;
% current in long straight wire  I1
  I1 = I10 .* t;
  dI1dt = I10;
   
% steady state current I2 in coil  dI2dt = 0 / time constant tau
   I2sat = (M/R) * dI1dt;
   tau = L/R;
% current in coil I2 as a function of time
   I2 = zeros(NT,1);
   K1 = 1*dt*M/L; K2 = -1*dt*R/L;
  
   for n = 2 : NT
       I2(n) = I2(n-1) + K1 * dI1dt + K2 * I2(n-1);
   end

% emf: wire emf1 / coil emf 2 / net emf
   emf1 =  M .* gradient(I1,dt);
   emf2 = -L .* gradient(I2,dt)';
   emf = emf1 + emf2;

% Electric for coil  Einduced
   EavgJ = rho .* I2(end) ./ (pi*a^2);
   EavgL = emf(end) / (4* sL);
  
% analytical solution
k1 = (M/L) * dI1dt; k2 = R/L;
I2a = (k1/k2).*(1 - exp(-k2.*t));

% ========================================================================
% COMMAND WINDOW OUTPUT 
% ========================================================================

disp('  ');
fprintf('coil saturation current  I2sat  =  %2.3e  A\n',I2sat);
disp('   ');
fprintf('coil saturation emf  emf_sat  =  %2.3e  V\n',emf(end));
disp('   ');
fprintf('time constant tau  =  %2.3e  s\n',tau);
disp('   ');
fprintf('electric field line integral EavgL =  %2.3e  V/m\n',EavgL);
disp('   ');
fprintf('electric field line integral EavgJ =  %2.3e  V/m\n',EavgJ);
disp('   ');
fprintf('resistance of coil R =  %2.3e  ohms\n',R);


% ======================================================================= 
% % GRAPHICS 
% % ======================================================================= 
   fs = 12;

   fig = figure(1);   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
   left_color = [0 0 1]; right_color = [1 0 0];
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
   xP1 = t.*1e3; yP1 = I2.*1e3; xP2 = xP1; yP2 = I1.*1e3;
   yyaxis left
   %[hG,h1,h2] = plotyy(xP1,yP1,xP2,yP2);
   plot(xP1,yP1,'b','lineWidth',2);
   xlabel('time  t  (ms)','fontsize',12);
   ylabel('coil  I_2  ( mA )','fontsize',12); 
   
   yyaxis right
   plot(xP2,yP2,'r','lineWidth',2);
   hold on
   xlabel('time  t  (ms)','fontsize',12);
   ylabel('wire  I_1  ( mA )','fontsize',12); 
  % hG(1).LineWidth = 2;
   %hG(1).Color = [0 0 1];
 %  hG(1).YColor = [0 0 1];
   grid on
   box on
   set(gca,'fontSize',fs);
   
  
figure(2)   % 22222222222222222222222222222222222222222222222222222222222
   set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]); 
   xP = t*1e3; yP = emf .* 1e6;
   plot(xP,yP,'m','lineWidth',2);
   hold on
   yP = emf1 .* 1e6;
   plot(xP,yP,'r','lineWidth',2);
   yP = emf2 .* 1e6;
   plot(xP,yP,'b','lineWidth',2)
   
   xlabel('time  t  (ms)','fontsize',12);
   ylabel('emf ( \muV )','fontsize',12); 
   legend('net emf','emf_1','emf_2');
   
   set(gca,'fontSize',fs);
   grid on
   

disp('   ');   
toc
