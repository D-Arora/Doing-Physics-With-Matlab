% climateG.m

clear 
close all
clc


% CONSTANTS ============================================================

  p0 = 1.01325e5;     % Sea level standard atmospheric pressure	[Pa]
  L_d = 9.76e-3;      % Temperature lapse rate  [K/m]
  c_pd = 1004.65805;  % Constant-pressure specific heat for dry air [J/(kg·K)]
  T0 = 288.16;	      % Sea level standard temperature [K]
  g0 = 9.80665;       % Earth-surface gravitational acceleration [m/s2]
  M = 0.02896968;     % Molar mass of dry air	[kg/mol]
  R0 = 8.314462618;   % Universal gas constant [J/(mol·K)]
  R_d = 2.8704e2;     % Gas constant dry air

  ME = 5.9722e24;   % mass of Earth  [kg]
  RE = 6.3743e6;    % radius of Eath  [m]
  G = 6.67408e-11;  % Gravitation constant  [m3/(kg.s2)]

% SETUP  ==============================================================
  z = linspace(0,1e4,999);    % altitude   [m]    
  g = (G*ME)./(RE+z).^2;

% Temperature variation with altitude  ================================
  T = T0 - L_d.*z;

% Pressure variation with altitude  ================================
  p = p0.*(1 - L_d.*z/T0).^(g0*M/(R0*L_d));  

% PDry air density variation with altitude  =========================
  rho = p./(R_d.*T);    

% Measurements  =======================================================
zM = [ 0:0.1:0.9 1:0.5:10];

gM = [9.8072 9.8069	9.8066 9.8062 9.8059 9.8056	9.8053 9.805 9.8047	...
    9.8044 9.8041 9.8025 9.801 9.7995 9.7979 9.7964	9.7948 9.7933 ...
    9.7917 9.7902 9.7887 9.7871	9.7856 9.784 9.7825	9.781 9.7794 9.7779	9.7764];

pM = [1013.25 1001.2 989.45	977.72	966.11	954.61	943.22	931.94	920.77 ...
    909.71 898.8 845.59	795	746.9 701.2	657.8 616.6	577.5 540.5	505.4 ...
    472.2 440.7 411.1 383 356.5	331.5 308 285.8	265];
pM = pM.*100;

TM = [288.15 287.5 286.85 286.2	285.55 284.9 284.25	283.6 282.95 282.3 ...
    281.65 278.4 275.15	271.91 268.66 265.41 262.17	258.92 255.68   ...
    252.43 249.19 245.94 242.7 239.46 236.22 232.97	229.73 226.49 223.25];

rhoM = [1.225 1.213 1.202 1.19 1.179 1.167 1.156 1.145 1.134 1.123 ...
    1.112 1.058 1.007 0.957	0.909 0.863	0.819 0.777	0.736 0.697	0.66 ...
    0.624 0.59	0.557 0.526	0.496 0.467	0.44 0.414];

zDk = [0 0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.9	1 1.5 2	2.5	3 3.5 4	4.5	5 ...
    5.5	6 6.5 7 7.5	8 8.5 9	9.5	10 11 12 13 14 15 16 17	18 19 20 21	22 ...
    23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39 ...
    40	41	42	43	44	45	46	47	48	49	50	55	60	65	70	75	80 ...
    85	90	95	100];

zD = 1e3*zDk;

TD = [288.15 287.5 286.85 286.2	285.55	284.9	284.25	283.6	282.95 ...
    282.3	281.65	278.4	275.15	271.91	268.66	265.41	262.17 ...
    258.92	255.68	252.43	249.19	245.94	242.7	239.46	236.22 ...
    232.97	229.73	226.49	223.25	216.78	216.65	216.65	216.65	 ...
    216.65	216.65	216.65	216.65	216.65	216.65	217.58	218.57	 ...
    219.57	220.56	221.55	222.54	223.54	224.53	225.52	226.51 ...
    227.5	228.49	230.97	233.74	236.51	239.28	242.05	244.82	...
    247.58	250.35	253.11	255.88	258.64	261.4	264.16	266.93	...
    269.68	270.65	270.65	270.65	260.77	247.02	233.29	219.59	...
    208.4	198.64	188.89	186.87	188.42	195.08];

   gD = (G*ME)./(RE+zD).^2;

   num = length(zD);
   PD = zeros(num,1);
   PD(1) = p0;
   for c = 1:num-1
      PD(c+1) = PD(c) * exp((-2*gD(c)./R_d)*(zD(c+1)-zD(c)) / (TD(c+1) + TD(c)));
   end

   rhoD = PD./(R_d.*TD');  
% 
% pD = zeros(num,1);
%    pD(1) = p0;
%    pD(2) = pD(1) - pD(1)*(g0/(R_d*TD(1)))*(zD(2) - zD(1));
%    for c = 2:num-1
%       pD(c+1) = pD(c-1) - pD(c)*(g0/(R_d*TD(c)))*(zD(c+1) - zD(c-1));
%    end
   % pD(pD<0) = 0;


% Pressure
  Hp = 7.4e3;
  zE = linspace(0,1e5,999);
  pE = p0.* exp(-zE/Hp);
% density
  Hp = 8.7e3;
  zE = linspace(0,1e5,999);
  rhoE = rhoD(1).* exp(-zE/Hp);

% LAPSE RATE  ==========================================================
  F = find(zD>1e4,1);
  F = 1:F;
  LR = fitlm(zD(F),TD(F)) 
  Lcoeff = LR.Coefficients.Estimate
  Lintercept = Lcoeff(1)
  Lslope = Lcoeff(2)*1e3

% dry air lapse rate
  zd = linspace(0,11,999);
  Td = -1e3.*L_d.*zd + Lintercept;

% GRAPHICS  ============================================================

figure(1)  % 11111111111111111111111111111111111111111111111111111111111
   pos = [0.05 0.05 0.35 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 10;
   box 

subplot(2,2,1)   % g
  xP = gD; yP = zD.*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  yticks(0:20:100)
  xlabel('g  [ m.s^{-2} ]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)

subplot(2,2,2)
  xP = TD; yP = zD.*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  yticks(0:20:100)
  xlabel('T  [ K ]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)

subplot(2,2,3)
  xP = pE./100; yP = zE.*1e-3;
    plot(xP,yP,'r','linewidth',2)
hold on
  xP = PD./100; yP = zD.*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  
  yticks(0:20:100)
  xlim([-1e2 11e2])
  xlabel('P [ hPa ]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)

subplot(2,2,4)
   xP = rhoE; yP = zE.*1e-3;
    plot(xP,yP,'r','linewidth',2)
hold on
  xP = rhoD; yP = zD.*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([-0.1 1.25])
  yticks(0:20:100)
  xticks(0:0.25:1.25)
  xlabel('\rho [ kg.m^{-3}')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)  

figure(2)  % 2222222222222222222222222222222222222222222222222222222222
   pos = [0.45 0.05 0.35 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 10;
   box 

subplot(2,2,1)   % g
  xP = gD(F); yP = zD(F).*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  yticks(0:2:10)
  xlabel('g  [ m.s^{-2} ]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)

subplot(2,2,2)
  xP = TD(F); yP = zD(F).*1e-3;
     plot(xP,yP,'b','linewidth',2)
  hold on   
  xP = Td; yP = zd;
     plot(xP,yP,'r','linewidth',2)
  grid on
  yticks(0:2:10)

  text(200,5,'dry air','Color','r')
  txt = sprintf('L_{dry} = %2.1f K/km \n',-1e3*L_d);
  text(200,2,txt,'Color','r')
  text(245,9.5,'environment','Color','b')
  txt = sprintf('L_{env} = %2.1f K/km \n',Lslope);
  text(245,7,txt,'Color','b')

  xlabel('T  [ K ]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)

subplot(2,2,3)
  xP = PD(F)./100; yP = zD(F).*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  yticks(0:2:10)
  xlim([-1e2 11e2])
  xlabel('P [ hPa ]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS)

subplot(2,2,4)
  xP = rhoD(F); yP = zD(F).*1e-3;
     plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([-0.1 1.25])
  yticks(0:2:10)
  xticks(0:0.25:1.25)
  xlabel('\rho [ kg.m^{-3}')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS) 

  figure(3)
   pos = [0.05 0.3 0.25 0.45];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   FS = 10;
   box 
   xP = TD; yP = zD./1e3;
   plot(xP,yP,'b','linewidth',2)
   hold on

   xP = [180 300]; yP = [11 11];
   plot(xP,yP,'r','linewidth',1)
   text(182,5,'troposphere')
   text(272,9,'tropopause','Color','r')

   xP = [180 300]; yP = [50 50];
   plot(xP,yP,'r','linewidth',1)
   text(182,25,'stratoposphere')
   text(272,48,'stratopopause','Color','r')

   xP = [180 300]; yP = [90 90];
   plot(xP,yP,'r','linewidth',1)
   text(182,58,'mesosphere')
   text(272,88,'mesopause','Color','r')
   text(202,95,'thermosphere')

%   xlim([-0.2 10.2e4])
   grid on
  xlabel('T [ K]')
  ylabel('z  [ km ]')
  set(gca','fontsize',FS) 

