% oscC00NB.m
% 190524  Matlab 2018b

% SIMULATION OF THE MOTION OF N COUPLED OSCILLATORS
% DYNAMICS OF CRYSTAL LATTICES
% All parameters are changed within the Script
% [default values and units]

% DOING PHYSICS ONLINE 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/oscC00NA.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney


clear
close all
clc

% SETUP ===============================================================

% Propagation constants units  pi/d: kC >> dkC
   kC = 111;
   dkC = 10.5;

% Time span
   tN = 21;
   tMax = 0.5;
   
   
% spring costants kS  [10]
   kS = 10;
% mass of atoms or particles m  [0.1]
   m = 0.1;
% separation distance d bewtween atoms or particles  [1]   
   d = 1;
% Parameters for X axis
   xMax = 0.2*d; xN = 1501;
   
% CALCULATIONS ========================================================
% Propgation constants  /  Wavelengths  wL1 wL2
  k = kC *pi/d; dk = dkC*pi/d;
  wL1= 2*pi/k;  wL2 = 2*pi/dk;

% X axis / time span
   x = linspace(0,xMax,xN);
   t = linspace(0,tMax,tN);

% Angular frequencies omega
   omega  = sqrt(4*kS/m) * abs(sin(k*d/2));
   domega = sqrt(4*kS/m) * abs(sin(dk*d/2));
   

% v0
  v0 = d*sqrt(kS/m);   
% Phase velocity
  vP = omega/k; 
% Group velocity
  vG = domega/dk; 

% Displacement
  e = zeros(xN,tN);       % superpostion
  e1 = zeros(xN,tN);      % sinusoidal signal
  e2 = zeros(xN,tN);      % envelope - modulation 

for cx = 1:xN
 for ct = 1:tN
     e1(cx,ct) = sin( omega *t(ct) - k *x(cx) );
    % e2(cx,ct) = cos( domega*t(ct) - dk*x(cx) )  ;
  e2(cx,ct) = cos( 2*omega*t(ct) - 2*k*x(cx) ) ;
 end
end
% e = e1 .* e2;
  e = e1 + e2;
if vP < 1e-6
    e1(:,:) = 0;
end

if vG < 1e-6 
    e2(:,:) = 0;
end

e = e1 .* e2;


% GRAPHICS ============================================================ 

% Title
% tm1 = '\lambda = ';
% tm2 = num2str(wL, '%2.2f  ');
% tm3 = '   k = ';
% tm4 = num2str(k1*d/pi, '%2.2f ');
% tm5 = '\pi/d';
% tm6 = '   \omega = ';
% tm7 = num2str(omega1, '%2.2f \n ');
% tm8 = '   v_P = ';
% tm9 = num2str(vP, '%2.2f  ');
% tm10 = '   v_G = ';
% tm11 = num2str(vG, '%2.2f  ');
% 
% tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9 tm10 tm11];



 figure(1)

   pos = [0.52 0.2 0.35 0.5];
   set(gcf,'Units','normalized')
   set(gcf,'Position',pos)
   set(gcf,'color','w')
   
for ct = 1: tN
   
    subplot('position',[0.1 0.4 0.8 0.5])
   % Sinusoidal signal
     xP = x(1:end); yP = e1(1:end,ct);
     plot(xP,yP,'r','linewidth',1.5)
     hold on
     xP = t(ct)*vP; yP = e1(1,1);
     hPlot = plot(xP,yP,'or');
     set(hPlot,'markersize',8,'markerfacecolor','r')
   
   % Envelope - modulation
     xP = x; yP = e2(:,ct);
     plot(xP,yP,'m','linewidth',1.5)
     xP = t(ct)*vG; yP = e2(1,1);
     hPlot = plot(xP,yP,'om');
     set(hPlot,'markersize',8,'markerfacecolor','m')
     
   % Displacement
     xP = x; yP = e(:,ct);
     plot(xP,yP,'b','linewidth',2)
   
     
   hold off
   ylim([-1.2 1.2])
   xlim([0 xMax])
   grid on
   set(gca,'fontsize',12)
   xlabel('x ')
   ylabel('left   e   right')
%    title(tm)
   text(10,1.5,num2str(t(ct),'%2.2f'))
   
   
   subplot('position',[0.1 0.05 0.85 0.2])
   
   xlim([0 100])
   ylim([0 80])
   
   fs = 14;
   txt = ['d = ' num2str(d,'%2.1f') '   m = ' num2str(m,'%2.1f') '   k_S = ' num2str(kS,'%2.1f') ];
   Htext = text(2,70,txt);
   set(Htext,'fontsize',fs)
   
   txt = ['sinusoidal signal  \lambda = ' num2str(wL1,'%2.2f') '   k = ' num2str(kC,'%2.1f') '\pi/d' ...
      '   \omega = ' num2str(omega,'%2.2f') '   v_P = ' num2str(vP,'%2.2f')];
   Htext = text(2,40,txt);
   set(Htext,'fontsize',fs,'color','r')
   
   txt = ['envelope  d\lambda = ' num2str(wL2,'%2.2f') '   dk = ' num2str(dkC,'%2.1f') '\pi/d' ...
       '   d\omega = ' num2str(domega,'%2.2f')   '   v_G = ' num2str(vG,'%2.2f')];
   Htext = text(2,10,txt);
   set(Htext,'fontsize',fs, 'color','m')
   
   set(gca,'fontsize',fs)
   axis off
   
   pause(0.1)
   
end


