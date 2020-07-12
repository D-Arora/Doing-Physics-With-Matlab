% op_005.m

% NUMERICAL ANAYLSIS OF POLARIZED EM WAVES
%  +Z direction of propagation
%  Electric field E(Ex, Ey) in a XY plane (Ez = 0)

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1003.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181118


clear
close all
clc


tic

% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flag = 1 (save)
   flag1 = 1;  flag2 = 2;
%  Enter file name
     ag_name1 = 'ag_A.gif';
     ag_name2 = 'ag_AA.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.1;  
%  Frame counter start
    nt = 1;


% INPUT SECTION  ======================================================

% Number of time steps [nT = 100] / Number of grid points for Z domain 
nT = 100;
nZ = 100;     
    
% INPUT flagP = 1 or 2
%  flagP = 1: Electric field components
%  flagP = 2: Jones Matrix  A, B and C must be real numbers
%  Electric field amplitude components E0x and E0y  [a.u.]
%  Electric field phases   phix and phiy  [rad]

flagP = 2;

switch flagP
  case 1
    E0x = 2;
    E0y = 1.5;
    phix = 0;
    phiy = 3*pi/4;

    case 2
    A = 2;
    B = 6;
    C = 8;
end
    

% CALCULATION SECTION  ================================================



switch flagP
    case 1
     % Relative phase difference of Ey w.r.t. Ex  [rad]
       phi = phiy-phix;   
     % JONES VECTOR V and Normalized Jones Vector VN
       A = E0x;
       B = E0y*cos(phi);
       C = E0y*sin(phi);
       V = [A; B + 1i*C];
       N = sqrt(A^2 + B^2 + C^2);
       AN = A/N; BN = B/N; CN = C/N;
       VN = V./N;   
       
    case 2
       E0x = A;
       phi = atan2(C,B);
       phix = 0;
       phiy = phi;
       E0y = sqrt(B^2 + C^2);
       V = [A; B + 1i*C];
       N = sqrt(A^2 + B^2 + C^2);
       AN = A/N; BN = B/N; CN = C/N;
       VN = V./N;   
end

% Time domain
  T = 1;                      % period of vibration  [a.u.]
  w = 2*pi/T;                 % angular frequency of vibration  [a.u.]
  t = linspace(0,2*T,nT);       % time  [a.u.]
  uT = exp(-1i*w*t);          % Wave function for time evolution 

  % spatial Z domain
  k = 2*pi;                   % Propagation constant
  z = linspace(0,3,nZ);       % Z domain grid points
  uZ = exp(1i*k*z);           % Spatial Z wavefunction
  
% Electric Field
  ECx = E0x * exp(1i*phix);   % Complex electric field amplitudes
  ECy = E0y * exp(1i*phiy);
  
  
  
  Ex = ECx .* uT;             % time dependent electric field components
  Ey = ECy .* uT;

% Magnitude of electric field vector as a function of time  
   E = sqrt(real(Ex).^2+real(Ey).^2);
% Max magnitude of electric field vector   
   Emax = max(E);
% Min magnitude of electric field vector
   Emin = min(E);
% Index for time when electric field vectro has max magnitude   
   k  = find(E == max(E),1);
% Orientation of major axis of ellipse w.r.t. X axis [deg]   
% alpha = atan2d(real(Ey(k)),real(Ex(k)));
  alpha = atand(real(Ey(k))/real(Ex(k)));
% Eccentricity 
  e = sqrt(Emax^2 - Emin^2);


    
% GRAPHICS SECTION  ===================================================  
  FS = 12;

  
figure(4)
  pos = [0.05 0.05 0.22 0.30];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
for cz = 1 :nT

xP = real(Ex)/Emax;  yP = real(Ey)/Emax;
plot(xP,yP,'k','linewidth',0.5)    
hold on    
xP = [0 real(Ex(cz))./Emax]; yP = [0 0];
plot(xP,yP,'b','linewidth',2)

xP = [0 0]; yP = [0 real(Ey(cz))./Emax];
plot(xP,yP,'r','linewidth',2)

xP = [0 real(Ex(cz))./Emax]; yP = [0 real(Ey(cz))./Emax];
plot(xP,yP,'k','linewidth',3)
xP = real(Ex(cz))./Emax; yP = real(Ey(cz))./Emax;
hPlot = plot(xP,yP,'ko');
set(hPlot,'markerfacecolor','k','markeredgecolor','k','markersize',8)


%xP = real(Ex(cc)); yP = real(Ey(cc));
%arrow([0,0], [xP,yP],'Width',5, 'EdgeColor',[0 0 1],'FaceColor',[0 0 1],'length',25);

hold off
grid on
axis equal
xlim([-1 1])
ylim([-1 1])
set(gca,'xtick',-1:0.5:1)
set(gca,'ytick',-1:0.5:1)
set(gca,'fontsize',12)
xlabel('E_x / E_{max}')
ylabel('E_y / E_{max}')

  tm1 = '\phi = ';
  tm2 = num2str(rad2deg(phi),'%2.0f  deg  \n');
  tm = [tm1 tm2];
  title(tm,'fontsize',14)
  
if flag1 == 1
       frame = getframe(4);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
         imwrite(imind,cm,ag_name1,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name1,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
end

pause(0.000001)

hold off
end  

figure(3)
  pos = [0.05 0.45 0.35 0.28];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  nt = 1;
  for ct = 1 : 5 : nT
      clf
  for cz = 1 :nZ
      hold on

      yP = real(Ex(ct).*uZ)./Emax;
      zP = zeros(nZ,1);
      xP = z;
      plot3(xP,yP,zP,'b','linewidth',2)
      
      zP = real(Ey(ct).*uZ)./Emax;
      yP = zeros(nZ,1);
      xP = z;
      plot3(xP,yP,zP,'r','linewidth',2)

      yP = real(Ex(ct).*uZ)./Emax;
      zP = real(Ey(ct).*uZ)./Emax;
      xP = z;
      plot3(xP,yP,zP,'k','linewidth',2)
      
      yP = real(Ex(ct).*uZ(end))./Emax;
      zP = real(Ey(ct).*uZ(end))./Emax;
      xP = z(end);
      hPlot = plot3(xP,yP,zP,'o');
      set(hPlot,'markerfacecolor','k','markeredgecolor','k','markersize',8)
      
      yP = [0 real(Ex(ct).*uZ(end))./Emax];
      zP = [0 real(Ey(ct).*uZ(end))./Emax];
      xP = [z(end) z(end)];
      plot3(xP,yP,zP,'k');
      
      
      
    
      xlim([0 3]); zlim([-1 1]); ylim([-1 1])     
      zlabel('E_y / E_{max}'); ylabel('E_x / E_{max}'); xlabel('z');
      grid on
      view(81,13)
     % view(143,23)
      set(gca,'fontsize',12)
  end
  
     
     if flag2 == 1
       frame = getframe(3);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
         imwrite(imind,cm,ag_name2,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name2,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
     end
   
   %  pause(0.00001)
  end
  
 %% 
figure(1)
  FS = 12;
  pos = [0.28 0.05 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = t; yP = real(Ex);
    plot(xP,yP,'b','linewidth',3)
  hold on
  yP = real(Ey);
    plot(xP,yP,'r','linewidth',1)
  yP = real(E);
    plot(xP,yP,'k','linewidth',2)

  set(gca,'fontsize',FS)
  xlabel('time  [a.u]')
  ylabel('electric field  [a.u]')
  legend('E_x','E_y','E','orientation','horizontal', 'location','northoutside')
  grid on
  box on

%%
figure(2)
  pos = [0.60 0.05 0.30 0.55];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xlim([0 120])
  ylim([0 120])
  h = 120; dh = -5;

  text(0,h,'POLARIZED LIGHT','fontsize',14)
  h = h+dh;
  text(0,h,'   Elliptical / Circular / Linear','fontsize',14)
%   h = h+dh;
%   text(0,h,'      RIGHT polarized:   0 < \phi < \pi    clockwise rotation','fontsize',14)
%   h = h+dh;
%   text(0,h,'      LEFT polarized:    -\pi < \phi < 0     anticlockwise rotation','fontsize',14)
%   
  dh = -8;
  h = h+1.5*dh;
  text(0,h,'ELECTRIC FIELD','fontsize',14)
  
  tm1 = 'E_{0x}  =  ';
  tm2 = num2str(E0x,'%2.2f  \n');
  tm = [tm1 tm2];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
  tm1 = '\phi_x  =  ';
  tm2 = num2str(phix/pi,'%2.2f  \n');
  tm3 = ' \pi  rad';
  tm = [tm1 tm2 tm3];
  text(50,h,tm,'fontsize',14)
  
  tm1 = 'E_{0y}  =  ';
  tm2 = num2str(E0y,'%2.2f  \n');
  tm = [tm1 tm2];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
  tm1 = '\phi_y  =  ';
  tm2 = num2str(phiy/pi,'%2.2f  \n');
  tm3 = ' \pi  rad';
  tm = [tm1 tm2 tm3];
  text(50,h,tm,'fontsize',14)
  
  tm1 = 'E_{max} = ';
  tm2 = num2str(Emax,'%2.2f  \n');
  tm = [tm1 tm2];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
  
  tm1 = 'E_{min} = ';
  tm2 = num2str(Emin,'%2.2f  \n');
  tm = [tm1 tm2];
  text(50,h,tm,'fontsize',14)
  
  tm1 = '\phi = \phi_y - \phi_x  =  ';
  tm2 = num2str(phi/pi,'%2.2f  \n');
  tm = [tm1 tm2 tm3];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
  
   
  tm1 = 'Orientation of ellipse w.r.t X axis   \alpha =  ';
  tm2 = num2str(alpha,'%2.2f  deg \n');
  tm = [tm1 tm2];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
 
  tm1 = 'Eccentricity of ellipse e =  ';
  tm2 = num2str(e,'%2.2f   \n');
  tm = [tm1 tm2];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
  
  
  h = h+1.5*dh;
  text(0,h,'JONES VECTOR   [A ;   B + i  C]','fontsize',14)
  tm1 = num2str(V,'%2.4f ');
  h = h+dh;
  text(10,h,tm1,'fontsize',14);
  
  h = h+1.5*dh;
  text(0,h,'JONES NORMALIZED VECTOR [A_N ;  B_N + i C_N]','fontsize',14)
  tm1 = num2str(VN,'%2.4f ');
  h = h+dh;
  text(10,h,tm1,'fontsize',14);
  
  tm1 = 'Normalization factor N =  ';
  tm2 = num2str(N,'%2.4f   \n');
  tm = [tm1 tm2];
  h = h+dh;
  text(5,h,tm,'fontsize',14)
axis off


toc
  




