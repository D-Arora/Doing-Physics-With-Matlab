% op_001.m

% Graphs and animation of a plane EM wave
% Input wavelength from 380 nm to 780 nm
% EM wave color coded to match wavelength by calling colorCode.m

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1001.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181103

close all
clear 
clc
tic

% SETUP ===============================================================
% S.I. units for all parameters
% Wavelength [m]    380 nm to 780 nm
   wL = 700e-9;
% Irradiance {W/m^2]
   Savg = 1e-3;
% Max z value
   zMax = 2000e-9;
% Number of elements for sapce domain   
   nZ = 500;
% Number of time steps   
   nT = 144;

% Constants
   c0 = 3.00e8;            % speed of light
   eps0 = 8.85418782e-12;  % permittivy of free space
   mu0  = 1.25663706e-6;   % permeability of free space
   
   
% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 0;
%  Enter file name
     ag_name = 'agOP1.gif';
%  Delay in seconds before displaying the next image  
    delay = 0;  
%  Frame counter start
    nt = 1;
   
   
% CALCULATIONS ========================================================
  f = c0/wL;
  T = 1/f;
  w = 2*pi*f;
  k = 2*pi/wL;
  E0 = sqrt(2*Savg/(c0*eps0));
  t = linspace(0,2*T,nT);
  z = linspace(0,zMax,nZ);
  thisColor = ColorCode(wL);
  
  UT = exp(-1i.*w.*t);
  US = E0.*exp(1i.*k.*z);

    
% GRAPHICS ===========================================================

figure(3)
  FS = 12;
  pos = [0.6 0.05 0.35 0.40];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = z.*1e9; 
  for cc = 1:nT  
    yP = real(US.*UT(cc));
    plot(xP,yP,'color',thisColor,'linewidth',2);
     
    xlabel('position  z  [ nm ]')
    ylabel('electric field  E  [ V.m^{-1} ]')
  
    tm1 = '\lambda = ';
    tm2 = num2str(wL*1e9,'%2.0f \n');
    tm3 = '  nm    T = ';
    tm4 = num2str(T,'%2.2e \n');
    tm5 = '  s     f = ';
    tm6 = num2str(f,'%2.2e \n');
    tm7 = '  Hz';
    tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
    title(tm)
    
    tm1 = 't / T =  ';
    tm2 = num2str(t(cc)/T,'%2.1f \n');
    tm = [tm1 tm2];
    text(850,0.92,tm,'fontsize',FS+2)
    
    box on
    grid on
    set(gca,'fontsize',FS) 
  
   if flagS == 1
       frame = getframe(3);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
   end
   pause(0.01)
  end  
  
 
figure(1)
  FS = 12;
  pos = [0.05 0.05 0.35 0.30];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
  tStep = [1, 19, 37];
    xP = z.*1e9;
    yP = real(US.*UT(tStep(1)));
    plot(xP,yP,'color',thisColor,'linewidth',2);
  hold on
    yP = real(US.*UT((tStep(2))));
    plot(xP,yP,'color',thisColor,'linewidth',1);
    
    yP = real(US.*UT((tStep(3))));
    plot(xP,yP,'color',thisColor,'linewidth',1);
    
  xlabel('position  z  [ nm ]')
  ylabel('electric field  E  [ V.m^{-1} ]')
  
  tm1 = 't_1 / T = 0';
  tm2 = '   t_2 / T = ';
  tm3 = num2str(t(tStep(2))/T,'%2.2f \n');
  tm4 = '    t_3 / T = ';
  tm5 = num2str(t(tStep(3))/T,'%2.2f \n');
  tm = [tm1 tm2 tm3 tm4 tm5];
  title(tm)
  set(gca,'xtick',0:250:2000);
  box on
  grid on
  set(gca,'fontsize',FS)
 
 
 figure(2)
  FS = 12;
  pos = [0.35 0.5 0.35 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
    zP = real(US.*UT(1))./max(real(US.*UT(1))); yP = zeros(nZ,1); xP = z.*1e9;
    plot3(xP,yP,zP,'color','b','linewidth',2);
    hold on
   yP = real(US.*UT(1))./max(real(US.*UT(1))); zP = zeros(nZ,1); 
    plot3(xP,yP,zP,'color','r','linewidth',2);
    Harea = area(xP,yP);
    set(Harea,'EdgeColor','r','FaceColor',[1 0.8 0.8]);
  view(33,33)
  xlabel('z  [ nm ]')
  ylabel('H [a.u.]','color','r')
  zlabel('E [a.u.]','color','b','rotation',0)
  
  box on
  grid on
  set(gca,'fontsize',FS) 
  
toc

 
  