% op_004.m

% JONES MATRICES: Matrix treatment of polarized light
%   An optical element is represented by 2x2 matrix
% Variable naming conventions
%   incident light 1 / emergent light 2
%   E electric field / components Ex and Ey   [a.u.]
%      0 amplitude / C complex amplitude
% normalized electric field |E| = 1 a.u.
% phi relative phase of Ey w.r.t. Ex [rad]
% Jones Vector: A B C
% Normalized Jones Vector  AN BN CN
% Jones Matrix [ a b ; c b]
% uT time dependence  (wavefunction)
% Calls external function  arrow.m

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1004.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181111


clear
close all
clc
tic

% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 1;
%  Enter file name
     ag_name = 'ag_A.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.1;  
%  Frame counter start
    nt = 1;

% INPUT SECTION  ======================================================
% Jones Vector parameters
A1 = 1;
B1 = 1; 
C1 = 0;

% JONES MATRIX
a = 1; b = 0; c = 0; d = 1i;


% Rotator
%beta = -45;
%a = cosd(beta); b = -sind(beta); c = sind(beta); d = cosd(beta);


% Number of time steps for animation  [nT = 100]
nT = 100;


% CALCULATION SECTION =================================================
% Time and temporal wavefunction   
   t = linspace(0,2*pi,nT);
   w = 1; 
   uT = exp(-1i*w*t);

% Incident electric field (1): Jones Vector
   [AN1, BN1, CN1, phi1] = ABC(A1, B1, C1);
   [E0x1, E0y1, ECx1, ECy1] = electicField(AN1, BN1 , CN1, phi1);
   Ex1 = ECx1 .* uT;
   Ey1 = ECy1 .* uT;

% JONES MATRIX
   JM = [a, b; c, d];

% Emergent electric field (2) 
   JV1 = [A1; B1 + 1i* C1];
   JV2 = JM * JV1;
 
   A2 = abs(JV2(1));
   B2 = real(JV2(2));
   C2 = imag(JV2(2));

  [AN2, BN2, CN2, phi2] = ABC(A2, B2, C2);

  [E0x2, E0y2, ECx2, ECy2] = electicField(AN2, BN2 , CN2, phi2);

   Ex2 = ECx2 .* uT;
   Ey2 = ECy2 .* uT;


% GRAPHICS SECTION ====================================================
figure(1)
  FS = 12;
  pos = [0.1 0.05 0.42 0.30];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

   subplot('position',[0.36 0.1 0.3 0.8]) 
      xlim([0 100])
      ylim([0 100])
      h = 100; dh = -10;
      
      text(0,h,'Jones Vector: incident light','fontsize',14)
      h = h+dh;
      hT = text(2,h,sprintf('A_1 = %2.2f  B_1 = %2.2f  C_1 = %2.2f',A1,B1,C1));
      set(hT,'fontsize',12)
      h = h+dh;
      
      text(0,h,'Jones Matrix','fontsize',14)
      h = h+dh;
      hT = text(2,h,num2str(JM,2));
      set(hT,'fontsize',12)
      h = h+2*dh;
      text(0,h,'Normalized Jones Vector:','fontsize',14)
      h = h+dh;
      text(0,h,'    emergent light','fontsize',14)
      h = h+dh;
      hT = text(2,h,sprintf('A_2 = %2.2f  B_2 = %2.2f  C_2 = %2.2f',AN2,BN2,CN2));
      set(hT,'fontsize',12)
      h = h+dh;
      
      axis off
  
  for cc = 1 : nT
     subplot('position',[0.05 0.3 0.3 0.5])
       plot([-E0x1 E0x1],[0 0],'r','linewidth',2)
       hold on
       plot([0 0],[-E0y1 E0y1],'r','linewidth',2)
       plot([0 abs(Ey1(1))*cos(phi1)],[0 abs(Ey1(1))*sin(phi1)],'m','linewidth',2)
       plot([0,0],[-1,1],'k')
       plot([-1,1],[0,0],'k')
       plot(real(Ex1),real(Ey1),'b','linewidth',1)
       xP = real(Ex1(cc)); yP = real(Ey1(cc));
       arrow([0,0], [xP,yP],'Width',5, 'EdgeColor',[0 0 1],'FaceColor',[0 0 1],'length',25);
       hold off
       axis equal
       grid on
       xlim([-1 1])
       ylim([-1 1])
       set(gca,'xtick',-1:0.5:1)
       set(gca,'ytick',-1:0.5:1)
       set(gca,'fontsize',14)
       xlabel('E_x');
       ylabel('E_y','rotation',0)
       tm1 = 'Incident:  ';
       tm2 = '\phi / \pi  =  ';
       tm3 = num2str(phi1 / pi,'%2.2f  \n');
       tm = [tm1 tm2 tm3];
       title(tm,'fontsize',14,'fontweight','normal') 
             
     subplot('position',[0.7 0.3 0.25 0.5])
       plot([-E0x2 E0x2],[0 0],'r','linewidth',2)
       hold on
       plot([0 0],[-E0y2 E0y2],'r','linewidth',2)
       plot([0 abs(Ey2(1))*cos(phi2)],[0 abs(Ey2(1))*sin(phi2)],'m','linewidth',2)
       plot([0,0],[-1,1],'k')
       plot([-1,1],[0,0],'k')
       plot(real(Ex2),real(Ey2),'b','linewidth',1)
       xP = real(Ex2(cc)); yP = real(Ey2(cc));
       arrow([0,0], [xP,yP],'Width',5, 'EdgeColor',[0 0 1],'FaceColor',[0 0 1],'length',25);
       hold off
       axis equal
       grid on
       xlim([-1 1])
       ylim([-1 1])
       set(gca,'xtick',-1:0.5:1)
       set(gca,'ytick',-1:0.5:1)
       set(gca,'fontsize',14)
       xlabel('E_x');
       ylabel('E_y','rotation',0)
       tm1 = 'Emergent   ';  
       tm2 = '\phi / \pi  =  ';
       tm3 = num2str(phi2 / pi,'%2.2f  \n');
       tm = [tm1 tm2 tm3];
       title(tm,'fontsize',14,'fontweight','normal') 
 
   if flagS == 1
       frame = getframe(1);
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
         
       pause(0.1)
  end
  
 toc  
 
 % ====================================================================
 % FUNCTIONS
 % ====================================================================
 
  function [A, B, C, phi] = ABC(Ain, Bin, Cin)
     N = sqrt(Ain^2 + Bin^2 + Cin^2);
     A = Ain / N;
     B = Bin / N;
     C = Cin / N;
     
     if B == 0 && C == 0
       phi = 0;
     else
       phi = atan2(C,B);
     end
  end
  
  function [E0X, E0Y, ECX, ECY] = electicField(A, B , C, phi)
      b = sqrt(B^2 + C^2);
      E0X = A ;
      E0Y = b; 
      ECX = A;
      ECY = b*exp(1i*phi);
  end
  
