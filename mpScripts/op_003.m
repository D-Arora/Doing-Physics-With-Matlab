% op_003.m

% JONES VECTORS: Matrix treatment of polarized light

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1003.htm
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
A1 = 2;
B1 = 6; 
C1 = 8;

% Number of time steps for animation
nT = 100;


% CALCULATION SECTION =================================================
% Normalization parameter
   N = sqrt(A1^2 + B1^2 + C1^2);

% Normalized Jones vector parameters
   A = A1 / N;
   B = B1 / N;
   C = C1 / N;

% Relative phase difference between vibration of Ey w.r.t. Ex
   if B1 == 0 && C1 == 0
      phi = 0;
   else
    phi = atan2(C,B);
   end

% Spatial components of electric field and complex electric field amplitude
   b = sqrt(B^2 + C^2);
   E_0x = A ;
   E_0y = b; 

   E_Cx = A;
   E_Cy = b*exp(1i*phi);

% Time and temporal wavefunction   
   t = linspace(0,2*pi,nT);
   w = 1; 
   u_T = exp(-1i*w*t);

% Electric field components (function of space and time)
   Ex = E_Cx .* u_T;
   Ey = E_Cy .* u_T;

% Orientation angle of ellipse
   K = 2*E_0x*E_0y*cos(phi)/(E_0x^2 - E_0y^2);
   alpha = rad2deg(atan(K))/2;

% Orientation for linear polarized light
   theta = rad2deg(atan(real(E_Cy/E_Cx)));


% GRAPHICS SECTION ====================================================
figure(1)
  FS = 12;
  pos = [0.5 0.05 0.35 0.4];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  subplot('position',[0.05 0.1 0.2 0.8])
     xlim([0 100])
     ylim([0 100])
     h = 100; dh = -9;
     text(0,h,'INPUTS: Jones Vector','fontsize',14)
     tm1 = 'A_1  =  ';
     tm2 = num2str(A1,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'fontsize',14)
     tm1 = 'B_1  =  ';
     tm2 = num2str(B1,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'fontsize',14)
     tm1 = 'C_1  =  ';
     tm2 = num2str(C1,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'fontsize',14)
     h = h+dh;
     text(0,h,'Normalized Jones Vector','fontsize',14)
     h = h+dh;
     text(0,h, '     A      (B + \iti C)','fontsize',14)
     tm1 = 'A  =  ';
     tm2 = num2str(A,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'fontsize',14)
     tm1 = 'B  =  ';
     tm2 = num2str(B,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'fontsize',14)
     tm1 = 'C  =  ';
     tm2 = num2str(C,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'fontsize',14)
     tm1 = 'E_{0x}  =  ';
     tm2 = num2str(E_0x,'%2.2f  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(0,h,tm,'color','r','fontsize',14)
     tm1 = 'E_{0y}  =  ';
     tm2 = num2str(E_0y,'%2.2f  \n');
     tm = [tm1 tm2];
     text(90,h,tm,'color','r','fontsize',14)
     
     tm1 = '\theta  =  ';
     tm2 = num2str(theta,'%2.2f  deg  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'color','r','fontsize',14)
     
     tm1 = '\alpha  =  ';
     tm2 = num2str(alpha,'%2.2f  deg  \n');
     tm = [tm1 tm2];
     h = h+dh;
     text(2,h,tm,'color','r','fontsize',14)
     
     axis off
  
  subplot('position',[0.37 0.2 0.65 0.65])
  
  for cc = 1:nT
     plot([-E_0x E_0x],[0 0],'r','linewidth',2)
     hold on
     plot([0 0],[-E_0y E_0y],'r','linewidth',2)
     plot([0 abs(Ey(1))*cos(phi)],[0 abs(Ey(1))*sin(phi)],'m','linewidth',2)
     plot([0,0],[-1,1],'k')
     plot([-1,1],[0,0],'k')
     plot(real(Ex),real(Ey),'b','linewidth',1)
     
     xP = real(Ex(cc)); yP = real(Ey(cc));
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
     tm1 = '\phi / \pi  =  ';
     tm2 = num2str(phi / pi,'%2.2f  \n');
     tm = [tm1 tm2];
     title(tm,'fontsize',14)
     
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
