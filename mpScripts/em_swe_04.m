% em_swe_01.m

% 14 nov 15
% Solving the scalar wave equation using the FDTD method
% Ian Cooper
% School of Physics, University of Sydney
% documentation: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts

close all
clear all
clc
tic

% ========================================================================
% INPUTS
% ========================================================================

% Number of spatial grid points
  Nx = 100;
% Number of time steps
  Nt = 300;
% length of simulation region
   L = 100;
% wave or propagation speed
   v = 10;
% Courant number
   S = 1;
% Pulse type   Gaussian flag = 1  /  Rectangular flag = 2    
   flag = 1;
% Gaussian pulse
   A = 0.5;         % pulse amplitude
   s = 2.5;        % pulse width
   np = 1000;       % pulse starting point given by index np
   
% Setup for saving images (im)
     % flag: 0 animated gif NOT saved or 1 saved
         f_gif = 1;  
     % file name for animated gif   
         ag_name = 'ag_pulse_001.gif'; 
     % the delay in seconds before displaying the next image  
         delay = 0.0; 
         
%  Set Boundary Conditions in main FTTD time loop                

% ========================================================================
% SETUP
% ========================================================================

% Initialise wavefunction and set boundary conditions (comments / uncommnet)  
   u = zeros(Nx,Nt); 

% Spatial grid spacing / time step / spatial grid
   hx = L / Nx;
   ht = S * hx / v;
   S2 = S^2;
   x = linspace(0, L , Nx);
   t = linspace(0,ht*Nt,Nt);

if flag == 1     % Gaussian pulse: 1st time step / 2nd time step
   
   for nx = 2 : Nx-1   % time step nt = 1
      % u(nx,1) = A .* exp(-0.5 .* ((x(nx)- x(np))./(s*hx)).^2);
       u(nx,1) = A .* exp(-0.5 .* ((nx-20)./s).^2) ...
                + A .* exp(-0.5 .* ((nx-80)./s).^2);
   end
   
   for nx = 2 : Nx - 1   % time step 2
   %  u(nx,2) = A .* exp(-0.5 .* ((x(nx)-x(np) - v*ht)./(s*hx)).^2);
      u(nx,2) = A .* exp(-0.5 .* ((nx-21)./s).^2) ...
               + A .* exp(-0.5 .* ((nx-79)./s).^2);
     %u(nx,2) = u(nx,1) - 0.5*S2*(u(nx+1,1) - 2*u(nx,1) + u(nx-1,1));
     %u(nx,2) =2*u(nx,1)-u(nx,1) + S2*(u(nx+1,1) - 2*u(nx,1)+u(nx-1,1)) ;  
   end
end

if flag == 2
% Rectangular pulse
   u(10:30,1) = 0.5;
   u(11:31,2) = 0.5;
end

% ========================================================================
% Solving the scalar wave equation FDTD   
% ========================================================================

for nt = 2 : Nt-1
    for nx = 2 : Nx-1
        u(nx,nt+1) = 2*u(nx,nt)- u(nx,nt-1)+ S2*(u(nx+1,nt) - 2*u(nx,nt) + u(nx-1,nt));     
    end
    
    % Set Boundary Conditions    comment / uncomment
      u(1, nt+1)  = 0;               % Fixed end
     % u(Nx, nt+1) = 0;               % Fixed end
      u(Nx,nt+1) = u(Nx-1,nt+1);      % Free end
    % u(1,nt+1) = u(2,nt+1);          % Free end
    %  u(Nx,nt+1) = u(Nx-1,nt);      % ABC
    %  u(1,nt+1)  = u(2,nt);         % ABC
end

% ========================================================================
% GRAPHICS
% ========================================================================
 
   figure(1)
      set(gcf,'units','normalized','position',[0.5,0.1,0.4,0.7]); 
      grid on
      col = 'b'; LW = 2; fs = 14;
      tm1 = 't  =  ';
      tm3 = '  a.u. ';
      t1  = [11 21 31 41 51];     % enter display time indices 
      t2 =  t(t1);                % display times
      subplot(5,1,1)
         xP = x; yP = u(:,t1(1))';
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1])
         tm2 = num2str(t2(1),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
         grid on; box on;
        % set(gca,'fontsize',fs);
        ylabel('u');
        text(27,0.8,'\rightarrow','fontsize',18);
        text(67,0.8,'\leftarrow','fontsize',18);
      subplot(5,1,2)
         xP = x; yP = u(:,t1(2));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1])
         tm2 = num2str(t2(2),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
         grid on; box on
         ylabel('u');
         text(37,0.8,'\rightarrow','fontsize',18);
        text(57,0.8,'\leftarrow','fontsize',18);
     subplot(5,1,3)
         xP = x; yP = u(:,t1(3));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1]);
         tm2 = num2str(t2(3),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
         grid on; box on;
         ylabel('u');
         text(55,0.8,'\rightarrow','fontsize',18);
        text(40,0.8,'\leftarrow','fontsize',18);
     subplot(5,1,4)
         xP = x; yP = u(:,t1(4));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1]);
         tm2 = num2str(t2(4),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
         grid on; box on;
         ylabel('u');
         text(57,0.8,'\rightarrow','fontsize',18);
         text(37,0.8,'\leftarrow','fontsize',18);
     subplot(5,1,5)
         xP = x; yP = u(:,t1(5));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1]);
         tm2 = num2str(t2(5),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);    
         grid on; box on;
         ylabel('u');
         xlabel('x');
         text(67,0.8,'\rightarrow','fontsize',18);
         text(27,0.8,'\leftarrow','fontsize',18);
         
   figure(2) % -----------------------------------------------------------
      set(gcf,'units','normalized','position',[0.05, 0.5, 0.4, 0.3]); 
      grid on
      col = 'b'; LW = 2; fs = 14;
      tm1 = 't  =  ';
      tm3 = '  a.u. ';
      
      xP = x; yP = u(:,end);
      plot(xP, yP,'color',col,'linewidth',LW);
      hold on
      col = [0.4 0.4 0.4]; LW = 2; 
      xP = x; yP = u(:,1);
      plot(xP, yP,'color',col,'linewidth',LW);
      
      axis([0 L -1.1 1.1])
      tm2 = num2str(t(end),'%2.1f\n');
      tm = [tm1 tm2 tm3];
      title(tm,'color','b');
      xlabel('x  [a.u.]'); ylabel('u  [a.u.]');
      set(gca,'fontsize',fs);
      set(gca,'Xtick',0:10:100);
      box on; grid on
      
 % ========================================================================        
 %   ANIMATION: SAVE FILE AS animated gig
 % ========================================================================
 
 figure(3)
    set(gcf,'units','normalized','position',[0.05,0.1,0.35,0.25]);
    set(gcf,'color','w');
    grid on
    col = 'b'; LW = 2; fs = 14;
    c = 1;
    tm1 = 't  =  ';
    tm3 = '  s ';
  
    xP = x; yP = u(:,1); col = [0 0 1]; LW = 1;
    plot(xP, yP,'color',col,'linewidth',LW);
    axis([0 L -1.1 1.1]);
    set(gca,'Xtick',0:10:100);
    LW = 4; col = [0 0 1];
    
for nt = 1: Nt
      xP = x; yP = u(:,nt);
      plot(xP,yP,'color',col,'linewidth',LW);
      axis([0 L -1.1 1.1])
      grid on
      set(gca,'fontsize',fs);
      xlabel('x  [ m ]'); ylabel('s  [mm.]');
      tm2 = num2str(t(nt),'%2.2f\n');
      tm = [tm1 tm2 tm3];
      title(tm);
      pause(0.001);
      %pause
      drawnow;
    
   if f_gif > 0 
      frame = getframe(3);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      % On the first loop, create the file. In subsequent loops, append.
      if nt==1
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
      else
        imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
      end
   end   

end
   
     

% ========================================================================

toc