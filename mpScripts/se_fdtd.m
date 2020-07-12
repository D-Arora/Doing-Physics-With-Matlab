% se_fdtd_01.m

% 12 dec 15
% Solving the Time Dependent Schrodinger Equation using the FDTD Method
% Motion of a wave packet
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
Nx = 1001;         % [1001]   must be an odd number - number of grid points
Nt = 10000;       % [10000]  number of time steps
%Nt = round(2*Nt);

me = 9.10938291e-31;    % electron mass
hbar = 1.054571726e-34; % hbar Planck's constant
e = 1.602176565e-19;    % elementary charge

U0 = -600;   U1 = -610; U2 = 10;   % max value for PE / limits for PE vs x plot
UyTick = [-600 0];
%UyTick = [0 100];

NxC = round(Nx/2);        % centre postion for Simpson's rule Prob1 and Prob2     

L = 4e-9;                 % [ 4e-9]   length

% Show figure 2 for animation  (0 no)  (1 yes) for flag1
   flag1 = 1;    
   fs = 14;            % fontsize
   time_step = 500;    % jump frames
   
% Save animation as a gif file (0 no)  (1 yes) for flag2
    flag2 = 0;    
% file name for animated gif   
    ag_name = 'ag_se_fdtd_20.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.25; 
% Frame to start
    frame1 = 0;

% Initialise the pulse
   nx0 = round(Nx/4);     % pulse centre   round(Nx/4)
   %wL = L/40;             % pulse wavelength   L/40
   s = L/25;              % pulse width        L/25
   wL = 1.6e-10;

   % Potential energy function
% (1) zero   (2) step   (3) electric field
% (4) hill U0 > 0  or well U) < 0   (5)  parabolic well
    flagU = 5;
    
% Graphics:   limits and tick spacing
  yL1 = -1e5; yL2 = 1e5;  yL = yL2/2;      % [1e5]  wavefunction
  yLL1 = 0;   yLL2 = 8e9; yLL = 2e9;       % [6e9]  prob density
    
  
% ========================================================================
% SETUP
% ========================================================================
   h = hbar*2*pi;
   KE = (h/wL)^2/(2*me*e);  % theoretical KE
   x = linspace(0,L,Nx);
   dx = x(2) - x(1);
   yR = zeros(1,Nx);  yI = zeros(1,Nx);
   t = zeros(Nt,1);
   C1 = 1/10;
   dt = C1 * 2 * me * dx^2 / hbar;
   C2 = e*dt/hbar;
   C3 = -hbar^2 / (2 * me * dx^2 * e);
   T = Nt * dt;
   %wL = 40*dx;
   
   
   nStep = round((1:9)*Nt/9);
   
   
% Potential energy function U  -------------------------------------------
   U = zeros(1,Nx);
    
   switch flagU
      case 2                   % potential step
         U(round(NxC):end) = U0; 
      
      case 1
         U = zeros(1,Nx);      % zero PE - free wave packet 
   
     case 3                    % accelerating / retarding potential
          U = -(U0/L).* x + U0;
          a = e*U0 / (me * L);
          F = me * a;
          v1 = 7.1287e+06;
          v2 = v1 + a * T;
          x1 = 9.9620e-10;
          x2 = x1 + v1 * T-dt + 0.5 * a * (T-dt)^2;
          W = F * (x2-x1)/e;
          
     case 4                    % potential hill / well
         w1 = -80; 
         U(NxC-2*40/2:NxC+2*40/2-w1) = U0;
         w = NxC+2*40/2-w1 - (NxC-2*40/2); 
         
     case 5                   % parabolic well
         a = 4*abs(U0)/L^2; b = -4*abs(U0)/L; c = 0;
         U = a .* x.^2 + b .* x + c;
         period = 2*pi*sqrt(me / (2*a*e));
         Nt_period = period/dt  ;
   end      
         

% ========================================================================
% INITIAL WAVE PACKET
% ========================================================================

for nx = 1 : Nx
  
   yR(nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2)*cos(2*pi*(x(nx)-x(nx0))/wL);
   yI(nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2)*sin(2*pi*(x(nx)-x(nx0))/wL);
   %yR(nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2);
   %yR(nx) = sin(2*pi*(x(nx))/wL);
   % yR(nx) = 1*sin(2*pi*(x(nx))/(2*L/5)) + 3*sin(2*pi*(x(nx))/(2*L/4));
end
  % yI(1:Nx) = 0;
  
% Normalize initial wave packet
    y2 = yR.^2 + yI.^2;
    A = simpson1d(y2,0,L);
    yR = yR ./ sqrt(A); yI = yI ./ sqrt(A);
    prob_density = yR.^2 + yI.^2;
    yR1 = yR; yI1 = yI; yP11 = prob_density;

% Kinetic energy KE
         fn = zeros(1,Nx-2);
    for nx = 2 : Nx-1
         fn(nx) = C3 * (yR(nx) - 1i * yI(nx)) * ...
         (yR(nx+1) - 2 * yR(nx) + yR(nx-1)+ 1i *(yI(nx+1) - 2 * yI(nx) + yI(nx-1)));
    end
    K1avg = simpson1d(fn(2:end),0,L);
    
% ========================================================================
% FDTD solution of Schrodinger Equation
% ========================================================================
cs = 1; probD = zeros(Nx,9);

for nt = 1 : Nt
   for nx = 2 : Nx - 1
      yR(nx) = yR(nx) - C1*(yI(nx+1)-2*yI(nx)+yI(nx-1)) + C2*U(nx)*yI(nx);
   end
   
   for nx = 2 : Nx-1
      yI(nx) = yI(nx) + C1*(yR(nx+1)-2*yR(nx)+yR(nx-1)) - C2*U(nx)*yR(nx);
   end

 if nt == nStep(cs)
    probD(:,cs) = yR.^2 + yI.^2;
    cs = cs +1 ;
 end
   
  if flag1 == 1
      
   if mod(nt,time_step) == 0
   figure(1) 
   fs = 14;
   set(gcf,'units','normalized','position',[0.1 0.1 0.42 0.8]);
   subplot(3,1,1) % ------------------------------------------------------
   [hAx,hLine1,hLine2] = plotyy(x,yR,x,U);
   
   set(hAx(1),'yLim',[yL1 yL2]);
   set(hAx(1),'yTick',yL1 : yL : yL2);
   set(hAx(1),'fontsize',fs,'YColor',[0 0 1]);
   set(hLine1,'linewidth',2,'Color',[0 0 1]);
      
   set(hAx(2),'yLim',[U1 U2]);
   if U0 == 0 
      set(hAx(2),'yTick',0);
   else
      set(hAx(2),'yTick',UyTick);
   end 
   set(hAx(2),'fontsize',fs,'YColor',[0 0 0]);
   set(hLine2,'linewidth',1,'Color',[0 0 0]);
   
   % Potential energy  U
    fn = U .* (yR.^2 + yI.^2);
    Uavg = simpson1d(fn,0,L);

% Kinetic energy KE
    for nx = 2 : Nx-1
         fn(nx) = C3 * (yR(nx) - 1i * yI(nx)) * ...
         (yR(nx+1) - 2 * yR(nx) + yR(nx-1)+ 1i *(yI(nx+1) - 2 * yI(nx) + yI(nx-1)));
    end
    Kavg = simpson1d(fn,0,L);
 
% Total energy   
    Eavg = Uavg + Kavg;
   
   xlabel('x  [m]','fontsize',fs);
   ylabel(hAx(1),'\psi_{real}','fontsize',2*fs,'Color',[0 0 1]);
   ylabel(hAx(2),'U  [eV]','fontsize',1*fs,'Color',[0 0 0]);
  
   tm1 = 't_{steps} = ';
   tm2 = num2str(nt,'%5.0f\n');
   tm3 = '  <E> = ';
   tm4 = num2str(real(Eavg),'%3.0f\n');
   tm5 = '  <K> = ';
   tm6 = num2str(real(Kavg),'%3.0f\n');
   tm7 = '  <U> = ';
   tm8 = num2str(real(Uavg),'%3.0f\n');
   tm9 = '   eV';
   
   tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9];
   
   h_title = title(tm,'fontsize',fs);
   set(gca,'fontsize',fs);
   box on; grid on;  
   
 subplot(3,1,2) % ------------------------------------------------------
   [hAx,hLine1,hLine2] = plotyy(x,yI,x,U);
   
   set(hAx(1),'yLim',[yL1 yL2]);
   set(hAx(1),'yTick',yL1 : yL : yL2);
   set(hAx(1),'fontsize',fs,'YColor',[0 0 1]);
   set(hLine1,'linewidth',2,'Color',[0 0 1]);
      
   set(hAx(2),'yLim',[U1 U2]);
   if U0 == 0 
      set(hAx(2),'yTick',0);
   else
      set(hAx(2),'yTick',UyTick);
   end 
   set(hAx(2),'fontsize',fs,'YColor',[0 0 0]);
   set(hLine2,'linewidth',1,'Color',[0 0 0]);
   
   xlabel('x  [m]','fontsize',fs);
   ylabel(hAx(1),'\psi_{imag}','fontsize',2*fs,'Color',[0 0 1]);
   ylabel(hAx(2),'U  [eV]','fontsize',1*fs,'Color',[0 0 0]);
   set(gca,'fontsize',fs);
   box on; grid on;  
   
   subplot(3,1,3) % ------------------------------------------------------
     [hAx,hLine1,hLine2] = plotyy(x,yR.^2 + yI.^2,x,U);
   
     set(hAx(1),'yLim',[yLL1 yLL2]);
     set(hAx(1),'yTick',yLL1 : yLL : yLL2);
     set(hAx(1),'fontsize',fs,'YColor',[0 0 1]);
     set(hLine1,'linewidth',2,'Color',[0 0 1]);
      
     set(hAx(2),'yLim',[U1 U2]);
    if U0 == 0 
      set(hAx(2),'yTick',0);
    else
      set(hAx(2),'yTick',UyTick);
    end 
    set(hAx(2),'fontsize',fs,'YColor',[0 0 0]);
    set(hLine2,'linewidth',1,'Color',[0 0 0]);
   
    xlabel('x  [m]','fontsize',fs);
    ylabel(hAx(1),'\psi^* . \psi','fontsize',2*fs,'Color',[0 0 1]);
    ylabel(hAx(2),'U  [eV]','fontsize',1*fs,'Color',[0 0 0]);
    set(gca,'fontsize',fs);
    box on; grid on;
   
   pause(0.001);
   drawnow;
        
      if flag2 > 0
         frame1 = frame1 + 1;
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
      % On the first loop, create the file. In subsequent loops, append.
         if frame1 == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
      end 
      
   end
  end
end


% ========================================================================
% EXPECTATION VALUES
% ========================================================================
    C3 = -hbar^2 / (2 * me * dx^2 * e);
% PROBABILITY
    fn = (yR - 1i*yI) .* (yR + 1i*yI);
    prob = simpson1d(fn,0,L);
    
    prob1 = simpson1d(fn(1:NxC),0,L/2);
    prob2 = simpson1d(fn(NxC:end),L/2,L);

% POSITION x
   fn = x .* (yR.^2 + yI.^2);
   xavg = simpson1d(fn,0,L) ;
   
% Potential energy  U
    fn = U .* (yR.^2 + yI.^2);
    Uavg = simpson1d(fn,0,L);

% Kinetic energy KE
    for nx = 2 : Nx-1
         fn(nx) = C3 * (yR(nx) - 1i * yI(nx)) * ...
         (yR(nx+1) - 2 * yR(nx) + yR(nx-1)+ 1i *(yI(nx+1) - 2 * yI(nx) + yI(nx-1)));
    end
    Kavg = simpson1d(fn,0,L);
 
% Total energy   
    Eavg = Uavg + Kavg;
 
% Momentum  p and speed v
 dyRdx = zeros(1,Nx); dyIdx = zeros(1,Nx);
 dyRdx(1) = (yR(2) - yR(1))/dx; dyRdx(Nx) = (yR(Nx) - yR(Nx-1))/dx;
 dyIdx(1) = (yI(2) - yI(1))/dx; dyIdx(Nx) = (yI(Nx) - yI(Nx-1))/dx;
 for nx = 2 : Nx-1
    dyRdx(nx) = (yR(nx+1) - yR(nx-1))/(2*dx);
    dyIdx(nx) = (yI(nx+1) - yI(nx-1))/(2*dx);
 end
   fn = - hbar .* (yI .* dyRdx - yR .* dyIdx + ...
         1i  .* yR .* dyRdx + yI .* dyIdx );
   pavg = simpson1d(fn,0,L) ; 
 
   vavg = pavg / me;
   
% Uncertainty - position
   fn = x .* x .* (yR.^2 + yI.^2);
   x2avg = simpson1d(fn,0,L); 
   delta_x = sqrt(x2avg - xavg^2);

% Uncertainty - momentum
    C3 = -hbar^2/dx^2;
    for nx = 2 : Nx-1
        fn(nx) = C3 * (yR(nx) - 1i * yI(nx)) * ...
        (yR(nx+1) - 2 * yR(nx) + yR(nx-1)+ 1i *(yI(nx+1) - 2 * yI(nx) + yI(nx-1)));
    end
    p2avg = simpson1d(fn,0,L);
    delta_p = sqrt(p2avg - pavg^2); 
 
% Uncertainty Principle
    dxdp = delta_x * delta_p;
    
    
% ========================================================================   
% Graphics   
% ========================================================================

figure(2)
   set(gcf,'units','normalized','position',[0.1 0.1 0.42 0.8]);
   
   subplot(3,1,1) % ------------------------------------------------------
   xP = x; yP = yR1;
   plot(xP,yP,'r','linewidth',1);
   hold on
   
   [hAx,hLine1,hLine2] = plotyy(x,yR,x,U);
   
   set(hAx(1),'yLim',[yL1 yL2]);
   set(hAx(1),'yTick',yL1 : yL : yL2);
   set(hAx(1),'fontsize',fs,'YColor',[0 0 1]);
   set(hLine1,'linewidth',2,'Color',[0 0 1]);
      
   set(hAx(2),'yLim',[U1 U2]);
   if U0 == 0 
      set(hAx(2),'yTick',0);
   else
      set(hAx(2),'yTick',UyTick);
   end 
   set(hAx(2),'fontsize',fs,'YColor',[0 0 0]);
   set(hLine2,'linewidth',1,'Color',[0 0 0]);
   
   xlabel('x  [m]','fontsize',fs);
   ylabel(hAx(1),'\psi_{real}','fontsize',2*fs,'Color',[0 0 1]);
   ylabel(hAx(2),'U  [eV]','fontsize',1*fs,'Color',[0 0 0]);
   
   tm1 = 't_{steps} = ';
   tm2 = num2str(nt,'%5.0f\n');
   tm3 = '  <E> = ';
   tm4 = num2str(real(Eavg),'%3.0f\n');
   tm5 = '  <K> = ';
   tm6 = num2str(real(Kavg),'%3.0f\n');
   tm7 = '  <U> = ';
   tm8 = num2str(real(Uavg),'%3.0f\n');
   tm9 = '   eV';
   
   tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9];
   h_title = title(tm,'fontsize',fs);
   set(gca,'fontsize',fs);
   box on; grid on;
   
 subplot(3,1,2) % ------------------------------------------------------
   xP = x; yP = yI1;
   plot(xP,yP,'r','linewidth',1);
   hold on
   
   [hAx,hLine1,hLine2] = plotyy(x,yI,x,U);
   
   set(hAx(1),'yLim',[yL1 yL2]);
   set(hAx(1),'yTick',yL1 : yL : yL2);
   set(hAx(1),'fontsize',fs,'YColor',[0 0 1]);
   set(hLine1,'linewidth',2,'Color',[0 0 1]);
      
   set(hAx(2),'yLim',[U1 U2]);
   if U0 == 0 
      set(hAx(2),'yTick',0);
   else
      set(hAx(2),'yTick',UyTick);
   end 
   set(hAx(2),'fontsize',fs,'YColor',[0 0 0]);
   set(hLine2,'linewidth',1,'Color',[0 0 0]);
   
   xlabel('x  [m]','fontsize',fs);
   ylabel(hAx(1),'\psi_{imag}','fontsize',2*fs,'Color',[0 0 1]);
   ylabel(hAx(2),'U  [eV]','fontsize',1*fs,'Color',[0 0 0]);
   set(gca,'fontsize',fs);
   box on; grid on;  
   
   subplot(3,1,3) % ------------------------------------------------------
   xP = x; yP = yP11;
   plot(xP,yP,'r','linewidth',1);
   hold on
   
   [hAx,hLine1,hLine2] = plotyy(x,yR.^2 + yI.^2,x,U);
   
   set(hAx(1),'yLim',[yLL1 yLL2]);
   set(hAx(1),'yTick',yLL1 : yLL : yLL2);
   set(hAx(1),'fontsize',fs,'YColor',[0 0 1]);
   set(hLine1,'linewidth',2,'Color',[0 0 1]);
      
   set(hAx(2),'yLim',[U1 U2]);
   if U0 == 0 
      set(hAx(2),'yTick',0);
   else
      set(hAx(2),'yTick',UyTick);
   end 
   set(hAx(2),'fontsize',fs,'YColor',[0 0 0]);
   set(hLine2,'linewidth',1,'Color',[0 0 0]);
   
   xlabel('x  [m]','fontsize',fs);
   ylabel(hAx(1),'\psi^* . \psi','fontsize',2*fs,'Color',[0 0 1]);
   ylabel(hAx(2),'U  [eV]','fontsize',1*fs,'Color',[0 0 0]);
   set(gca,'fontsize',fs);
   box on; grid on;


figure(3) %-----------------------------------------------------------
   fs = 14;
   set(gcf,'units','normalized','position',[0.54 0.1 0.4 0.8]);
   axis([0 100 0 100]);
   
   tm = 'Gaussian Wave Packet';
   d = 6; tx = 0; ty = 106;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'wavelength,   \lambda  =  ';
   tm2 = num2str(wL,'%2.2e\n');
   tm3 = '   m';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'kinetic energy  KE = h^2 / (2 \lambda^2 m_e) =  ';
   tm2 = num2str(KE,'%3.2f\n');
   tm3 = ' eV';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs)
   
   tm1 = 'initial kinetic energy  <KE_{initial}>  =  ';
   tm2 = num2str(real(K1avg),'%3.2f\n');
   tm3 = '   eV';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'width,   s  =  ';
   tm2 = num2str(s,'%2.2e\n');
   tm3 = '   m';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'spatial grid points,  N_x  =  ';
   tm2 = num2str(Nx,'%2.0f\n');
   tm = [tm1 tm2];
   d = 6; tx = 10; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'Evolution time:   time steps,  N_t  =  ';
   tm2 = num2str(Nt,'%2.0f\n');
   tm3 = '     time, t  =  ';
   tm4 = num2str(Nt*dt,'%2.2e\n');
   tm5 = '  s';
   tm = [tm1 tm2 tm3 tm4 tm5];
   tx = 0; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'probability,  Prob  =  ';
   tm2 = num2str(prob,'%2.4f\n');
   tm = [tm1 tm2];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'probability left side,  Prob_1  =  ';
   tm2 = num2str(prob1,'%2.4f\n');
   tm = [tm1 tm2];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'probability right side,  Prob_2  =  ';
   tm2 = num2str(prob2,'%2.4f\n');
   tm = [tm1 tm2];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
    
   tm1 = 'position,  <x>  =  ';
   tm2 = num2str(xavg,'%2.4e\n');
   tm = [tm1 tm2];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'potential energy,  <U>  =  ';
   tm2 = num2str(real(Uavg),'%3.2f\n');
   tm3 = '   eV';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'kinetic energy,  <K>  =  ';
   tm2 = num2str(real(Kavg),'%3.2f\n');
   tm3 = '   eV';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'total energy,  <E>  =  ';
   tm2 = num2str(real(Kavg)+Uavg,'%3.2f\n');
   tm3 = '   eV';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'momentum,  <p>  =  ';
   tm2 = num2str(real(pavg),'%3.3e\n');
   tm3 = '   kg.m/s';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'speed,  <v>  =  ';
   tm2 = num2str(real(vavg),'%3.3e\n');
   tm3 = '   m/s';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm = 'Heisenberg Uncertainty Principle   \Deltax . \Delta p > hbar / 2 = 0.53e-34 ';
   tx = 0; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '\Deltax  =  ';
   tm2 = num2str(delta_x,'%3.2e\n');
   tm3 = '   m';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '\Deltap  =  ';
   tm2 = num2str(real(delta_p),'%3.2e');
   tm3 = '   kg.m/s';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '\Deltax . \Deltap  =  ';
   tm2 = num2str(real(delta_x * delta_p),'%3.2e');
   tm3 = '   kg.m^2/s';
   tm = [tm1 tm2 tm3];
   tx = 10; ty = ty - d;
   text(tx,ty,tm,'fontsize',fs);
   
   axis off
   
   figure(4)  % ----------------------------------------------------------
   set(gcf,'units','normalized','position',[0.05 0.05 0.42 0.80]);
   set(gcf,'Color',[1 1 1]);
   for cc = 1 : 9
     subplot(3,3,cc) % ------------------------------------------------------
     xP = x; yP1 = probD(:,cc); yP2 = U;
     [hAx,hLine1,hLine2] = plotyy(x,yP1,x,yP2);
     set(hLine1,'linewidth',2,'Color',[0 0 1]);
     set(hAx(1),'yLim',[yLL1 yLL2]);
     set(hAx(1),'yTick',[]);
     set(hAx(2),'yLim',[U1 U2]);
     set(hAx(2),'yTick',[]); 
     set(hAx(2),'fontsize',fs,'YColor',[1 1 1]);
     set(hLine2,'linewidth',1,'Color',[0 0 0]);
     set(gca,'fontsize',fs);
     box on;
     set(hAx(2),'box','off');
     axis off
   end

   
   % =====================================================================
   
   
%   w
%  100*prob2
   
   toc