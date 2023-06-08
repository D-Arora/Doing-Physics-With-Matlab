% QMGAE.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMGAE.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230531   Matlab R2021b

% Solving the [1D] Time Dependent Schrodinger Equation using the FDTD Method
% The motion of a free electron in a uniform electric field
% Expectation values: probability density, position, momentum, velocity,
%                     kinetic energy, potential energy, total energy


% http://eng-web1.eng.famu.fsu.edu/~dommelen/quantum/style_a/packets.html
% https://www.mdpi.com/2073-8994/14/10/2215

close all
clear 
clc

global e hbar me

tic


% INPUT PARAMETERS   ==================================================
  Nx = 999;                % Number of grid points [999 must be an ODD]
  Nt = 22200;              % Number of time steps [22200]
  L = 10e-9;               % Width of X domain [10e-9 m] 
% Wavepacket and potential energy parameters
   nx1 = round(0.6*Nx);    % Pulse centre index [round(0.6*Nx]
   s = L/25;               % pulse width [L/20]
   wL0 = 3e-10;            % nominal wavelength [m]
   U0 = -100;              % Potential energy function [-100 eV]
   
% CONSTANTS   ========================================================
me = 9.10938291e-31;     % electron mass
h    = 6.62607015e-34;   % Planck's constant
hbar = 1.054571726e-34;  % Reduced Planck's constant
e = 1.602176565e-19;     % elementary charge

% ====================================================================   
% Setup for saving images animated gif file and AVI files

% Save animation; flagA = 1 YES / flagA = 0 NO
     flagA1 = 1; 
%  Enter file name
     ag_name1 = 'agA1.gif';
%  Delay in seconds before displaying the next image  
    delay1 = 0.; delay2 = 0.4;  
%  Frame counter start
    cs = 1;
%  ====================================================================

% SETUP  ==============================================================
  x = linspace(0,L,Nx); dx = x(2) - x(1);  % X grid
  k1 = -hbar^2/(2*me);
  C1 = 1/10;                               % Stability constant
  dt = C1 * 2 * me * dx^2 / hbar;          % Time step
  t = 0:dt:(Nt-1)*dt;                      % Time grid
  C2 = e*dt/hbar;
  C3 = -hbar^2 / (2 * me * dx^2 * e);
  
% Initialize arrays: Wavefunction y yR yI / prob density pd  
  psiR = zeros(Nt,Nx);      % real part wavefunction
  psiI = zeros(Nt,Nx);      % imaginary part wvefunction
  psi = psiR + 1i.*psiI;    % complex wavefunction
  pd = zeros(Nt,Nx);        % probability density
  phase = zeros(Nt,Nx);     % phase

  Prob = zeros(1,Nt);       % probability
  Uavg = zeros(1,Nt);       % expectation value: potential energy
  Kavg = zeros(1,Nt);       % expectation value: kinetic energy
  Eavg = zeros(1,Nt);       % total energy
  xavg = zeros(1,Nt);       % expectation value: position
  pavg = zeros(1,Nt);       % momentum
  vavg = zeros(1,Nt);       % velocity
  

% POTENTIAL ENERGY FUNCTION  ======================================
  U = -(U0/L).*x + U0;

% INITIAL WAVE PACKET (t = 0) =====================================
  xC = x(nx1);     % centre of wavepacket at t = 0
  yR = exp(-0.5.*((x-xC)./s).^2).*cos((2*pi.*(x-xC))./wL0);
  yI = exp(-0.5.*((x-xC)./s).^2).*sin((2*pi.*(x-xC))./wL0);
% Normalize wavefunction  
   [yR, yI,pd(1,:)] = normalize(yR,yI,L);
   psiNorm1 = simpson1d(pd(1,:),0,L);
   psiR(1,:) = yR;
   psiI(1,:) = yI;
    

% Solve Schrodinger Equation: FDTD Method  ============================
for nt = 1 : Nt-1
   for nx = 2 : Nx - 1
      yR(nx) = yR(nx) - C1*(yI(nx+1)-2*yI(nx)+yI(nx-1)) + C2*U(nx)*yI(nx);
    end
   psiR(nt+1,:) = yR;
   
   for nx = 2 : Nx-1
      yI(nx) = yI(nx) + C1*(yR(nx+1)-2*yR(nx)+yR(nx-1)) - C2*U(nx)*yR(nx);
   end
    psiI(nt+1,:) = yI;
end 

for nt = 2:Nt
  psiI(nt,:) = 0.5*(psiI(nt,:) + psiI(nt-1,:));
end

% EXPECTATION VALUES  ==============================================
for nt = 1:Nt
    yR = psiR(nt,:); yI = psiI(nt,:);
    pd(nt,:) = yR.^2 + yI.^2;
    phase(nt,:) = angle(yR + 1i.*yI);
    [Prob(nt),xavg(nt),Uavg(nt),Kavg(nt),Eavg(nt),pavg(nt),vavg(nt)]...
     = exValues(yR,yI,U,x,dx,L);    
end

% UNCERTAINITIES and HEISENBERG UNCERTAINITY PRINCIPLE  ==============
% At t = 0  (initial values)
   yR = psiR(1,:); yI = psiI(1,:);
   [deltaX1, deltaP1, deltaXP1] = exValuesUC(yR,yI,x,dx,L);
% At t = Nt * dt (final values)
   yR = psiR(Nt,:); yI = psiI(Nt,:);
   [deltaX2, deltaP2, deltaXP2] = exValuesUC(yR,yI,x,dx,L);   

   DX = zeros(Nt,1); DP = DX; DXP = DX;
for nt = 1:Nt
    yR = psiR(nt,:); yI = psiI(nt,:);
    [DX(nt), DP(nt), DXP(nt)] = exValuesUC(yR,yI,x,dx,L);   
end
  
% Classical calculations from simulation results
%   nx2 = find(pd(Nt,:) == max(pd(Nt,:)));
%   v_class = (x(nx2)-x(nx1))/t(Nt);       % velocity [m/s]
%   p_class = me*v_class;                  % momentum [N.s]                   
%   K_class = 0.5*me*v_class^2/e;          % K.E. [eV]

% Ehrenfest's Theorem
  F = U0*e/L;         % Classical: force acting on electron [N]
  a = F/me;           % Classical: acceleration of electron [m/s^2]
  p0 = h/wL0;         % Quantum: initial nominal momentum   [N.s]
  v1 = p0/me;         % Classical: initial velocity v(t = 0)  [m/s]
  x1 = xC;            % Initial centre position of wavepacket x(t = 0)  [m]
  v2 = v1 + a.*t(end);   % Classical: final velocity v(t = end)  [m/s]
  X2 = x1 + v1*t(end) + 0.5*a*t(end)^2;   % Classical final position x(t = end) [m]
  p1 = me.*v1;         % Classical; momentum p(t = 0)  [N.s]
  p2 = me.*v2;         % Classical: momentum p(t = end)  [N.s]
  K1 = p1^2./(2*me*e);   % Classical: kinetic energy K(t = 0)  [eV]
  K2 = p2^2./(2*me*e);   % Classical: kinetic energy K(t = end)  [eV]
  W = K2 - K1;           % Classical: work done on electron by field [eV]
  Jp = p2 - p1;          % Classical: impulse = change in momentum  [N.s]
  JF = F*t(end);         % Classical: impulse = force * time [N.s]
  W = K2 - K1;           % Classical: work = change in KE  [eV]
  wL1 = wL0;             % Initial wavelength  [m]
  wL2 = abs(h/p2);            % Final wavelength [m]

  fit1 = polyfit(t,vavg,1);   % Fit v = v + at
  a1 = fit1(1);       % <a>     
  fit2 = polyfit(t,xavg,2);   % Fit x = 0.5 a t^2 + u t + x0
  a2 = 2*fit2(1);      % <a>    
  v02 = fit2(2);      % initial velocity 




%%

% =====================================================================
% GRAPHICS
% =====================================================================
  FS = 12;
  
figure(1)   % Animation ----------------------------------------------
  pos = [0.05 0.05 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
  xP = x.*1e9; 
  YLIM = [-0.6e5 0.6e5];
    
for cc = 1:200: Nt
subplot(3,1,1)
  yP = psiR(cc,:);
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = psiR(1,:);
  plot(xP,yP,'k','linewidth',0.1)
  hold off
  grid on
  box on
  set(gca,'fontsize',FS)
  ylim(YLIM)
  ylabel('\psi_{real}')
subplot(3,1,2)
  yP = psiI(cc,:);
  plot(xP,yP,'r','linewidth',2)
  hold on
  yP = psiI(1,:);
  plot(xP,yP,'k','linewidth',0.1)
  hold off
  grid on
  box on
  set(gca,'fontsize',FS)
  ylim(YLIM)
  ylabel('\psi_{imag}')
subplot(3,1,3)
  yP = psiR(cc,:).^2 + psiI(cc,:).^2;
  plot(xP,yP,'k','linewidth',2)
  hold on
  yP = psiR(1,:).^2 +psiI(1,:).^2;
  plot(xP,yP,'k','linewidth',1)
  hold off
  ylabel('Prob. Density [ nm^{-1} ]')
  grid on
  box on
  set(gca,'fontsize',FS)
  ylim([0 2e9])
  xlabel('x  [ nm ]')
  hold off
  pause(0.00001)
     
   if flagA1 == 1
     frame = getframe(1);
     im = frame2im(frame);
     [imind,cm] = rgb2ind(im,256);
%  On the first loop, create the file. In subsequent loops, append.
     if cs == 1
        imwrite(imind,cm,ag_name1,'gif','DelayTime',delay1,'loopcount',inf);
     else
        imwrite(imind,cm,ag_name1,'gif','DelayTime',delay1,'writemode','append');
     end
        cs = cs+1;
   end
  
end

%%
% PARAMETER LIST -----------------------------------------------------
  figure(2) 
    fs = 14;
    set(gcf,'units','normalized','position',[0.41 0.05 0.4 0.7]);
    axis([0 100 -20 100]);
   
    tm = 'Simulation parameters';
    tx = 1; ty = 106;
    text(tx,ty,tm,'fontsize',fs,'color','b');
    
    tm = sprintf('Nx = %0.0f   Ny = %0.0f   L = %0.2e m \n',Nx,Nt,L);
    d = 6; tx = 4; ty = 106-d;
    text(tx,ty,tm,'fontsize',fs);

    tm = sprintf('Centre index = %0.0f   x_c = %0.2f nm ', nx1,x(nx1)/1e-9);
    d = 6; tx = 4; ty = ty-d;
    text(tx,ty,tm,'fontsize',fs);

    tm = sprintf('s = %0.2e m   \\lambda_0 = wL0 = %0.2e m   ',s,wL0);
    tx = 4; ty = ty-d;
    text(tx,ty,tm,'fontsize',fs);

    tm = 'Simulation: (1) initial values  (2) final values (C) computed values';
    tx = 1; ty = ty-d;
    text(tx,ty,tm,'fontsize',fs,'color','b');

    tm = sprintf('t_1 = %0.2f fs   t_2 = %0.2f  fs',t(1), t(end)/1e-15);
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs);

    tm = sprintf('<x_1> = %0.2f nm   <x_2> = %0.2f  nm   x_{2C} = %2.2f nm'...
        ,xavg(1), xavg(end), X2/1e-9);
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs);

    tm = sprintf('<v_1> = %0.2e m.s^{-1} <v_2> = %0.2e m.s^{-1} v_{2C} = %0.2e m.s^{-1}'...
        ,vavg(1), vavg(end),v2 );
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<a> = %0.2e m.s^{-2}    a_C = %0.2e m.s^{-2}',a1, a);
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<p_1> = %0.2e N.s   <p_2> = %0.2e N.s',pavg(1), pavg(end));
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<p_2> - <p_1> = %0.2e N.s   J_C = %0.2e N.s',-pavg(1)+pavg(end),Jp);
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<K_1> = %0.2f eV   <K_2> = %0.2f eV',Kavg(1), Kavg(end));
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<K_2> - <K_1> = %0.2f eV   W = %0.2f eV',-Kavg(1)+Kavg(end),W);
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)
    
    tm = sprintf('<U_1> = %0.2f eV   <U_2> = %0.2f eV',Uavg(1), Uavg(end));
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<U_2> - <U_1> = %0.2f eV',-Uavg(1)+Uavg(end));
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('<E_1> = %0.2f eV   <E_2> = %0.2f eV',Eavg(1), Eavg(end));
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    tm = sprintf('\\lambda_1 = %0.2e nm   \\lambda_2 = %0.2e nm',wL1, wL2);
    tx = 4; ty = ty-d-2;
    text(tx,ty,tm,'fontsize',fs)

    axis off
     

%%
figure(5)  % time evoltion plots: dx, x, a  --------------------- 
   FS = 14;
   set(gcf,'units','normalized','position',[0.3 0.05 0.25 0.7]);
   xP = t./1e-15;
subplot(4,1,3)
  yP = DX./1e-9';
  plot(xP,yP,'b','linewidth',2)
  xlim('tight')
  grid on
  box on
  ylabel('\Deltax  [ nm ]')
  set(gca,'fontsize',FS)
subplot(4,1,1)
  yP = xavg./1e-9;
  plot(xP,yP,'b','linewidth',2)
  xlim('tight')
  grid on
  box on
  ylabel('<x> [ nm ]')
  set(gca,'fontsize',FS)
subplot(4,1,2)
  yP = vavg;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim('tight')
  box on
  ylabel('<v>  [ m/s ] ')
  txt = sprintf('<a> = %0.2e', a1);
  title(txt,'FontWeight','normal')
  set(gca,'fontsize',FS)
subplot(4,1,4)
  yP = DP';
  plot(xP,yP,'b','linewidth',2)
  xlim('tight')
  grid on
  box on
  ylabel('\Deltap  [ N.s ]')
  xlabel('t  [ fs ]')
  set(gca,'fontsize',FS)

%%
figure(6)  % time evolution plots: E U K & U(x)  -----------------
   FS = 14;
   set(gcf,'units','normalized','position',[0.50 0.05 0.25 0.4]);
   xP = t./1e-15;
subplot('Position',[0.2 0.55, 0.7,0.4])
  yP = Uavg';
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = Kavg';
  plot(xP,yP,'r','linewidth',2)
  yP = Eavg';
  plot(xP,yP,'k','linewidth',2)
  legend('<U>','<K>','<E>','Orientation','horizontal','location',...
      'southwest','box','off')
  xlim('tight')
  grid on
  box on
  xlabel('time t [ fs ]')
  ylabel('energy  [ eV ]')
  set(gca,'fontsize',FS)

subplot('Position',[0.2 0.15, 0.7,0.2])
  xP = x./1e-9; yP = U;
  plot(xP,yP,'b','linewidth',2)
  xlim('tight')
  grid on
  box on
  xlabel('position x [ nm ]')
  ylabel('U [ eV ]')
  set(gca,'fontsize',FS)
 
% ==================================================================
toc
  

% ===================================================================
% FUNCTIONS
% ===================================================================
  
 function [f,g,pd] = normalize(F,G,L)
    M = F.^2 + G.^2;
    A = simpson1d(M,0,L);
    f = F ./ sqrt(A); g = G ./ sqrt(A);
    pd = f.^2 + g.^2;
 end


 function [Prob,xavg,Uavg,Kavg,Eavg,pavg,vavg] = exValues(pR,pI,U,x,dx,L)
 % EXPECTATION VALES
 global e hbar me
    y = pR + 1i.*pI;
    
    fn = conj(y).*y;
    Prob = simpson1d(fn,0,L);     % probability
  
    fn = conj(y).*x.*y;
    xavg = real(1e9*simpson1d(fn,0,L));  % position
   
    fn = conj(y).*U.*y;
    Uavg = real(simpson1d(fn,0,L));    % potential energy
  
    fn = 4*del2(y,dx);
    fn = conj(y).*fn;
    Kavg = (-hbar^2/(2*me*e))*real(simpson1d(fn,0,L));  % kinetic energy
  
    Eavg = Kavg + Uavg;   % total energy
  
    fn = gradient(y,dx);
    fn = conj(y).*fn;
    pavg =  real((-1i*hbar)*simpson1d(fn,0,L));  % momentum
    
    vavg = pavg/me;
 end
 

 function [deltaX, deltaP, deltaXP] = exValuesUC(yR,yI,x,dx,L)
 % EXPECTATION VALES: UNCERTAINITIES
 global hbar 
    y = yR + 1i.*yI;
    
    fn = conj(y).*x.*y;
    X1 = simpson1d(fn,0,L);    % position
   
    fn = conj(y).*x.*x.*y;
    X2 = simpson1d(fn,0,L);    % position^2
  
    fn = gradient(y,dx);
    fn = conj(y).*fn;
    P1 =  (-1i*hbar)*simpson1d(fn,0,L);  % momentum
  
   fn = 4*del2(y,dx);
   fn = conj(y).*fn;
   P2 =  hbar^2*simpson1d(fn,0,L);  % momentum^2
    
%     fn = gradient(y,dx);
%     fn = conj(y).*gradient(fn);
%     P2 =  hbar^2*simpson1d(fn,0,L);  % momentum^2
    
    deltaX = abs(sqrt(X2 - X1^2));
    deltaP = abs(sqrt(P2 - P1^2));
   
    deltaXP = deltaX*deltaP;
end 
 
 