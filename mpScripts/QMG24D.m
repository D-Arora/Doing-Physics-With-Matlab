% QMG24C.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMGAC.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230531   Matlab R2021b

% Solving the Time Dependent Schrodinger Equation using the FDTD Method
% The motion of a Guassian wave packet is animated.
% Expectation values: probability density, position, momentum, velocity,
%                     kinetic energy, potential energy, total energy

close all
clear 
clc

global e hbar me

tic

% ========================================================================
% INPUT PARAMETERS 
% ========================================================================
Nx = 701;                %  ODD number - number of grid points [701]
Nt = 7001;              %  Number of time steps [30000]
L = 5e-9;                %  Width of X domain [5e-9] 

me = 9.10938291e-31;     % electron mass
h    = 6.62607015e-34;   % Planck's constant
hbar = 1.054571726e-34;  % hbar Planck's constant
e = 1.602176565e-19;     % elementary charge

% =======================================================================   
% Setup for saving images animated gif file and AVI files
% =======================================================================
% Save animation; flagA = 1 YES / flagA = 0 NO
     flagA1 = 1; flagA2 = 0;
%  Enter file name
     ag_name1 = 'agA1.gif'; ag_name2 = 'agB.gif';
%  Delay in seconds before displaying the next image  
    delay1 = 0.; delay2 = 0.4;  
%  Frame counter start
    cs = 1;

% =======================================================================
% SETUP
% ======================================================================
  x = linspace(0,L,Nx); dx = x(2) - x(1);
  k1 = -hbar^2/(2*me);
  C1 = 1/10;
  dt = C1 * 2 * me * dx^2 / hbar;
  t = 0:dt:Nt*dt;
  C2 = e*dt/hbar;
  C3 = -hbar^2 / (2 * me * dx^2 * e);
  
% Initialize arrays  
% Wavefunction y yR yI / prob density pd  
  psiR = zeros(Nt,Nx);      % real part wavefunction
  psiI = zeros(Nt,Nx);      % imaginary part wvefunction
  psi = psiR + 1i.*psiI;    % complex wavefunction
  pd = zeros(Nt,Nx);        % probability density
  phase = zeros(Nt,Nx);     % phase

  U = zeros(1,Nx);          % potential energy
  Prob = zeros(1,Nt);       % probability
  Uavg = zeros(1,Nt);       % expectation value: potential energy
  Kavg = zeros(1,Nt);       % expectation value: kinetic energy
  Eavg = zeros(1,Nt);       % total energy
  xavg = zeros(1,Nt);       % expectation value: position
  pavg = zeros(1,Nt);       % momentum
  vavg = zeros(1,Nt);       % velocity
  

% ========================================================================
% INITIAL WAVE PACKET
% ========================================================================
   nx1 = round(Nx/5);     % pulse centre   [round(Nx/6)]
   s = L/25;              % pulse width    [L/20]
   wL0 = 2*1.5e-10;          % wavelength
  
  yR = exp(-0.5.*((x-x(nx1))./s).^2).*cos((2*pi).*(x-x(nx1))./wL0);
  yI = exp(-0.5.*((x-x(nx1))./s).^2).*sin((2*pi).*(x-x(nx1))./wL0);
 % yI = zeros(1,Nx);
% Normalize wavefunction  
   [yR, yI,pd(1,:)] = normalize(yR,yI,L);
   psiNorm1 = simpson1d(pd(1,:),0,L);
   psiR(1,:) = yR;
   psiI(1,:) = yI;
    
% =====================================================================
% Solve Schrodinger Equation: FDTD Method
% =====================================================================
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
% =====================================================================

for nt = 2:Nt
  psiI(nt,:) = 0.5*(psiI(nt,:) + psiI(nt-1,:));
end

% EXPECTATION VALUES
for nt = 1:Nt
    yR = psiR(nt,:); yI = psiI(nt,:);
    pd(nt,:) = yR.^2 + yI.^2;
    phase(nt,:) = angle(yR + 1i.*yI);
    [Prob(nt),xavg(nt),Uavg(nt),Kavg(nt),Eavg(nt),pavg(nt),vavg(nt)]...
     = exValues(yR,yI,U,x,dx,L);    
end

% UNCERTAINITIES and HEISENBERG UNCERTAINITY PRINCIPLE
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
  nx2 = find(pd(Nt,:) == max(pd(Nt,:)));
  v_class = (x(nx2)-x(nx1))/t(Nt);       % velocity [m/s]
  p_class = me*v_class;                  % momentum [N.s]                   
  K_class = 0.5*me*v_class^2/e;          % K.E. [eV]
    

% Phase and group velocities
  f0 = Eavg(1)*e/h;         % peak frequency E = hf  [Hz]
  vPhase = wL0*f0;          % Phase velocity [m/s]
  vGroup = pavg(1)/me;      % Group velocity [m/s]
  ratioGP = vGroup/vPhase;  % vGroup / vPhase

%%  
% Fourier transform PSI at t = 0   K = 1/wL
  xMin = 0; xMax = L; 
  KMax = 2/wL0; KMin= -KMax; nK = 2001;
  K = linspace(KMin,KMax,nK);
    hP = psiR(1,:); HP = zeros(1,nK);
   for c = 1:nK
     g = hP.* exp(1i*2*pi*K(c)*x);
     HP(c) = simpson1d(g,xMin,xMax);
   end
  %   HP = HP./max(HP);
     psd = conj(HP).*HP;

   hP = psiR(Nt,:); HP = zeros(1,nK);
   for c = 1:nK
     g = hP.* exp(1i*2*pi*K(c)*x);
     HP(c) = simpson1d(g,xMin,xMax);
   end
 %    HP = HP./max(HP);
     psd2 = conj(HP).*HP;  
    
%%

% =====================================================================
% GRAPHICS
% =====================================================================
  FS = 12;
  
figure(1)
  pos = [0.05 0.05 0.3 0.6];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
 % set(gcf,'color','w'); 
  xP = x.*1e9; 
  YLIM = [-0.6e5 0.6e5];
    
for cc = 1:100: Nt
subplot(3,1,1)
  yP = psiR(cc,:);
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = psiR(1,:);
  plot(xP,yP,'b','linewidth',1)
  hold off
  grid on
  box on
  set(gca,'fontsize',FS)
  ylim(YLIM)
  ylabel('\psi_{real}')

  tm = sprintf('t = %0.2f fs  <E> = %0.2f  <K> = %0.0f  <U> = %0.0f eV \n',...
      t(cc)/1e-15,Eavg(cc),Kavg(cc),Uavg(cc));
  h_title = title(tm,'fontsize',FS,'fontweight','normal');
  
subplot(3,1,2)
  yP = psiI(cc,:);
  plot(xP,yP,'r','linewidth',2)
  hold on
  yP = psiI(1,:);
  plot(xP,yP,'r','linewidth',1)
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
  ylim([0 3e9])
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
   set(gcf,'units','normalized','position',[0.41 0.05 0.35 0.85]);
   axis([0 100 -20 100]);
   
   tm = 'GAUSSIAN  WAVE  PACKET';
   tx = 15; ty = 106;
   text(tx,ty,tm,'fontsize',fs,'color','b');
   
   tm1 = 'wavelength  \lambda_0 = ';
   tm2 = num2str(1e9*wL0,'%2.2f\n');
   tm3 = ' nm';
   tm = [tm1 tm2 tm3];
   d = 5; tx = 4; ty = 106-d;
   text(tx,ty,tm,'fontsize',fs);
     
   tm1 = 'width s = ' ;
   tm2 = num2str(1e9*s,'%2.2f\n');
   tm3 = ' nm';
   tm4 = '   spatial grid points N_x = ';
   tm5 = num2str(Nx,'%2.0f\n');
   tm = [tm1 tm2 tm3 tm4 tm5];
   tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);

   d = 7;
   tm = sprintf('Expectation values');
   tx = 0; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs,'color','b');
   d = 3;
   tm1 = 'Initial values  t = 0 s';
   tm2 = '    Final values  t = ';
   tm3 = num2str(t(Nt),'%2.2e\n');
   tm4 = '  s';
   tm = [tm1 tm2 tm3 tm4];
   tx = 0; ty = ty-d-2;
   text(tx,ty,tm,'fontsize',fs);

   tm1 = 'Probability_i = ';
   tm2 = num2str(Prob(1),'%2.5f\n');
   tm = [tm1 tm2];
   d = 5; tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = 'Probability_f = ';
   tm2 = num2str(Prob(Nt),'%2.5f\n');
   tm = [tm1 tm2];
   tx = 46; 
   text(tx,ty,tm,'fontsize',fs);
 
   tm1 = '<K_i> = ';
   tm2 = num2str(Kavg(1),'%2.2f\n');
   tm3 = '  eV';
   tm = [tm1 tm2 tm3];
   tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<K_f> = ';
   tm2 = num2str(Kavg(Nt),'%2.2f\n');
   tm3 = '  eV';
   tm = [tm1 tm2 tm3];
   tx = 46; 
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<U_i> = ';
   tm2 = num2str(Uavg(1),'%2.2f\n');
   tm3 = '  eV';
   tm = [tm1 tm2 tm3];
   tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<U_f> = ';
   tm2 = num2str(Uavg(Nt),'%2.2f\n');
   tm3 = '  eV';
   tm = [tm1 tm2 tm3];
   tx = 46; 
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<E_i> = ';
   tm2 = num2str(Eavg(1),'%2.2f\n');
   tm3 = '  eV';
   tm = [tm1 tm2 tm3];
   tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<E_f> = ';
   tm2 = num2str(Eavg(Nt),'%2.2f\n');
   tm3 = '  eV';
   tm = [tm1 tm2 tm3];
   tx = 46; 
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<p_i> = ';
   tm2 = num2str(pavg(1),'%2.2e\n');
   tm3 = '  N.s';
   tm = [tm1 tm2 tm3];
   d = 5; tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<p_f> = ';
   tm2 = num2str(pavg(Nt),'%2.2e\n');
   tm3 = '  N.s';
   tm = [tm1 tm2 tm3];
   d = 6; tx = 46; 
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<v_i> = ';
   tm2 = num2str(vavg(1),'%2.2e\n');
   tm3 = '  m/s';
   tm = [tm1 tm2 tm3];
   tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
   tm1 = '<v_f> = ';
   tm2 = num2str(vavg(Nt),'%2.2e\n');
   tm3 = '  m/s';
   tm = [tm1 tm2 tm3];
   tx = 46; 
   text(tx,ty,tm,'fontsize',fs);
   
   tm = 'Heisenberg Uncertainty Principle: (\Deltax.\Deltap) / hbar >= 0.5';
   tx = 0; ty = ty-d-2;
   text(tx,ty,tm,'fontsize',fs,'color','b');
   
   tm = ' Initial values                      Final values';
   tx = 0; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs,'color','b');
    
   tm1 = '  \Deltax = ';
   tm2 = num2str(1e9*deltaX1,'%2.2f\n');
   tm3 = '  nm';
   tm = [tm1 tm2 tm3];
   tx = 0; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs)
   
   tm1 = '\Deltax = ';
   tm2 = num2str(1e9*deltaX2,'%2.2f\n');
   tm3 = '  nm';
   tm = [tm1 tm2 tm3];
   tx = 50; 
   text(tx,ty,tm,'fontsize',fs)
   
   tm1 = '  \Deltap = ';
   tm2 = num2str(deltaP1,'%2.2e\n');
   tm3 = '  N.s';
   tm = [tm1 tm2 tm3];
   tx = 0; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs)
   
   tm1 = '\Deltap = ';
   tm2 = num2str(deltaP2,'%2.2e\n');
   tm3 = '  N.s';
   tm = [tm1 tm2 tm3];
   tx = 50; 
   text(tx,ty,tm,'fontsize',fs)
   
   tm1 = '  (\Deltax.\Deltap) / hbar = ';
   tm2 = num2str(deltaX1*deltaP1/hbar,'%2.2f\n');
   tm = [tm1 tm2];
   tx = 0; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs)
   
   tm1 = '(\Deltax.\Deltap)/hbar = ';
   tm2 = num2str(deltaX2*deltaP2/hbar,'%2.2f\n');
   tm = [tm1 tm2];
   tx = 50; 
   text(tx,ty,tm,'fontsize',fs)
   
   tm = sprintf('Classical values: v = \\Deltax / \\Deltat   p = m v   K = p^2 / 2m');
   tx = 0; ty = ty-d-3;
   text(tx,ty,tm,'fontsize',fs,'color','b')
   tm = sprintf('v = %2.2e  m/s   p = %2.2e N.s   K = %2.2f eV \n', ...
        v_class,p_class,K_class);
   d = 7; tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   tm = sprintf('v_{Group} = %2.2e  m/s \n',vGroup);
   d = 6; tx = 4; ty = ty-d;
   text(tx,ty,tm,'fontsize',fs);
   
    tm = sprintf('Quantum values  p = h / \\lambda   E = hf   v = f \\lambda \n');
    tx = 0; ty = ty-d-3;
    text(tx,ty,tm,'fontsize',fs,'color','b')
    tm = sprintf('f_0 = %2.2e Hz   \\lambda_0 = %2.2f nm \n',f0,wL0/1e-9);
    d = 6; tx = 4; ty = ty-d;
    text(tx,ty,tm,'fontsize',fs);

    tm = sprintf('v_{Phase} = %2.2e m/s      v_{Phase} / v_{Group} = %2.2f \n',vPhase,vPhase/vGroup);
    d = 6; tx = 4; ty = ty-d;
    text(tx,ty,tm,'fontsize',fs);
   
    axis off
   
  

%%
figure(5) 
   FS = 14;
   set(gcf,'units','normalized','position',[0.51 0.05 0.3 0.56]);
   xP = t(1:end-1)./1e-15;
subplot(3,1,1)
  yP = DX./1e-9';
  plot(xP,yP,'b','linewidth',2)
  xlim('tight')
  grid on
  box on
  ylabel('\Deltax  [ nm ]')
  xlabel('t  [ fs ]')
  set(gca,'fontsize',FS)
subplot(3,1,2)
  yP = DP;
  plot(xP,yP,'b','linewidth',2)
  ylim([5 10].*1e-24)
  xlim('tight')
  grid on
  box on
  ylabel('\Deltap  [ N.s ]')
  xlabel('t  [ fs ]')
  set(gca,'fontsize',FS)
subplot(3,1,3)
  yP = DXP./hbar;
  plot(xP,yP,'b','linewidth',2)
  hold on
  plot([xP(1) xP(end)],[0.5 0.5],'r','linewidth',2)
  grid on
  xlim('tight')
  ylim([0 20])
  box on
  ylabel('\Deltax \Deltap / hbar  ')
  xlabel('t  [ fs ]')
  set(gca,'fontsize',FS)
   
 figure(3) % WARNING:Takes a long time to plot the graph, don't know why?
 % Color mapping for phase phi using script ColorCode.m
   m = (380-660) / (2*pi);
   b = 380 - m*pi;
   cs = 1;
 
   pos = [0.5 0.05 0.350 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   box on
   hold on
   xlabel('x  [nm]')
   ylabel('Prob. density  [nm_{-1}]')
  
  for ct = [1 Nt] %1 : 400: Nt
    % clf
    % hold off 
     indexcx = find(pd(ct,:) > pd(ct,:)/20);
     for cx = indexcx(1:end-1) %= 1:2:Nx-1
      hold on
      wLc = (m * phase(ct,cx) + b)*1e-9;
      thisColor = ColorCode(wLc);
      xP = 1e9*x(cx:cx+1); yP = pd(ct,cx:cx+1);
      h_area = area(xP,yP);
      set(h_area,'FaceColor',thisColor);
      set(h_area,'EdgeColor',thisColor);
      ylim([0 3e9])
      xlim([0 1e9*L])
      set(gca,'fontsize',FS)
      xlabel('x  [nm]')
      ylabel('Prob. density  [nm^{-1}]')
      grid on
      box on
     % hold off
     end
  %pause(0.0001)  
  %   hold off
   
%    if flagA2 == 1
%      frame = getframe(3);
%      im = frame2im(frame);
%      [imind,cm] = rgb2ind(im,256);
% %  On the first loop, create the file. In subsequent loops, append.
%      if cs == 1
%         imwrite(imind,cm,ag_name2,'gif','DelayTime',delay2,'loopcount',inf);
%      else
%         imwrite(imind,cm,ag_name2,'gif','DelayTime',delay2,'writemode','append');
%      end
%         cs = cs+1;
%    end
%    hold off
  end
%  --------------------------------------------------------------------  

%%
figure(6)    % FOURIER TRANSFORM
set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.1 0.25 0.20]);
  set(gcf,'color','w');
  FS = 14;

  xP = h.*K; yP = psd;
  plot(xP,yP,'b','LineWidth',4)
  hold on
  yP = psd2;
  plot(xP,yP,'r','LineWidth',2)
  xlim([0,h*KMax])
  grid on; box on
  xlabel('p  [ N.s ]'); ylabel('psd')
  set(gca,'FontSize',FS)
  pPeak = max(h*K(psd==max(psd)));
  txt = sprintf('        p_{peak} = %0.4e   ',pPeak);
  title(txt,'FontWeight','normal')
%%  
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
    Uavg = real(simpson1d(fn,0,L))/e;    % potential energy
  
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
 
 