% atmDDP.m

% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% IAN COOPER
%    matlabvisualphysics@gmail.com
% LINK: DOWNLOAD MATLAB SCRIPTS 
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb
% MATLAB VERSIOB R2021b
% 211204

% Damped Driven Pendulum  DDP
% SI units are for all variables
%  theta - angular displacement / omega - angular frequency / gamma - drive strength
%  ode45 solves the equation of motion

%  Ouput:  Phase space plots
%          
% Supporting documentation
%  https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/atmDDP.htm 




clear
close all
clc


tic
% SI units

% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Time 
  tMax = 60;
  tStart = 40;
  N = 9999;
  if tStart > tMax; tStart = 0; end

% Drive strength
  gamma = 1.078;

% Initial conditions: displacement X = theta / w = velocity
  X0 = -pi/2;
  w0 = 0;

% Flag: Calculate and plot Fourier Transform  flagFT = 1
  flagFT = 1;
% Display outpu parameters in Command Window  flagP = 1;
  flagP = 1;

% SETUP ====================================================
% Time  
  tMin = 0;
  tSpan = [tMin tMax];
  %tSpan = linspace(tMin,tMax,N);
 
% Driving strength / period / angular frequency
  TD = 1;
  wD = 2*pi/TD;

% Natutal frequency and period
  wN = 1.5*wD; TN = 2*pi/wN;

% Damping   
  beta = wN/4;


% Initial conditions vector
  s0 = [X0 w0];

% Coefficient vector
  K(1) = wN^2;
  K(2) = 2*beta;
  K(3) = gamma*wN^2;
  K(4) = wD;
  

% SOLVE ODE  ======================================================  
% Solve ODE
  opts = odeset('RelTol',1e-10);  
  [t, sol] = ode45(@(t,s) EM(t,s,K), tSpan,s0); 

% Angular displacement and angular velocity
  X = sol(:,1); w = sol(:,2);

% Restrict angular displacment -pi to + pi
  % X = atan2(tan(X),1);

% Frequencies  natural / driving
%  F = K(3).*cos(K(4)*t);
  fN = 1/TN;
  fD = 1/TD;

% Find oscillation period from peaks 
  [pks, locs] = findpeaks(X,t);
  TP = (locs(end) - locs(end-2))/2;
 
if flagFT ==1
    
  % FOURIER TRANSFORM - angulat displacement (theta) =================
   nT = length(t);
   nF = 2001;
   fMax = 5; fMin = -fMax;
   f = linspace(fMin,fMax, nF);

   H = zeros(1,nF);
   h = X';

   for c = 1:nF
     g = h.* exp(1i*2*pi*f(c)*t');
     H(c) = simpson1d(g,tMin,tMax);
   end
   % One-sided power spectral density PSD  Ph(f)
   Ph = 2.*conj(H).*H;
end


% OUPUTS ===========================================================
if flagP == 1
   txt = sprintf('theta(0) = %2.3f  \n',X0);
    disp(txt) 
  txt = sprintf('omega(0) = %2.3f  \n',w0);
    disp(txt) 
  txt = sprintf('Driving period TD = %2.3f  \n',TD);
    disp(txt);
  txt = sprintf('Driving freq wD = %2.3f  \n',wD);
    disp(txt)
  txt = sprintf('Natural period TN = %2.3f  \n',TN);
    disp(txt)
  txt = sprintf('Natural freq  wN = %2.3f  \n',wN);
    disp(txt)
  if flagP == 1
  txt = sprintf('period from peaks  TP = %2.3f  \n',TP);
    disp(txt)
  end
  txt = sprintf('beta = %2.3f  \n',beta);
    disp(txt) 
  txt = sprintf('gamma = %2.3f  \n',gamma);
    disp(txt)   
end

% GRAPHICS  ===========================================================

figure(1)   % drive signal
   txtG = sprintf('\\gamma = %2.5f   \n',gamma);
   pos = [0.05 0.5 0.35 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

   xP = 0:0.01:10; yP = cos(2*pi.*xP);
   plot(xP,yP,'r','linewidth',2)
   grid on
   ylabel('F   [ rad.s^{-2} ]')
   xlabel('t  [ s ]')
   title(txtG,'fontweight','normal')
   set(gca,'fontsize',14)


figure(2)   % angular displacement  X = theta
   pos = [0.05 0.05 0.35 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

   xP = t; yP = X./pi;
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('\theta   [ rad / \pi ]')
   xlabel('t  [ s ]')
   title(txtG,'fontweight','normal')
   set(gca,'fontsize',14)

figure(3)   % angular frequency  w = omega
   pos = [0.45 0.05 0.35 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

   xP = t; yP = w;
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('\omega   [ rad.s^{-1} ]')
   xlabel('t  [ s ]')
   title(txtG,'fontweight','normal')
   set(gca,'fontsize',14)

figure(4)   % phase space
   pos = [0.45 0.5 0.35 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

  
   nStart = find(t>=tStart,1);
   nEnd = length(t);
   nR = nStart:nEnd;

   xP = X(nR)./pi; yP = w(nR);
   plot(xP,yP,'b','linewidth',2)
   hold on
   HP = plot(xP(1),yP(1),'go');
   set(HP,'markerfacecolor','g','markersize',10)
%   xlim([-1.1 1.1])
%   xticks(-1:0.25:1)
   xtickformat('%2.2f')
   grid on
   xlabel('\theta   [ rad / \pi ]')
   ylabel('\omega   [ rad.s^{-1} ]')
   title(txtG,'fontweight','normal')
   set(gca,'fontsize',14)


if flagFT == 1
figure(5)   % Fourier Transform
   pos = [0.3 0.3 0.35 0.40];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on

  

subplot(2,1,1)
   xP = f; yP = Ph./max(Ph);
   plot(xP,yP,'b','lineWidth',2);
   hold on
   xP = [fN fN]; yP = [0 1];
   plot(xP,yP,'r','lineWidth',2);
   xP = [fD fD]; yP = [0 1];
   plot(xP,yP,'m','lineWidth',2);
  % title('One-sided PSD');
   ylabel('PSD  [ a.u.]');
   xlabel('f  [ Hz ]');
   grid on 
   txt = sprintf('F_N = %2.2f Hz  \n',fN);
   hT = text(3.7, 0.7,txt);
   set(hT,'fontsize',14,'color','r')
   txt = sprintf('F_D = %2.2f Hz  \n',fD);
   hT = text(3.7, 0.2,txt);
   set(hT,'fontsize',14,'color','m')
   set(gca,'fontsize',14)
   title(txtG,'fontweight','normal')
   xlim([0 fMax])

subplot(2,1,2)
   yPmin = min(log10(Ph./max(Ph))); 
   xP = [fN fN]; yP = [yPmin 0.4];
   plot(xP,yP,'r','lineWidth',1);
   hold on
   xP = [fD fD]; yP = [yPmin 0.4];
   plot(xP,yP,'m','lineWidth',1);

   xP = f; yP = log10(Ph./max(Ph));    
   plot(xP,yP,'b','lineWidth',2);
    
   ylabel('log( PSD )');
   xlabel('f  [ Hz ]');
   grid on 
   
   set(gca,'fontsize',14)
   xlim([0 fMax])

end


 
    toc

%  FUNCTIONS  =========================================================

function sDot = EM(t,s,K)

  X = s(1);
  w = s(2);
  sDot(1) = w;
  sDot(2) = -K(1)*sin(X) - K(2)*w + K(3)*cos(K(4)*t);

  sDot = sDot';


end

