%CN03.m

% RC FILTER CIRCUIT
% Finite Difference Method for RC circuit anlysis

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180121

clear all
close all
clc

% INPUTS: S.I. UNITS ====================================================
  R  = 1.0e3;    % [1e3]
  RO = 1.0e6;    % [1e6]
  C  = 1e-6;     % [1e-6]
   
  VS = 10;       % [10]
  

% =======================================================================  
% INPUTS: select voltage source
   flagV = 1; 
%  ======================================================================   

% Time constant and time step
   tau = R*C;
   dt = tau/100;
   
   switch flagV
    case 1  % step function  off to on
       % step height  [10]
         VS = 10;              
       % max time interval  [10e-3] 
         tMax = 10e-3;        
       % percentage off time  [10]
         pOFF = 10;         
                
         t = 0:dt:tMax;
         N = length(t);
         nOFF = round(N*pOFF/100);
         vS = zeros(1,N);
         vS(nOFF:end) = VS;
     
    case 2  % step function  on to off
       % step height  [10]
         VS = 10;              
       % max time interval  [10e-3] 
         tMax = 10e-3;        
       % percentage on time  [50]
         pON = 50;         
       
         t = 0:dt:tMax;
         N = length(t);
         nON = round(N*pON/100);
         vS = zeros(1,N);
         vS(1:nON) = VS;
         
    case 3    % pulses ON / OFF
      % pulse height  [10]
         VS = 10;              
       % max time interval  [10e-3] 
         tMax = 10e-3;        
       % percentage off time  [10]
         pOFF = 50;         
       % Number of periods     [3]
         nT = 3;
       
       % Calculations  
         t = 0:dt:tMax;
         N = length(t);  
                  numT = round(N/nT-0.5);
         nON = round(numT*pOFF/100);
         nOFF = numT-nON;
         vS = zeros(1,N);
         for c = 1 : nT
           vS(nOFF + numT*(c-1) : nOFF + numT*(c-1)+nON) = VS;
         end

    case 4     % sinusoidal input
       f = 10;
       w = 2*pi*f;
       T = 1/f;
       nT = 5; 
       tMax = nT*T;
       t = 0:dt:tMax;
       N = length(t);
       vS = VS .* sin(w*t);
       
       fc = (1/(2*pi))* (R + RO)/(R*RO*C);    % cut-off frequency
       
     case 5     % sinusoidal input - rectified
       dt = tau/500;
       f = 200;
       w = 2*pi*f;
       T = 1/f;
       nT = 8;
       tMax = nT*T;
       t = 0:dt:tMax;
       N = length(t);
       vS = abs(VS .* sin(w*t));

     case 6     % Input: superposition of sinusoidal signals
       f = 100;
       w = 2*pi*f;
       T = 1/f;
       nT = 3;
       tMax = nT*T;
       t = 0:dt:tMax;
       N = length(t);
       vS = VS .* sin(w*t) + (0.25*VS) .* sin(20*w*t);

   end
  
   
% CALCULATIONS ==========================================================  
  
%  Initialise values 
  vR = zeros(1,N);
  vC = zeros(1,N);
  iR = zeros(1,N);
  iRO = zeros(1,N);
  iC = zeros(1,N);
  qC = zeros(1,N);
  
%  Time Step #1
   vR(1)  = vS(1);
   vC(1)  = 0;            % assume capacitor initially uncharged
   iR(1)  = vR(1)/R;
   iC(1)  = iR(1);
   iRO(1) = 0;

   k = dt/C;
% Time Steps #2 to #N
  for c = 2 : N
    vC(c)  = vC(c-1) + k*iC(c-1);
    vR(c)  = vS(c) - vC(c);
    iR(c)  = vR(c)/R;
    iRO(c) = vC(c)/RO;
    iC(c)  = iR(c) - iRO(c);
    qC(c)   = qC(c-1) + iC(c-1)*dt;
  end
 
% Peak values 
  VS  = max(vS);
  VR  = max(vR);
  IR  = max(iR);
  VC  = max(vC);
  IC  = max(iC);
  IRO =  max(iRO);

  % Powers and energy
     iS = iR;
     pS = vS .* iS;
     pO = vC .* iRO;
     pC = vC .* iC;
     uS = zeros(1,N); uO = zeros(1,N); uC = zeros(1,N);
     for c = 2 : N
         uS(c) = uS(c-1) + pS(c)*dt;
         uO(c) = uO(c-1) + pO(c)*dt;
         uC(c) = uC(c-1) + pC(c)*dt;
     end
     

% =====================================================================
%   DISPLAY RESULTS
% =====================================================================
     disp('Circuit Elements  ') 
     fprintf(' Resistance  R =  %3.2f   ohms \n',R);
     fprintf(' Rutput resistance  Ro =  %3.2f  ohms  \n',RO);
     fprintf(' Capacitance  C  = %3.2e  F  \n',C);
     disp('  ');
     fprintf(' Time constant  tau  =  %3.2e s  \n',tau);
     fprintf(' time step      dt   =  %3.2e s  \n',dt);
     fprintf('          tau / dt   =  %3.2f    \n',tau/dt);
     if flagV == 4
        fprintf('  fS   =  %3.2f    \n',f);
        fprintf('  fC   =  %3.2f    \n',fc);
        fprintf('  Av   =  %3.2f    \n',max(vC)/max(vS));
     end
     
     
     % ======================================================================
% GRAPHICS
% ======================================================================

figure(1)    
   FS = 14;
   pos = [0.07 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
   subplot(3,1,1)
     xP = 1e3.*t; yP = vS;
     plot(xP,yP,'b','linewidth',2');
     grid on
     %xlabel('t  [ ms ]');
     ylabel('v_S  [ V ]');
     
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.1e');
     tm3 = ' \Omega   R_o = ';
     tm4 = num2str(RO, '%3.1e');
     tm5 = ' \Omega   C = ';
     tm6 = num2str(C, '%3.1e');
     tm7 = ' F';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     
         
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]);
     set(gca,'yLim',[-VS-0.2 VS+0.2]);
     title(tm,'fontweight','normal','fontsize',12);
     set(gca,'yTick',-10:5:10);
     
   subplot(3,1,2)
     xP = 1e3.*t; yP = vC;
     plot(xP,yP,'r','linewidth',2');
     grid on
     %xlabel('t  [ ms ]');
     ylabel('v_C =  v_O  [ V ]');
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]); 
     set(gca,'yLim',[-VS-0.2 VS+0.2]);
     set(gca,'yTick',-10:5:10);
     tm1 = ' \tau  =  R C  =  ';
     tm2 = num2str(tau, '%3.2e');
     tm3 = '   s';
     tmx = [tm1 tm2 tm3];
     title(tmx,'fontweight','normal','fontsize',12)
     
   subplot(3,1,3)
     xP = 1e3.*t; yP = vR;
     plot(xP,yP,'m','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('v_R  [ V ]');
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]);
     set(gca,'yLim',[-VS-0.2 VS+0.2]);
     set(gca,'yTick',-10:5:10);
   
figure(2)    
   FS = 14;
   pos = [0.37 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
  subplot(3,1,1)
     xP = 1e3.*t; yP = pS;
     plot(xP,yP,'b','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('p_S  [ W ]');
      
     tm1 = 'R = ';
     tm2 = num2str(R, '%3.1e');
     tm3 = ' \Omega   R_o = ';
     tm4 = num2str(RO, '%3.1e');
     tm5 = ' \Omega   C = ';
     tm6 = num2str(C, '%3.1e');
     tm7 = ' F';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]);
     title(tm,'fontweight','normal','fontsize',12);
          
     subplot(3,1,2)
     xP = 1e3.*t; yP = pO;
     plot(xP,yP,'r','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('p_O  [ W ]');
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]); 
     
     subplot(3,1,3)
     xP = 1e3.*t; yP = pC;
     plot(xP,yP,'m','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('p_C  [ W ]');
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)])
  

  figure(3)    
   FS = 12;
   pos = [0.67 0.05 0.28 0.72];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   
  subplot(3,1,1)
     xP = 1e3.*t; yP = uS;
     plot(xP,yP,'b','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('u_S  [ J ]');
      
     tm1 = '      R = ';
     tm2 = num2str(R, '%3.1e');
     tm3 = ' \Omega   R_o = ';
     tm4 = num2str(RO, '%3.1e');
     tm5 = ' \Omega   C = ';
     tm6 = num2str(C, '%3.1e');
     tm7 = ' F ';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]);
     title(tm,'fontweight','normal','fontsize',12);
     
     subplot(3,1,2)
     xP = 1e3.*t; yP = uO;
     plot(xP,yP,'r','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('u_O  [ J ]');
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]); 
     
     subplot(3,1,3)
     xP = 1e3.*t; yP = uC;
     plot(xP,yP,'m','linewidth',2');
     grid on
     xlabel('t  [ ms ]');
     ylabel('u_C  [ J ]');
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)])

  
 figure(4)
   FS = 14;
   pos = [0.50 0.05 0.28 0.32];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
      
   yyaxis left 
   xP = 1e3.*t; yP = 1e3.*iC;
   plot(xP,yP,'b','linewidth',2');
   grid on
   xlabel('time  t  [ms]');
   ylabel('i_C  [ mA ]');
   tm1 = '      R = ';
     tm2 = num2str(R, '%3.1e');
     tm3 = ' \Omega   R_o = ';
     tm4 = num2str(RO, '%3.1e');
     tm5 = ' \Omega   C = ';
     tm6 = num2str(C, '%3.1e');
     tm7 = ' F ';
     tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7];
   
     ax = gca;
     ax.YAxis(1).Color = 'b';
     set(gca,'fontsize',FS);
     set(gca,'xLim',[0 1e3*max(t)]);
     title(tm,'fontweight','normal','fontsize',12);
     yyaxis right
     xP = 1e3.*t; yP = qC;
     plot(xP,yP,'m','linewidth',2');
     ylabel('q_C  [ C ]');
     ax.YAxis(2).Color = 'm';
   
  
   
  


