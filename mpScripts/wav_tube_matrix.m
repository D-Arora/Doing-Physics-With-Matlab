% wav_tube_matrix.m
% Eigenvalues and eigenfunctions for standing waves in tubes 
% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% 20140718


clear all
close all
clc


% INPUTS =================================================================
N = 2000; num = N;      % number of data points

soln = 4;
Neig = 8;              % no. of eigenvalues to be found  must be EVEN
% Setup ==================================================================
cS = 343;               % speed of sound in air

% Specify type of instrument ---------------------------------------------
disp('  ');
disp('     Select type of instrument');
disp('  ');
disp('          1:   Cylindical pipe:     Closed   Open');
disp('   ');
disp('          2:   Cylindrical pipe:    Open    Open');
disp('   ');
disp('          3:   Conical pipe (cone): Open     Open');
disp('   ');
disp('          4:   Conical pipe (cone):  Closed   Open');
disp('   ');
disp('          5:   Oboe:                 Open   Open');
disp('   ');
disp('          6:   Trumpet:              Closed   Open');
disp('   ');
disp('          8:   Human vocal tract "ah ...": two pipe model');
disp('   ');
disp('          9:   Human vocal tract "ah ..."');
disp('   ');

flag_inst = input('     Enter 1, 2, ... , 9:   ');
disp('   ');

switch flag_inst

case 1   % Uniform organ pipe: Closed Open  111111111111111111111111111111
   disp('   ');
   L = input('   Length of cylindical pipe (m),  L  =  ');
   disp('   ');
   R = input('   Radius of cylindical pipe (m), R =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end == 1
       x_max = L + 0.6 * R;
   else
       x_max = L;
   end
   x_min = 0;
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);  dx2 = dx^2;   
   A = (pi*R^2).* ones(1,num);    
   
    p(1) = 1; U(1) = 0;  flagBC = 1;
   
 case 2   % Uniform organ pipe: Open Open 22222222222222222222222222222222
   disp('   ');
   L = input('   Length of cylindical (m),  L  =  ');
   disp('   ');
   R = input('   Radius of cylindical pipe (m), R =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end == 1
       x_max = L + 0.6 * R;
       x_min = -0.6 * R;
   else
       x_min = 0; x_max = L;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1); dx2 = dx^2;    
   A = (pi*R^2).* ones(1,num);    
   
   p(1) = 0; U(1) = 1/A(1);  flagBC = 0;
   
   
 case 3   % Conical pipe (cone): Open Open 3333333333333333333333333333
   disp('   ');
   L = input('   Length of conical pipe (m),  L  =  ');
   disp('   ');
   R_min = input('   Min radius of conical pipe (m), R_min =  ');
   disp('   ');
   R_max = input('   Max radius of conical pipe (m), R_max =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R_max; x_min = -0.6 * R_min;
   else
       x_max = L; x_min = 0;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1); dx2 = dx^2;
   R = ((R_max-R_min)/L) .*x + R_min;
   A = (pi .* R.^2);   
   
   p(1) = 0; V(1) = 1/A(1); flagBC = 0;
   
case 4   % Conical pipe (cone): Closed Open 444444444444444444444444444
   disp('   ');
   L = input('   Length of conical pipe (m),  L  =  ');
   disp('   ');
   R_min = input('   Min radius of conical pipe (m), R_min =  ');
   disp('   ');
   R_max = input('   Max radius of conical pipe (m), R_max =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end == 1
       x_max = L + 0.6 * R_max; 
   else
       x_max = L; x_min = 0;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1); dx2 = dx^2;
   R = ((R_max-R_min)/L) .*x + R_min;
   A = (pi .* R.^2);   
   
   p(1) = 1; V(1) = 0; flagBC = 1;
 
case 5   % Oboe 55555555555555555555555555555555555555555555555555555555
    load x_oboe;
    load y_oboe;
    x = x_oboe;
    R = y_oboe;
    
    A = pi .* R.^2;
    A = A';
    
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end == 1
       x = x .* (max(x) + 0.6 * max(R))/max(x); 
   end
    
   dx = x(2)-x(1); dx2 = dx^2;
   x_min = min(x);
   x_max = max(x);
    
   p(1) = 0; V(1) = 0; flagBC = 1; 
   
case 6   % Trumpet 666666666666666666666666666666666666666666666666666
    load x_trumpet;
    load y_trumpet;
    x = x_trumpet;
    R = y_trumpet;
        
    A = pi .* R.^2; A = A';
    
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x = x .* (max(x) + 0.6 * max(R))/max(x); 
   end
    
   dx = x(2)-x(1); dx2 = dx^2;
   x_min = min(x);
   x_max = max(x);
    
   p(1) = 1; V(1) = 0;  flagBC = 1; 
   
   % scaled trumpet
   x = 0.6253 .* x ./ x_max;
   R = 1.82e-2 .* R ./ R(end);
   
    p(1) = 0; V(1) = 0; flagBC = 1; 
   
 case 7  % Organ pipe - hole:  open open 7777777777777777777777777777

     disp('   ');
     disp('   length of organ pipe L  = 0.800 m ');
     disp('   ');
     disp('   radius of organ pipe R  = 0.012 m ');
     disp('   ');
     disp('   hole located 0.160 m from end of organ pipe');
     disp('   ');
     disp('   hole radius  R_hole  = 0.006 m');
     disp('   ');
     
     L = 0.800;  R =  0.012;  
     
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R;
   else
       x_min = 0; x_max = L;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);     
   A = (pi*R^2).* ones(1,num);    
   
   p(1) = 0; V(1) = 1/A(1);   
   
   
 case 8   % Human vocal tract "ah ...": two pipe model 8888888888888888
     disp('   ');
     disp('   length of pipe 1, L1  = 0.090 m ');
     disp('   ');
     disp('   radius of pipe 1, R1  = 5.642e-3 m ');
     disp('   ');
     disp('   length of pipe 2, L2 = 0.080 m');
     disp('   ');
     disp('   radius of pipe 2, R2 = 1.4927e-2 m');
     disp('   ');
     
   L1 = 0.09; R1 = 5.642e-3;
   L2 = 0.08; R2 = 1.4927e-2;  
   L = L1+L2; R = R1;
   
   
   x_min = 0; x_max = L;
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1); dx2 = dx^2;    
   A = (pi*R^2).* ones(1,num);    
   for c = 1 : num                 % can rewrite using logical functions
   if x(c) > L1
       R = R2;
       A(c) = (pi*R^2);  
   end
   end
   
   p(1) = 1; V(1) = 0;   flagBC = 1;
   
 case 9   % Human voice tract: "ah ..."  999999999999999999999999999999
    load x_voice;
    load y_voice;
    x = x_voice;
    R = y_voice;
    A = pi .* R.^2; A = A';
    
   disp('   ');
   x_min = 0; x_max = max(x);
    
   dx = x(2)-x(1); dx2 = dx^2;
   x_min = min(x);
   x_max = max(x);
    
   p(1) = 1; V(1) = 0;  flagBC = 1;
end

% MATRICES ===============================================================

% Area Matrix
dAdx = gradient(A,dx);                 % gradient of A wrt x
dAdx_A = dAdx ./A;
AREA_matrix = zeros(N-2,N-2);
for cc = 1 : N-2;
AREA_matrix(cc,cc) = dAdx_A(cc+1);
end

% Pressure gradient matrix
offP = ones(N-3,1); 
dP_matrix = (zeros(N-2) + diag(offP,1) - diag(offP,-1));
if flagBC == 1; dP_matrix(1,1) = -1; end;
dP_matrix = dP_matrix ./ (2*dx);

AP_matrix = AREA_matrix * dP_matrix;

% Make Second Derivative Matrix
off = ones(N-3,1);                 
SD_matrix = (2*eye(N-2) - diag(off,1) - diag(off,-1));
if flagBC == 1; SD_matrix(1,1) = 1; end;
SD_matrix = SD_matrix ./dx2;

Matrix = SD_matrix - AP_matrix ; 
 

% FIND EIGEN - VALUES & FUNCTIONS ========================================

SIGMA = 'sm';  
 [e_funct e_values] = eigs(Matrix,Neig,SIGMA);

% All Eigenvalues 1, 2 , ... N-2    wave number k
k2 = zeros(1,Neig);
for cc = 1 : Neig
 k2(cc) = e_values(cc,cc);
end

% Corresponding Eigenfunctions 1, 2, ... ,n
clear y
p = zeros(N,Neig);               % axial pressure
for cn = 1 : Neig
if flagBC == 1; yBC = e_funct(1,cn); else yBC = 0; end;    
p(:,cn) = [yBC; e_funct(:,cn); 0];
end % for

wL = (2*pi)./(sqrt(k2))
f = (cS ./ wL)
fN = f./ min(f)
flag_order = 0;
if fN(1) > fN(end), flag_order = 1; end;

figure(1)
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.1 0.1 0.5 0.55]);
subplot(5,1,1)
plot(x,R,'m','linewidth',2)
hold on
plot(x,-R,'m','linewidth',2)
if flagBC == 1
    plot([x(1) x(1)], [-R(1) R(1)],'m','linewidth',2); 
end
grid on
ylabel('width  (m)','fontsize',12);
xlabel('axial position  (m)','fontsize',12);
set(gca,'YLim',[-1.1*max(R),1.1*max(R)]);

for cc = 1 : 4
   if flag_order == 0, co = cc; f0 = f(1);
      else co = Neig + 1 - cc; f0 = f(end);
   end  
   yp = p(:,co)./ max(abs(p(:,co)));
   subplot(5,1,cc+1);
   plot(x,yp,'linewidth',2)
   hold on
   plot(x,-yp,'linewidth',2)
   set(gca,'Ylim',[-1.1 1.1]);
   ylabel('p  (a.u.)','fontsize',12);
   tx1 = '\lambda  =  ';
   tx2 = num2str(wL(co),'%3.2f');
   tx3 = '  m       ';
   tx4 = 'f  =  ';
   tx5 = num2str(f(co),'%5.1f');
   tx6 = '  Hz       ';
   tx7 = 'f_N / f_1  =  ';
   tx8 = num2str(f(co) / f0,'%2.2f');
   tx = [tx1 tx2 tx3 tx4 tx5 tx6 tx7 tx8];
   xlabel(tx);
   
end

figure(2)
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.3 0.1 0.5 0.55]);
subplot(5,1,1)
plot(x,R,'m','linewidth',2)
hold on
plot(x,-R,'m','linewidth',2)
if flagBC == 1
    plot([x(1) x(1)], [-R(1) R(1)],'m','linewidth',2); 
end
grid on
ylabel('width  (m)','fontsize',12);
xlabel('axial position  (m)','fontsize',12);
set(gca,'YLim',[-1.1*max(R),1.1*max(R)]);

for cc = 5 : 8
   if flag_order == 0, co = cc; f0 = f(1);
      else co = Neig + 1 - cc; f0 = f(end);
   end  
   yp = p(:,co) ./ max(abs(p(:,co)));
   subplot(5,1,cc-3);
   plot(x,yp,'linewidth',2)
   hold on
   plot(x,-yp,'linewidth',2)
   set(gca,'Ylim',[-1.1 1.1]);
   ylabel('p  (a.u.)','fontsize',12);
   tx1 = '\lambda  =  ';
   tx2 = num2str(wL(co),'%3.2f');
   tx3 = '  m       ';
   tx4 = 'f  =  ';
   tx5 = num2str(f(co),'%5.1f');
   tx6 = '  Hz       ';
   tx7 = 'f_N / f_1  =  ';
   tx8 = num2str(f(co) / f0,'%2.2f');
   tx = [tx1 tx2 tx3 tx4 tx5 tx6 tx7 tx8];
   xlabel(tx);
   
end

