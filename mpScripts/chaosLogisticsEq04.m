% chaosLogisticsEq04.m

% LOGISTIC DIFFERENCE EQUATION x(n+1) = 4*r*x*(1-x)
% Finding the fixed points
% Inputs: growth factor r (0 <r < 1)
%         dynamics period P (P = 1, 2, ... , 8)  

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220727 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%  https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%  https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/chaosLogisticEq.htm
% Scripts
%  https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%  https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

clear
close all
clc
tic

syms x 

% Input growth factor r   >>>>>
  r = 0.962;

% Input dynamics period  P = 1, 2, 3, ... , 8 >>>>>
  P = 6;


% ITERATIONS =========================================================
% [1D] mapping functions
 f1 = 4*r*x*(1-x);
 f2 = 4*r*f1*(1-f1);
 f3 = 4*r*f2*(1-f2);
 f4 = 4*r*f3*(1-f3);
 f5 = 4*r*f4*(1-f4);
 f6 = 4*r*f5*(1-f5);
 f7 = 4*r*f6*(1-f6);
 f8 = 4*r*f7*(1-f7);
 
 if P == 1; f = f1; end
 if P == 2; f = f2; end
 if P == 3; f = f3; end
 if P == 4; f = f4; end
 if P == 5; f = f5; end
 if P == 6; f = f6; end
 if P == 7; f = f7; end
 if P == 8; f = f8; end

% Solve f(xe) = xe
  eq1 = f - x;
% Fixed points  
  xe = vpasolve(eq1);
% Gradient df/dx  
  df = gradient(f,[x]);
% Symbolic to numeric
  f_dash = subs(df(1),x,{xe});
% Print results to Command window
  output(r, xe, f_dash, P)

toc


% FUNCTIONS ========================================================

function output(r, xe, f_dash, P)
   
   fprintf('r = %2.6f   \n',r)
   fprintf('Dynamic period P = %2.0f  \n',P)
   disp('  ')
    for c = 1: length(xe)
       fprintf('xe = %2.6f   f_dash = %2.3f  \n',xe(c), f_dash(c))
    end
              
    disp('  ')
    flagS = 0;
    for c = 1: length(xe)
      if abs(f_dash(c)) <= 1
      if flagS == 0; disp('Stable fixed points'); end
      fprintf('   xe = %2.6f   f_dash = %2.3f  \n',xe(c), f_dash(c))
      flagS = 1;
      end
    
    end
    disp('  ')  
end