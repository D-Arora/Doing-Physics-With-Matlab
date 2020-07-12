% CR001.m

% Modelling a simple voltage divider circuit

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% 171220 


clear all
close all
clc

% INPUTS   ===============================================================
   emf = 24;
   R = [150 300 50 250 100];

% CALCULATIONS   ======================================================== 
  Vm = [0;0;-emf];
  Rm(1,1) = R(1) + R(3)+ R(5);
  Rm(1,2) = R(5);
  Rm(1,3) = R(1);
  
  Rm(2,1) = R(5);
  Rm(2,2) = R(2) + R(4) + R(5);
  Rm(2,3) = -R(2);
  
  Rm(3,1) =  -R(1);
  Rm(3,2) =  R(2);
  Rm(3,3) =  -(R(1) + R(2));
  
  
  Im = Rm \ Vm;
  
  
  IR(1) = Im(1) - Im(2);
  IR(2) = Im(1) - Im(3);
  IR(3) = Im(2);
  IR(4) = Im(3);
  IR(5) = Im(2) - Im(3);
  
  V = IR .* R;
  P = IR.^2 .* R;
  
% DISPLAY RESULTS = =====================================================
   fprintf('emf =  %3.2f  V \n\n',emf)
   disp(' R [ohms]   IR [mA]    V [V]      P [mW]')
   for c = 1 : 5
    fprintf('  %3.0f       %3.3f      %3.3f      %3.2f\n',R(c),1e3*IR(c),V(c),1e3*P(c));   
   end
   
  Im
  
  A =[300,100,150;100,650,-300;-150,300,-450]
  B = [0;0;-24]
  C =1000.*( A\B)
  
  A =[450,-150,-300;-150,300,-100;-300,-100,650]
  B = [0;0;24]
  C =1000.*( A\B)
%    for c = 1 : 4
%    T1 = 'R';
%    T2 = num2str(c);
%    T = strcat(T1,T2);
%    fprintf(T);
%    fprintf(' = %3.0f ohms \n',R(c));
%    
%    T1 = 'IR';
%    T2 = num2str(c);
%    T = strcat(T1,T2);
%    fprintf(T);
%    fprintf(' = %3.2f mA \n',1e3*IR(c));
%    
%    T1 = 'VR';
%    T2 = num2str(c);
%    T = strcat(T1,T2);
%    fprintf(T);
%    fprintf(' = %3.2f V \n',V(c));
%    
%    T1 = 'PR';
%    T2 = num2str(c);
%    T = strcat(T1,T2);
%    fprintf(T);
%    fprintf(' = %3.2f V \n',P(c));
%    
%    disp('  ')
%    end
%    fprintf('R2 = %3.2f ohms \n',R(2))
%    fprintf('R3 = %3.2f ohms \n',R(3))
%    fprintf('R4 = %3.2f ohms \n',R(4))
   
   