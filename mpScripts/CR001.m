% CR001.m

% Modelling a simple voltage divider circuit

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 171220


clear all
close all
clc

% INPUTS   ===============================================================
   emf = 10;
   R = [0 100 200 400];

% CALCULATIONS   ======================================================== 
  Vm = [emf;0];
  Rm(1,1) = R(1) + R(2 )+ R(3);
  Rm(1,2) = -R(3);
  Rm(2,1) =  R(3);
  Rm(2,2) = -(R(3) + R(4));
  
  Im = Rm \ Vm;
  
  IR(1) = Im(1);
  IR(2) = Im(1);
  IR(3) = Im(1)- Im(2);
  IR(4) = Im(2);
  
  V = IR .* R;
  P = IR.^2 .* R;
  
% DISPLAY RESULTS = =====================================================
   fprintf('emf =  %3.2f  V \n\n',emf)
   disp(' R [ohms]   IR [mA]    V [V]      P [mW]')
   for c = 1 : 4
    fprintf('  %3.0f       %3.2f      %3.3f      %3.2f\n',R(c),1e3*IR(c),V(c),1e3*P(c));   
   end
   
  
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
   
   