% svr2D.m


clear 
close all
clc

% Number of rows and columns
  num = 10;
% Number of time steps
  nT = 5;
  
% Population
  N = num*num;
% State of an individual
  Q = ones(num,num);
% Initial state of infected individuals
  Q(5,5) = -5;
  Q(6,7) = -5;
% Weights between neighbours   5x5 matirx
  wC = 1.4;
  w01 = 1; w11 = 0.7; w02 = 0.5; w12 = 0.4; w22 = 0.3;
  
  W = [w22 w12 w02 w12 w22; w12 w11 w01 w11 w12; w02 w01 0 w01 w02; ...
        w12 w11 w01 w11 w12; w22 w12 w02 w12 w22];

for t = 1:nT    
% Find injected individuals
  [x, y] = find(Q < 0);
  numV = length(x);
  
% Virus matrices
  V = zeros(num,num); VS = zeros(num,num); 
  for c = 1: numV
       
%       if Rx < 1; Rx = 1; end
%       if Rx > num; Rx = num1; end    
      if x(c)-2 > 0 && x(c) < num-1 && y(c)-2 > 0 && y(c) < num-1 
        Rx = x(c)-2:x(c)+2;
        Ry = y(c)-2:y(c)+2;
        V(Rx,Ry) = W; 
        VS = VS + V; 
      end
  end
      
      VS(x,y) = 0;
  
  [x, y] = find(Q == 1);
  numV = length(x);clc
 
  
  for cx = 1:num
   for cy = 1:num
      if VS(cx,cy) > 1.4
          Q(cx,cy) = -5;
      end    
   end
  end
  
end  
  
  
Q

  
  
  

  