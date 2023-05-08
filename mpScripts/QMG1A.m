% QMG1A.m

clear; close all; clc

% PROBABILITY: PEOPLE IN A ROOM
% j   age of people from 10 to 22 
  j = 10:22; j = j';
% L   number of ages  
  L = length(j);
% Nj   number of people aged j: distribution function
  Nj = round( 50.*exp(-(18-j).^2./12) + 3.*randn(L,1) );
  Nj(Nj < 0) = 0;
% N total number of people
  N = sum(Nj);
% P   Probability of person of age j
  Pj = Nj/ N;
% P  Probability of all ages
  P = sum(Pj);
% Pmax   Most probable age
  ind = find(Pj == max(Pj));
  jMax = j(ind);
% avgAge   average (mean) age
  jAvg = sum(j.*Pj) ; 
% medianAge   median age
  cumS = cumsum(Nj);
  ind = find(cumS > N/2,1);
  jMedian = j(ind)-1;
% avgSq   average of the squares
  j2Avg = sum(j.^2.*Pj);
% variance and standard deviation (sigma)
  variance =  j2Avg - jAvg^2;
  sigma = sqrt(variance);

% DISPLAY
  fprintf('\n Total number of people, N = %2.0f  \n',N);
  table(j, Nj,Pj)
  fprintf('\n Average age, jAvg = %2.0f  \n',jAvg);
  fprintf('\n Most probable age, jMax = %2.0f  \n',jMax);
  fprintf('\n Median age, jAvg = %2.0f  \n',jMedian);
  fprintf('\n Average of the squares of the ages, j2Avg = %2.3f  \n',j2Avg);
  fprintf('\n Variance = %2.3f  \n',variance);
  fprintf('\n standard deviation, sigma = %2.3f  \n',sigma);
 
% GRAPHICS    
   figure(1)
   pos = [0.02 0.05 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   bar(j',Nj')
   xlabel('age  j')
   ylabel('N_j')
   set(gca,'fontsize',14)
