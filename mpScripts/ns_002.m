% ns_002.m

clear all
close all
clc
tic

% ========================================================================
% INPUTS 
% =================================================================
% No. time steps / time steps [s] / time interval
  timeStepS = 0.001;                  % 1 msec
  spikesPerS = 50;                    % 50 spikes per second, on average
  durationS = 1000;                  % 1 sec simulation
  t = 0:timeStepS:durationS; 
  Nt = length(t);
  % dt = 0.001;
  % t1 = 0;
  % t2 = 1;
  % t = 0:dt:t2;
% No . spike trains / threshold for firing of neutron
   Nst = 1;
   vT = spikesPerS * timeStepS;

% ========================================================================
% SETUP
% ========================================================================
% times
  % t = linspace(t1,t2,Nt);

%  random numbers (0 to 1) / spike trains
   Nrand = rand(Nst,Nt);
   spikes = vT > Nrand;
   
% time for spikes
   tInc = t(2) - t(1);
   nSpike = find(spikes);
   tSpikes = nSpike .* tInc; L = length(tSpikes);
   spikeIntervals = tSpikes(2:L) - tSpikes(1:L-1);

   tInc = durationS/50;
   for cc = 1: 50 
      cum_counts(cc) = size(tSpikes(tSpikes < tInc*cc),2) ;
   end
   counts = cum_counts(2:cc) - cum_counts(1:cc-1);
% ========================================================================   
% GRAPHICS
% spike trains gap / height
dy = 1;
yH = 5;
xP = [t(1) t(1)]; yP = [dy yH];

% figure(1)   %-------------------------------------------------------------
% set(gcf,'units','normalized','position',[0.1 0.1 0.42 0.8]);
% hold on
% for cst = 1 : Nst
%    for ct = 1 : Nt
%       if spikes(cst,ct) == 1
%          xP = [t(ct) t(ct)];
%          plot(xP,yP,'b','linewidth',1);
%       end
%    end
%    yP(1) = yP(2) + dy; yP(2)= yP(1) + yH;   
% end
% axis off
% text(0.35, -2,'time t   [0  to 1.00 s]','fontsize',14);
% title('Spike Trains','fontsize',14);

figure(2);   % ------------------------------------------------------------                                              
%x = [1:binSize:100];
%intervalDist = hist(spikeIntervals(spikeIntervals < 100), x);
%intervalDist = intervalDist / sum(intervalDist) / binSize; % normalize by dividing by spike number
%bar(x, intervalDist);
hist(spikeIntervals,100);
hold on


figure(99)
hist(counts,50)

toc