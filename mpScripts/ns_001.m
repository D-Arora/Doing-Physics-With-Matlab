function ns_001
  % ns_001.m

  % 04 jan 2016
  % Solving the Time Dependent Schrodinger Equation using the FDTD Method
  % Motion of a wave packet
  % Ian Cooper
  % School of Physics, University of Sydney
  % documentation: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
  % mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts

  
  close all
  clc
  global t;
  setup
end

% ------------------------------------------------------------------------
function setup
  global t;
  timeStepS = 0.001;                  % 1 msec
  spikesPerS = 50;                    % 50 spikes per second, on average
  durationS = 1.000;                  % 1 sec simulation
  t = 0:timeStepS:durationS;	    % a vector with each time step	

  vt = rand(size(t));
  spikes = (spikesPerS*timeStepS) > vt;
  rasterPlot(spikes, timeStepS);
end

% ------------------------------------------------------------------------
function rasterPlot(spikes, timeStepS)
  figure(1);
  axes('position', [0.1, 0.1, 0.8, 0.8]);
  axis([0, length(spikes) - 1, 0, 1]);
  trains = size(spikes, 1); 
  ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
  ticHeight = (1.0 - (trains + 1) * ticMargin) / trains;

  %for train = 1:trains
    spikeTimes = find(spikes == 1)
        
     for cc = 1:length(spikeTimes)
       line([spikeTimes(cc), spikeTimes(cc)], [0.01 0.99]);
     end
  %end
end

