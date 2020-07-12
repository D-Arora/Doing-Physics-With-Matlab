
function poissonTutorial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% poissonTutorial.m
%
% NOTE: For compatibility with different web browsers,
% this file has been renamed 'poissonTutorial.txt'. In
% order to use it in MATLAB, you should rename it as
% 'poissonTutorial.m' and then open it up with the 
% MATLAB editor.
%
% This tutorial introduces the Poisson model of stochastic
% neuronal firing. It will demonstrate the type of variability
% seen in the response of sensory neurons, and explore the
% properties of that variability
%
% This tutorial is closely based on a tutorial created by 
% David Heeger for the Cold Spring Harbor course on Computational
% Neuroscience: Vision, which was subsequently extended by 
% Greg Horwitz.
%
% INSTRUCTIONS:
%
% While other approaches are possible, it is recommended that you
% put a breakpoint on the first line of executable code below, and
% then use the Matlab debugger to step through the lines while
% reading the text that is interspersed with code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;          % clear, but don't move the windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTRODUCTION:
%
% It is not important that you memorize the equations for a Poisson
% process, but you should know something about its general properties
% how they relate to neuronal spiking.
%
% A Poisson process is a member of the class of stochastic (random)
% point (binary) processes.  It is formally defined by the following 
% set of equations (here, N(t) is the number of events that have occurred 
% by time 't').
%
% 1) N(0)=0 
%    This means that there are no spikes at the begining of time.
% 2) The process is stationary
%    This means that the rate of spike firing doesn't change over time,
% although neuronal responses can be modeled using Poisson processes with
% different rates during different epochs (e.g. background and stimulus periods).
% 3) The process has independent increments
%    This means that the number of spikes occurring in any interval is independent
% of the number of spikes occurring in any other (non-overlapping)
% interval. That is, what happened in the previous interval doesn't affect
% the current interval
% 3) Prob{N(h)=1} = lamba*h + o(h)
%    Translation: The probability of one event in a "brief" window of
% time is a constant, "lambda" multiplied by the duration of the window, "h".
% "o(h)" is the probability of more than one spike in a brief window, and it
% is vanishingly small. That is, spike counts in brief intervals are all zeros
% and ones.  For neurons, this applies if we use ~1 ms intervals.
%
% The Poisson model occurs very naturally in situations in which a
% large number of independent components all have the capacity 
% for generating an event, but none of them are very likely to do so
% at a given time.  For example, it does an excellent job of describing
% radioactive decay.
%
% Early neurophysiologists were struck by the irregularity spike trains
% recorded from neurons in vivo and thus the Poisson model of neural
% firing was born.  Even these early neurophysiologists were aware
% that the firing rates of real neurons are not stationary and don't have strictly
% independent increments.  Still, these violations aside, the Poisson
% model does a remarkable job of describing firing variability of cortical
% neurons.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATING POISSON SPIKING:
%
% One of the easiest ways to generate a Poisson spike train is to rely
% on the approximation:
%
%    Prob{1 spike during (t, t + timeStepS)} = r(t) * timeStepS
%
% where r(t) is the instantaneous firing rate (in spikes per second) and timeStepS 
% is the time step.  The right half of equation is essentially "lamba * h" from 
% equation 3 above, with the assumption that "o(h)" is zero -- timeStepS is short 
% enough so that there is no chance that the neuron would fire more than one spike 
% in any given interval.

% Let's begin by choosing a time step and by choosing an average
% firing rate.

timeStepS = 0.001;                  % 1 msec
spikesPerS = 50;                    % 50 spikes per second, on average
durationS = 1.000;                  % 1 sec simulation
times = 0:timeStepS:durationS;	% a vector with each time step		

% Now we choose random numbers, one for each time step, unformly distributed between 
% 0 and 1. We will use these to decide whether a spike has occurred at each time step.

vt = rand(size(times));

% Finally, create a vector of ones and zeros depending on whether the probability of firing
% (spikesPerS*timeStepS) is greater than the corresponding random number.  We
% have a function rasterPlot() that will plot the result as raster, using lines
% to show which time steps contain a spike.  This should appear in Figure 1,
% which you can move to a convenient position on your screen.  You might want to
% expand the window to see spikes that closely spaced.

spikes = (spikesPerS*timeStepS) > vt;
rasterPlot(spikes, timeStepS);

% The main point is that spikes occur at highly irregular times.  You might seen
% some signs of a pattern, but there is no reliable pattern in the occurrence of
% spikes.  If we make a new set of spike times, and plot that, you will see that
% it is different.  To make things simpler, we have a function makeSpikes() that
% will create spike trains in the manner described above, taking arguments for
% the time step, spike rate, and sample duration.  The function is at the end of this
% file if you want to see the code.  If you run the next two lines
% of code, a new spike train will appear in Figure 1.  

spikes = makeSpikes(timeStepS, spikesPerS, durationS);
rasterPlot(spikes, timeStepS);

% Because the Poisson process is stochastic, each spike train will look
% different, even though they have the same average number of spikes.  The
% variability is more easily seen if we display a set of spike trains together
% in the plot, in format know as a raster plot (each spike train is one line in
% the raster).

spikes = makeSpikes(timeStepS, spikesPerS, durationS, 25);
rasterPlot(spikes, timeStepS);

% Your eye might pick out patterns in the raster, but there is not reliable
% structure.  Each spike is independent of others in its train, and each spike
% train is independent of all others. When testing conditions are held constant
% in neurophysiological experiments, the spiking of individual neurons often
% looks like this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTER-SPIKE INTERVAL DISTRIBUTIONS:
%
% An interesting property of the Poisson process is that the interspike
% intervals are distributed exponentially.  This is a straightforward
% consequence of the definition of the Poisson process.  If the probability of a
% spike in each interval is p, the probability of seeing an interspike interval of 
% 1 interval is p.  The probability that there will not be a spike in the first
% interval is (1-p), so the probability of seeing an interspike interval of 2 is
% (1-p) * p.  The probability of seeing an interspike interval of 3 is (1-p)^2 *
% p.  Because the probability of going an extra interval without a spike
% decreases the probability by a factor of (1-p), the probability distribution
% for interspike interval will decrease by a fixed fraction for successive
% intervals, which will be an exponentially decaying function.  

% Let's plot the interspike intervals for our Poisson process. First, get the 
% spike times.

durationS = 1000;                                       % get lots of intervals
spikes = makeSpikes(timeStepS, spikesPerS, durationS);  % make a spike train
spikeTimes = find(spikes) * timeStepS * 1000;           % get times when spikes occurred (ms)

% We can get the inter-spike intervals by subtraction the times of successive
% spikes.

spikeIntervals = spikeTimes(2:length(spikeTimes)) - spikeTimes(1:length(spikeTimes) - 1);

% Compute a histogram of the spike intervals, normalized to unit volume.

figure(2);                                              % use Figure 2
binSize = 1;                                            % 1 ms bins
x = [1:binSize:100];
intervalDist = hist(spikeIntervals(spikeIntervals < 100), x);
intervalDist = intervalDist / sum(intervalDist) / binSize; % normalize by dividing by spike number
bar(x, intervalDist);

% As promised, the histogram of interspike intervals looks like an exponential 
% probability distribution.  For verificiation, we'll superimpose the exponential
% probability density function.
 
y = exppdf(x, 1 / (spikesPerS * timeStepS));            % exponential function, scaled to predicted max
axis([min(x) max(x) 0 max(y) * 1.1]);                   % make sure everything shows on the plot
xlabel('Interspike interval');
ylabel('Probability');
hold on;
plot(x, y, 'r');                                        % add exponential function 
hold off;

% The exponential interspike interval distribution for a Poisson process
% points to one way that real neurons will not be perfectly Poisson. Real neurons will
% have interspike intervals that depart from a perfectly exponential
% distribution.
%
% QUESTION 1: Why will real neurons not have a perfectly exponential
% interspike interval distribution?
%
% Before leaving the interval distribution, it is worth mentioning that although
% a Poisson process will have an exponential interval distribution, the
% existence of an exponential interval distribution does NOT necessarily mean that it arose
% from a Poisson process.  You could get an exponential distribution in a
% situation where the interval had some special dependence on each other, 
% which is not allowed for a Poisson process. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE COUNT VARIANCE:
%
% The Poisson model of spiking is important because it makes it possible
% to explore the variability of neuronal spiking and how that affect the
% reliability of neuronal signals. One notable property of a Poisson process 
% is that for a given spike rate mean of the spike counts over repeated counts
% is equal to the variance of those spike counts. By "variance", we mean the
% formal statistical definition of variance, which is the square of the standard
% deviation of a distribution.  NB: The equality of the mean and variance holds
% only for spike COUNTS, and is generally NOT the case for the spike RATES or
% interspike intervals of Poisson processes.
%
% This relationship between spike count mean and spike count variance can be
% seen by plotting some distributins of spike counts. Let's see what spike
% counts are for 1 second counting intervals when the average rate is 25
% spikes/s.

timeStepS = 0.001;                  % set the spike train parameters
spikesPerS = 20;                    % 20 spike/s rate
durationS = 1.000;                  % 1 s simulation
startS = 0.001;                     % count from for the full 1 second
endS = 1.000;
spikes = makeSpikes(timeStepS, spikesPerS, durationS, 250);  % make 250 spike trains
counts = countSpikes(spikes, timeStepS);                     % get the 250 spike counts
plotSpikeCounts(counts);            % plot the 250 spike counts (this function is below)

% The counts are random, but you should see that they range from as low as ~10
% to as high as ~30 -- a broad range around the expected 20 spikes.  The plot in
% the lower half of the figure has a point showing the spike count mean for this distribution
% (x-axis) and the spike variance (y-axis).  The point should fall close to the
% identity line.
%
% Let's add some more distributions with mean counts that are lower (2 and ~6).

for i = 1:2
    spikesPerS = spikesPerS / sqrt(10);
    spikes = makeSpikes(timeStepS, spikesPerS, durationS, 250);
    counts(i+1,:) = countSpikes(spikes, timeStepS, startS, endS);
end
plotSpikeCounts(counts);

% The lower plot shows that the equality between variance and mean holds
% regardless of what the mean count is.  
%
% The plots of the three distributions show how the Poission distribution differs from the 
% normal (also known as bell-shaped or Gaussian) distribution.  When a Poisson distribution has 
% a large mean, as with the blue distribution (mean 20 spikes), it is similar to a normal
% distribution.  For small means, however, the Poisson distribution differs
% greatly.  The red distribution (mean 2) is not clipped: Poisson distributions
% do not include negative numbers (there are no negative counts).  When the mean
% counts are small, the bins for 0 and 1 contain most of the counts.
%
% Variance is a somewhat abstract quantity, but you should be able to convince
% yourself that the means are approximately equal to the variances by
% considering the standard deviations, for which you should have some intuition.
% The standard deviation is the square root of the variance, so Poisson distributions
% with means of 2, 6 and 20 should have standard deviations of about 1.5, 2.5
% and 4.5.  You should see that the red, green and blue distributions have
% roughly these standard deviations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQUIVALENCE OF TIME, NEURONS AND RATE FOR SPIKE COUNTS:
%
% The last point we will treat here is how Poisson processes produce an equivalence between
% time, neurons and rate.  Because the occurrence of every spike is independent
% of every other spike, adding two Poisson distributions will give a
% distribution that is also Poisson.  Thus, if a neuron is firing an average of 10 spikes a second
% and you repeatedly count spikes in pairs of 1 second periods, the distribution of (spikes from period 1) +
% (spikes from period 2) will be a Poisson distribution with a mean count of 20
% spikes. If you instead had two neurons that were both firing an average of 10 
% spikes a second and you repeated counted spikes from both over 1 second, their combined counts 
% would also yield a Poisson distribution with a mean of 20.  Both distributions
% wuold be the same as that obtained by counting for 1 second periods from a
% neuron that was firing an average of 20 spikes per second.  It doesn't matter
% if doubled the mean by counting spike for twice as long, counting from
% twice as many (identical) neurons, or from a neuron firing at twice the rate,
% the result is the same.
%
% Let's look at this directly.  The following code will plot 4 distributions.
% Each includes 250 counts. The first (blue) distribution is spikes counts for 1 s from a 
% neuron firing an average of 10 spikes/s.  The second distribution is spike counts for 1 s
% come a neuron firing an average of 20 spikes/s.  The third distribution is
% spike counts for 2 s from a neuron firing an average of 10 spikes/s.  The
% final distribution is counts from 1 s combined from 2 neurons firing 10
% spikes/s. You should be able to see that the three final distributions all
% have a mean of 20 and are indistinguishable. 

counts = zeros(4, 250);

spikesPerS = 10;
durationS = 1.000;                  
spikes = makeSpikes(timeStepS, spikesPerS, durationS, 250);
counts(1, 1:250) = countSpikes(spikes, timeStepS);

spikesPerS = 20;
durationS = 1.000;                  
spikes = makeSpikes(timeStepS, spikesPerS, durationS, 250);
counts(2, 1:250) = countSpikes(spikes, timeStepS);

spikesPerS = 10;
durationS = 2.000;                  
spikes = makeSpikes(timeStepS, spikesPerS, durationS, 250);
counts(3, 1:250) = countSpikes(spikes, timeStepS);

spikesPerS = 10;
durationS = 1.000;                  
spikes1 = makeSpikes(timeStepS, spikesPerS, durationS, 250);
spikes2 = makeSpikes(timeStepS, spikesPerS, durationS, 250);
counts(4, 1:250) = countSpikes(spikes1 + spikes2, timeStepS);

plotSpikeCounts(counts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the end of the Poisson Tutorial.  The code below includes only
% functions that were used in the tutorial.  You should continue by running
% neuronalPerformance and following along in the Lesson.pdf in that folder. 


function counts = countSpikes(spikes, timeStepS, startS, endS)

if (nargin < 4)
    endS = length(spikes) * timeStepS;
end
if (nargin < 3)
    startS = 0;
end
trains = size(spikes, 1);
counts = zeros(1, trains);
startBin = startS / timeStepS + 1;
endBin = floor(endS / timeStepS);

for train = 1:trains
    counts(train) = sum(spikes(train, startBin:endBin));
end
end

function plotSpikeCounts(counts)

t = counts';

figure(3);
subplot(2, 1, 1);
hist(t, 0:max(max(t)) * 1.1);
s = axis;
axis([0 s(2) s(3) s(4)]);
subplot(2, 1, 2);
xlabel('Spike Counts');

m = mean(t);
v = var(t);
plot(m, v, 'bo');
s = axis;
axisLimit = max(s(2), s(4));
axis([0 axisLimit 0 axisLimit]);
hold on;
line([0, axisLimit], [0, axisLimit]);
hold off;
axis('square');
xlabel('Mean Spike Count');
ylabel('Spike Count Variance');

end

function rasterPlot(spikes, timeStepS)

figure(1);

times = [0:timeStepS:timeStepS * (length(spikes) - 1)];
axes('position', [0.1, 0.1, 0.8, 0.8]);
axis([0, length(spikes) - 1, 0, 1]);
trains = size(spikes, 1); 
ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
ticHeight = (1.0 - (trains + 1) * ticMargin) / trains;

for train = 1:trains
    spikeTimes = find(spikes(train, :) == 1);
    yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight);
    for i = 1:length(spikeTimes)
        line([spikeTimes(i), spikeTimes(i)], [yOffset, yOffset + ticHeight]);
    end
end

xlabel('Time (s)')
title('Raster plot of spikes');
end

function spikes = makeSpikes(timeStepS, spikesPerS, durationS, numTrains)

if (nargin < 4)
    numTrains = 1;
end
times = [0:timeStepS:durationS];
spikes = zeros(numTrains, length(times));
for train = 1:numTrains
    vt = rand(size(times));
    spikes(train, :) = (spikesPerS*timeStepS) > vt;
end
end

