% Lab 2: Build an integrate-and-fire model neuron and observe its spiking
% for various levels of injected current
clear all; %clear all variables
close all; %close all open figures
%DEFINE PARAMETERS
dt = 0.1; %time step [ms]
t_end = 500; %total time of run [ms]
t_StimStart = 100; %time to start injecting current [ms]
t_StimEnd = 400; %time to end injecting current [ms]
E_L = -70; %resting membrane potential [mV]
V_th = -55; %spike threshold [mV]
V_reset = -75; %value to reset voltage to after a spike [mV]
V_spike = 20; %value to draw a spike to, when cell spikes [mV]
R_m = 10; %membrane resistance [MOhm]
tau = 10; %membrane time constant [ms]
%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect = 0:dt:t_end; %will hold vector of times
V_vect = zeros(1,length(t_vect)); %initialize the voltage vector
 %initializing vectors makes your code run faster!
V_plot_vect = zeros(1,length(t_vect)); %pretty version of V_vect to be plotted, that displays a spike
 % whenever voltage reaches threshold
%INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m
PlotNum=0;
I_Stim_vect = 1.43:0.04:1.63; %magnitudes of pulse of injected current [nA]
for I_Stim = I_Stim_vect %loop over different I_Stim values
 PlotNum = PlotNum + 1;
 i = 1; % index denoting which element of V is being assigned
 V_vect(i)= E_L; %first element of V, i.e. value of V at t=0
 V_plot_vect(i) = V_vect(i); %if no spike, then just plot the actual voltage V
 I_e_vect = zeros(1,t_StimStart/dt); %portion of I_e_vect from t=0 to t=t_StimStart
 I_e_vect = [I_e_vect I_Stim*ones(1,1+((t_StimEnd-t_StimStart)/dt))]; %add portion from
 % t=t_StimStart to t=t_StimEnd
 I_e_vect = [I_e_vect zeros(1,(t_end-t_StimEnd)/dt)]; %add portion from
 % t=t_StimEnd to t=t_end
 NumSpikes = 0; %holds number of spikes that have occurred
 for t=dt:dt:t_end %loop through values of t in steps of dt ms
 V_inf = E_L + I_e_vect(i)*R_m; %value that V_vect is exponentially
 %decaying towards at this time step

 %next line does the integration update rule
 V_vect(i+1) = V_inf + (V_vect(i) - V_inf)*exp(-dt/tau);
 %if statement below says what to do if voltage crosses threshold
 if (V_vect(i+1) > V_th) %cell spiked
 V_vect(i+1) = V_reset; %set voltage back to V_reset
 V_plot_vect(i+1) = V_spike; %set vector that will be plotted to show a spike here
 NumSpikes = NumSpikes + 1; %add 1 to the total spike count
 else %voltage didn't cross threshold so cell does not spike
 V_plot_vect(i+1) = V_vect(i+1); %plot the actual voltage
 end
 i = i + 1; %add 1 to index, corresponding to moving forward 1 time step
 end
 AveRate_vect(PlotNum) = 1000*NumSpikes/(t_StimEnd - t_StimStart); %gives average firing
 %rate in [#/sec = Hz]
 %MAKE PLOTS
 figure(1)
 subplot(length(I_Stim_vect),1,PlotNum)
 plot(t_vect, V_plot_vect);
 if (PlotNum == 1)
 title('Voltage vs. time');
 end
 if (PlotNum == length(I_Stim_vect))
 xlabel('Time in ms');
 end
 ylabel('Voltage in mV');
end %for I_Stim
%COMPARE R_AVE TO R_ISI
I_threshold = (V_th - E_L)/R_m; %current below which cell does not fire
I_vect_long = (I_threshold+0.001):0.001:1.8; %vector of injected current for producing theory plot
r_isi = 1000./(tau*log((V_reset - E_L - I_vect_long*R_m)./(V_th - E_L - I_vect_long*R_m)));
figure(2)
plot(I_vect_long,r_isi)
hold on
plot(I_Stim_vect,AveRate_vect,'ro')
title('Comparison of r_{isi} vs. I_e and r_{ave} vs. I_e')
xlabel('Injected current (nA)')
ylabel('r_{isi} or r_{ave} (Hz)')
hold off %to ensure that doesn't keep this data the next time you want to plot something