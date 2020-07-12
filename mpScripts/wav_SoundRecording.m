% wav_SoundRecording.m
% Generate and save sound files for two frequency inputs

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170511
% Ignore Warning about clipping


clear all
close all
clc

% Frequency inputs

   f1 = 8000;  
   filename = 'wav_S8000.wav';
   f2 = f1;
   
   %f2 = 3003;

% Calculate Waveform   
  fs = 8*22050;                   % sample frequency (Hz)
  d = 4.0;                      % duration (s)
  n = fs * d;                   % number of samples
  t = (1:n) / fs;               % sound data preparation
  % s = sin(2 * pi * f1 * s);   % pure tone
  s = sin(2 * pi * f1 * t)+ sin(2 * pi * f2 * t);
  s = s./max(s);

  sound(s, fs);                   % Generate sound
  pause(d + 0.5);                 % waiting for sound end

% Save wav file to disk
   
   audiowrite(filename,s,fs);
   
   
   