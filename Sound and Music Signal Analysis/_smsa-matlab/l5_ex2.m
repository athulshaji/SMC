%% Exercise 2:
% Implement the NLS estimator (2.32) in Christensen and Jakobsson - Multi-Pitch Estimation (2009) 
% (assume complex signals, use, e.g., analytic.m) and test the estimator on your signals from exercise 1. 
% Describe the differences you observe for the various noise levels.
clear; clc
% addpath('../../_smsa-matlab')                                            % add the path of the smsa library

%% Generate signal and add noise
desired_snr = 10;                                                           % desired SNR relation between signal and noise
fs = 44100;                                                                % sampling frequency
f = 440;                                                                   % fundamental frequency
L = 10;                                                                     % number of harmonics
A = 0.8;                                                                   % amplitude of the signal
duration = 1;                                                              % duration of the signal in seconds
N = fs * duration;                                                         % size of signal in samples
t = 0:1/fs:duration-1/fs;                                                  % time array of the signal

signal = zeros(size(t));                                                   % initialize signal
for n = 1:L
    signal = signal + (A/n) * sin(2*pi*f*n*t);                             % create signal
end
signal = signal/max(signal);
x = awgn(signal,desired_snr,'measured');                                   % add gaussian noise to the signal
x = analytic(x)';                                                          % make the signal complex
%% 
% L = 5;                                                                     % number of harmonics
w_array = 400:500;
M = 1000;
sgmt = x(1:M);
estimated_w = zeros(1,length(w_array));
for i = 1:length(w_array)
    w = w_array(i);
    Z = l5_createVandermondeMat(w,L,M);
    estimated_w(i) = power(norm(Z'*sgmt',2),2);
end
[val idx] = max(estimated_w);
w_array(idx)