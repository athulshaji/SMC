%% Exercise 1:
% Implement one of the methods in Christensen - Pitch estimation for dummies (2015). 
% Generate a synthetic harmonic signal (assume a fixed model order) and add white Gaussian noise 
% resulting in different SNRs (0 dB, 10 dB, 20 dB), 
% and test your implementation using the signals you generated (calculate the MSE). 
clear; clc
%% Generate signal and add noise
desired_snr = 30;                                                           % desired SNR relation between signal and noise
fs = 8000;                                                                % sampling frequency
f = 440;                                                                   % fundamental frequency
L = 10;                                                                     % number of harmonics
A = 0.8;                                                                   % amplitude of the signal
duration = 5;                                                              % duration of the signal in seconds
N = fs * duration;                                                         % size of signal in samples
t = 0:1/fs:duration-1/fs;                                                  % time array of the signal

signal = zeros(size(t));                                                   % initialize signal
for n = 1:L
    signal = signal + (A/n) * sin(2*pi*f*n*t);                             % create signal
end
signal = signal/max(signal);
x = awgn(signal,desired_snr,'measured');                                   % add gaussian noise to the signal

%% Pitch estimation with Comb filtering
maxperiod = 50;
[mse, estimatedperiod] = l5_combfilter(x,maxperiod);
estimatedfreq = fs/estimatedperiod