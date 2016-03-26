% Exercise 1:
% Implement the MAP estimator (2.78) in Christensen and Jakobsson - Multi-Pitch Estimation (2009), 
% and use it on a synthetic signal for different numbers of harmonics, segment lengths and SNRs.
% Plot the results.
clear;clc

%% Generate signal and add noise
desired_snr = 10;                                                          % desired SNR relation between signal and noise
fs = 8000;                                                                 % sampling frequency
f = 440;                                                                   % fundamental frequency
L = 10;                                                                    % number of harmonics
A = 0.8;                                                                   % amplitude of the signal
duration = 0.5;                                                            % duration of the signal in seconds
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
w_array = [400,500];%2*pi*[0.001 0.1];
M = 1000;
N = M;
sgmt = x(1:M)';
% estimated_w = zeros(1,length(w_array));
estim_w_array = [];
for w = w_array(1):w_array(end)
    Z = l5_createVandermondeMat(w,L,M);
    sigma_n = power(norm(sgmt - Z*pinv(Z)*sgmt,2),2)/N;
    estim_w = N * log(sigma_n) + 1.5*log(N)+w*log(N);
    
    estim_w_array = [estim_w_array estim_w];
end
[val idx] = min(estim_w_array);
estim_w_array(idx)