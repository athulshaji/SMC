%% Exercise 2:
% 
% Use the joint_anls.m in the toolbox (link on course webpage) to estimate 
% the pitch and the model order of your signal from exercise 1. 
% Also try the estimator on a real signal from the IOWA database (available 
% at http://theremin.music.uiowa.edu) as input signal. Describe your results. 
% An example of how to use the joint_anls.m function can be found in the toolbox (example1.m).
clear;clc
addpath('pitch/pitch')

%% Generate signal and add noise
desired_snr = 10;                                                          % desired SNR relation between signal and noise
fs = 8000;                                                                 % sampling frequency
f = 440;                                                                   % fundamental frequency
L = 10;                                                                    % number of harmonics
A = 0.8;                                                                   % amplitude of the signal
duration = 0.05;                                                            % duration of the signal in seconds
N = fs * duration;                                                         % size of signal in samples
t = 0:1/fs:duration-1/fs;                                                  % time array of the signal

signal = zeros(size(t));                                                   % initialize signal
for n = 1:L
    signal = signal + (A/n) * sin(2*pi*f*n*t);                             % create signal
end
signal = signal/max(signal);
x = awgn(signal,desired_snr,'measured');                                   % add gaussian noise to the signal
x = analytic(x)';                                                          % make the signal complex

%% Estimate the pitch and model order of a synthesized signal with joint_anls.m
%   x        input signal (assumed complex)
%   w0_lim   search invertval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 4096)
w0_lim = 2*pi*[1e-3 0.25];
F = 4096;
[w0,L] = joint_anls(x, 2*pi*[1e-3 0.25], 4096);
%% Estimate the pitch and model order of a real signal with joint_anls.m
%   x        input signal (assumed complex)
%   w0_lim   search invertval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 4096)
[x ~] = audioread('audio/Vibraphone.dampen.mf.C6F6.aif');
x = analytic(x)';                                                          % make the signal complex
w0_lim = 2*pi*[1e-3 0.25];
F = 4096;
[w0,L] = joint_anls(x, 2*pi*[1e-3 0.25], 4096)