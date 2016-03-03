%% Exercise 1:
% Implement one of the methods in Christensen - Pitch estimation for dummies (2015). 
% Generate a synthetic harmonic signal (assume a fixed model order) and add white Gaussian noise 
% resulting in different SNRs (0 dB, 10 dB, 20 dB), 
% and test your implementation using the signals you generated (calculate the MSE). 

%% Generate signal and add noise
desired_snr = 0;                                                           % desired SNR relation between signal and noise
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

%% Pitch estimation with Comb filtering
lowestfreq = 100;                                                          % lowest frequency considered
highestfreq = 1000;                                                        % highest frequency considered
tau = round(fs/highestfreq):1:round(fs/lowestfreq);                        % all periods considered

MSE = zeros(length(tau),1);                                                % initialize Mean Squared Error array
for i = 1:length(tau)                                                      % for all the periods
    x_delayed = x(tau(i)+1:end);                                           % delay signal
    e = x(1:end-tau(i)) - x_delayed;                                       % compute the error between the n sample
    MSE(i) = (1/(N-tau(i))) * sum(power(e,2));                             % compute mean squared error for this tau (period)
end
[y,i] = min(MSE);                                                          % find the minimum mean squared error
estimatedfreq = fs/tau(i)                                                  % estimated frequency: period to frequency