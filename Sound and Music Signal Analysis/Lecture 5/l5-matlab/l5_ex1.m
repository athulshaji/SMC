%% Exercise 1:
% Implement one of the methods in Christensen - Pitch estimation for dummies (2015). 
% Generate a synthetic harmonic signal (assume a fixed model order) and add white Gaussian noise 
% resulting in different SNRs (0 dB, 10 dB, 20 dB), 
% and test your implementation using the signals you generated (calculate the MSE). 


%% Generate signal and add noise
desired_snr = 0;                                                           % desired SNR relation between signal and noise
fs = 44100;                                                                % sampling frequency
f = 440;                                                                   % fundamental frequency
A = 0.8;                                                                   % amplitude of the signal
duration = 1;                                                              % duration of the signal in seconds
N = fs * duration;                                                         % size of signal in samples
t = 0:1/fs:duration-1/fs;                                                  % time array of the signal

signal = A * sin(2*pi*f*t);                                                % create signal
x = awgn(signal,desired_snr,'measured');                                   % add gaussian noise to the signal

%% Pitch estimation
lowestfreq = 100;                                                          % lowest frequency considered
highestfreq = 1000;                                                        % highest frequency considered
tau = round(fs/highestfreq):10:round(fs/lowestfreq);                       % all periods considered

MSE = zeros(length(tau),1);                                                % initialize Mean Squared Error array
for i = 1:length(tau)                                                      % for all the periods
    e = zeros(1,N-tau(i));                                                 % initialize error array
    for n = 1:N-tau(i)                                                     % for all the samples of the signal ?
        e(n) = x(n) - x(n+tau(i));                                         % compute the error between the n sample
    end                                                                    % and the sample delayed by tau (period)
    MSE(i) = (1/(N-tau(i))) * sum(power(e,2));                             % compute mean squared error for this tau (period)
end
[y,i] = min(MSE);                                                          % find the minimum mean squared error
estimatedfreq = fs/tau(i)                                                  % estimated frequency: period to frequency