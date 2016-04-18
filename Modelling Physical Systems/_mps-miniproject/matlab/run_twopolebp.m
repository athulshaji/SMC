% generate noise
fs = 44100;                                                                % sampling frequency
duration = 1;                                                              % duration of noise (seconds)
% N = floor(fs/freqHz);           % compute delay line length
N = fs*duration;                                                           % duration of noise (samples)
x = 2 * rand(1,N) - 1;                                                     % generate noise centered at zero
% soundsc(x,fs)

% filter noise
fm = 4000;                                                                 % mode frequency
B = 0;                                                                     % bandwidth
y = twopole_bandpass(x, fm,fs,B);

plot(y)