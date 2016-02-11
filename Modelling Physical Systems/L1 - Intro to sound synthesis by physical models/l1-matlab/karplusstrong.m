% define input parameters
freqHz = 500;
iterations = 200000;
fs = 44100;
lpf_a = 0.5;                                            % low pass filter factor (decaying factor)

% Karplus strong function - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

N = floor(fs/freqHz);                                   % compute delay line length
x = 2 * rand(1,N) - 1;                                  % generate noise centered at zero

delayline = [zeros(1,N+1)];                             % initialize delayline
dloffset = length(delayline)-1;                         % compute delay offset

if iterations > length(x)
    diff = iterations - length(x);
    x = [x zeros(1,diff)];
end

y = 0;                                                  % initialize new sample
signal = zeros(1,iterations);                           % initialize output signal
n = 1;                                                  % initialize iterator

while n < iterations
    y = x(n) + lpf_a * delayline(N) + (1-lpf_a) * delayline(N+1);       % karplus strong with variable low pass filter
    
    delayline = [y, delayline(1:dloffset)];             % shift delay line and write new sample
    signal(n) = y;                                      % append new sample to output signal
    
    n = n+1;                                            % increase iterator
end
figure;plot(signal);title('Simple Karplus Strong algorithm')   % show signal
% soundsc(signal, fs)                                     % play it