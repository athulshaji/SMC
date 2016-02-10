% define input parameters
freqHz = 500;
iterations = 20000;
fs = 44100;
lpf_a = 0.9;                                            % low pass filter factor (decaying factor)

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
signal = 0;                                             % initialize output signal
n = 1;                                                  % initialize iterator

while n < iterations
    
    y = x(n) + 0.5*(delayline(N) + delayline(N+1));     % karplus strong delay
    
    y = lpf_a * y + (1 - lpf_a) * delayline(1);         % apply low pass filter to control decay
    
    delayline = [y, delayline(1:dloffset)];             % shift delay line and write new sample
    signal = [signal y];                                % append new sample to output signal
    
    n = n+1;                                            % increase iterator
end
% plot(signal);title('Simple Karplus Strong algorithm')   % show signal
soundsc(signal, fs)                                     % play it