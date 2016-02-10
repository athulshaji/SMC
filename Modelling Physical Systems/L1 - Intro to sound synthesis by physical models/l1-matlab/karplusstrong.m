freqHz = 1000;
iterations = 20000;
fs = 44100;

N = floor(fs/freqHz);           % compute delay line length
x = 2 * rand(1,N) - 1;          % generate noise centered at zero

delayline = [zeros(1,N+1)];     % initialize delayline
dloffset = length(delayline)-1; % compute delay offset

if iterations > length(x)
    diff = iterations - length(x);
    x = [x zeros(1,diff)];
end

y = 0;                          % initialize new sample
signal = 0;                     % initialize output signal
n = 1;                          % initialize iterator

while n < iterations
    y = x(n) + 0.5*(delayline(N) + delayline(N+1));
    delayline = [y, delayline(1:dloffset)];
    signal = [signal y];
    n = n+1;
end
plot(signal);title('Simple Karplus Strong algorithm')

% soundsc(signal, fs)