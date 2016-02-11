% define input parameters
freqHz = 500;
iterations = 50000;
fs = 44100;
lpf_a = 0.1;                                            % low pass filter factor (decaying factor)
apf_g = 0.5;                                            % all pass filter factor (harmonicity factor)
bi_filename = 'p_dobro_2.wav';              
[bodyimpulse,vs] = audioread(bi_filename);              % read body impulse response            

% Karplus strong function - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

N = floor(fs/freqHz);                                   % compute delay line length (samples)

% compute exciter signal:
% exciterSignal = 2 * rand(1,N) - 1;                      % generate noise centered at zero
NPp = N*0.5;                                             % NPp = pick positon in samples
exciterSignal = [[0:NPp]/NPp,(N-[(NPp+1):N])/(N-NPp)]*2-1;% triangular funciton
exciterSignal = conv(exciterSignal,bodyimpulse,'same');   % convolve exciter signal with body impulse:


delayline = [zeros(1,N+1)];                             % initialize delayline
dloffset = length(delayline)-1;                         % compute delay offset

if iterations > length(exciterSignal)                   % exciter signal has to be long 
    diff = iterations - length(exciterSignal);          % as the number of iterations
    exciterSignal = [exciterSignal zeros(1,diff)];      %
end

y = 0;                                                  % initialize new sample
ylpf = 0;                                               % initialize new sample low pass filtered
signal = zeros(1,iterations);                           % initialize output signal
n = 1;                                                  % initialize iterator

while n < iterations
    prev_ylpf = ylpf;
    ylpf = exciterSignal(n) + lpf_a * delayline(N) ...  % karplus strong with variable low pass filter
        + (1-lpf_a) * delayline(N+1);       
    
    y = - apf_g * ylpf + prev_ylpf + apf_g * delayline(1); % apply all pass filter to control harmonicity
    
    delayline = [y, delayline(1:dloffset)];             % shift delay line and write new sample
    signal(n) = y;                                      % append new sample to output signal
    
    n = n+1;                                            % increase iterator
end
figure;plot(signal);title('Simple Karplus Strong algorithm')   % show signal
soundsc(signal, fs)                                     % play it