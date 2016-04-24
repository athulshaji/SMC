clear;clc
nmodes      = 3;                                                           % number of modes
f0          = 440;                                                         % fundamental frequency
fm          = f0*(1:nmodes);                                               % frequency for each mode
iterations  = 20000;
fs          = 44100;                                                       % sampling frequency
fb_gain     = 0.97;                                                        % feedback gain. < 1

% create excitation - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N = ceil(fs/f0);
x = 2 * rand(1,N) - 1;                                                     % generate noise centered at zero

if iterations > length(x)
    x = [x zeros(1,iterations - length(x))];
end

% Band pass filter parameters - - - - - - - - - - - - - - - - - - - - - - -
B           = 0.01*ones(1,nmodes);
phi         = zeros(1,nmodes);
R           = zeros(1,nmodes);
theta       = zeros(1,nmodes);
A0          = zeros(1,nmodes);
d           = zeros(1,nmodes);
delayline   = zeros(nmodes,floor(fs/fm(1)));
for m = 1:nmodes
    [phi(m),R(m),theta(m),A0(m),d(m)] = cmpt_modeparameters(B(m),fm(m),fs);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
signal  = zeros(1,iterations);                                             % initialize output signal
y       = zeros(1,nmodes);                                                 % initialize new sample
xm1     = zeros(1,nmodes);                                                 % x[n-1]
xm2     = zeros(1,nmodes);                                                 % x[n-2]
ym1     = zeros(1,nmodes);                                                 % y[n-1]
ym2     = zeros(1,nmodes);                                                 % y[n-2]
out     = zeros(1,nmodes);                                                 % output samples
for n = 1:iterations
    
    for m = 1:nmodes                                                       % for each mode
        y(m) = A0(m) * (x(n) - xm2(m) + 2 * R(m) * cos(theta(m)) * ym1(m) - R(m) * R(m) * ym2(m)); % band pass filter differential equation
        
        out(m) = y(m) + fb_gain * delayline(m,d(m));                       % output sample: sum sample to delayed sample
        delayline(m,1:d(m)) = [out(m),delayline(m,1:d(m)-1)];              % store output sample to circular delay
        
        % update variables for next iteration
        xm2(m) = xm1(m);                                                   % x[n-2]
        xm1(m) = x(n);                                                     % x[n-1]
        
        ym2(m) = ym1(m);                                                   % y[n-2]
        ym1(m) = y(m);                                                     % y[n-1]
    end
    
    signal(n) = sum(out);
    
end
plot(signal)
soundsc(signal,fs)