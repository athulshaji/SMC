clear;clc

fm        = 440;

iterations  = 20000;
fs          = 44100;
fb_gain     = 0.97;

N = ceil(fs/fm);
x = 2 * rand(1,N) - 1;                                                     % generate noise centered at zero

fm_2 = fm * 4;

if iterations > length(x)
    x = [x zeros(1,iterations - length(x))];
end

% Band pass filter parameters - - - - - - - - - - - - - - - - - - - - - - -
B = 0.9;
[phi,R,theta,A0,d,delayline] = cmpt_modeparameters(B,fm,fs);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
signal = zeros(1,iterations);                                              % initialize output signal
y = 0;                                                                     % initialize new sample
xm1 = 0;                                                                   % x[n-1]
xm2 = 0;                                                                   % x[n-2]
ym1 = 0;                                                                   % y[n-1]
ym2 = 0;                                                                   % y[n-2]

for n = 1:iterations
    
%     y = A0 * ((x(n)+delayline(d+1)) - xm2 + 2 * R * cos(theta) * ym1 - R * R * ym2);        % differential equation
    
    y = A0 * (x(n) - xm2 + 2 * R * cos(theta) * ym1 - R * R * ym2);        % differential equation
    
    out = y + fb_gain * delayline(d);
    delayline = [out,delayline(1:d-1)];
    
    % update variables for next iteration
    xm2 = xm1;                                                             % x[n-2]
    xm1 = x(n);                                                            % x[n-1]
    
    ym2 = ym1;                                                             % y[n-2]
    ym1 = y;                                                               % y[n-1]
    
    
% 	x(n+d) = x(n+d) + delayline(d);
    signal(n) = out;
%     signal(n) = delayline(1);
end
plot(signal)











