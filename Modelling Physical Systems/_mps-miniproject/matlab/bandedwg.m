clear;clc
nmodes      = 3;                                                           % number of modes
f0          = 234.8;%440;                                                         % fundamental frequency
% fm          = f0*[1,1.97,2.97];%(1:nmodes);                                               % frequency for each mode
fm          = f0*[1,1.99,2.98];%(1:nmodes);                                               % frequency for each mode
% decaytimes  = [2.4,1.7,0.9];
decaytimes  = [3.3,2.5,1.7];
% fm          = [443];%[139,322,527,689,1033,1324,1518];
% fm          = [527,689,1033,1324,1518];
% iterations  = 100000;

fs          = 44100;                                                       % sampling frequency
fb_gain     = 0.995;                                                        % feedback gain. < 1
% a = 0.9;
% create excitation - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% noise:
N = ceil(fs/fm(1));
x = 2 * rand(1,N) - 1;

% cosine impact:
Fm = 0.6;
forceduration = 0.002;
t = 0:1/fs:forceduration;
x = Fm * (1 - cos(2*pi*t/forceduration));
% 
% % recorded impulse:
% % ximpulse = audioread('littlehang1note_443.wav')';% x = x(1:47500)';
ximpulse = audioread('hangimpulse.wav')';ximpulse = ximpulse(1:0.4*fs);
x = conv(x,ximpulse);

% compute number of iterations
iterations  = floor(max(decaytimes)*fs);%100000;

if iterations > length(x)
    x = [x zeros(1,iterations - length(x))];
else
    iterations = length(x);
end

% Band pass filter parameters - - - - - - - - - - - - - - - - - - - - - - -
% compute Q factors (B)
B           = decaytimes .* (2 * pi * fm);
% B           = [7,5,5].*ones(1,nmodes);
phi         = zeros(1,nmodes);
R           = zeros(1,nmodes);
theta       = zeros(1,nmodes);
A0          = zeros(1,nmodes);
d           = zeros(1,nmodes);
delayline   = zeros(nmodes,floor(fs/fm(1)));
for m = 1:nmodes                                                           % compute parameters for each mode :
    phi(m) = 2 * pi * (fm(m)/fs);                                          % resonance frequency
    % parameters for the poles
%     R(m) = 0.99-B(m)/2;                                                    % filter bandwidth
    R(m) = exp(-pi*B(m)*(1/fs));
    theta(m) = acos((2 * R(m) * cos(phi(m))) / (1 + R(m) * R(m)));         %
    A0(m) = (1 - R(m) * R(m)) * sin(theta(m));                             % gain
    d(m) = floor(fs/fm(m));                                                % compute delay line length
end
% m = 1;freqz([1,0,-1],[1,-2*R(m)*cos(theta(m)),R(m)*R(m)]);
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
        
%         if ~mod(n,10)
%             plot(signal);
%             title(num2str(n))
%             drawnow
%         end
    end
    
    signal(n) = sum(out);
    
end
tt = 0:1/fs:(length(signal)-1)/fs;
figure;plot(tt,signal)
soundsc(signal,fs)

% better impact