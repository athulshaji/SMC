clear;clc

fm_1        = 440;

iterations  = 20000;
fs          = 44100;

N = 40;
x = 2 * rand(1,N) - 1;                                                     % generate noise centered at zero

fm_2 = fm_1 * 4;

if iterations > length(x)
    x = [x zeros(1,iterations - length(x))];
end

% Band pass filter parameters - - - - - - - - - - - - - - - - - - - - - - -
B_1 = 0.9;
[phi_1,R_1,theta_1,A0_1] = cmpt_modeparameters(B_1,fm_1,fs);
delaylength1 = floor(fs/fm_1);                                              % compute delay line length
delayline_1 = zeros(1,delaylength1+1);                                       % initialize delayline
dloffset_1 = length(delayline_1)-1;                                            % compute delay line offset

B_2 = 0.4;
[phi_2,R_2,theta_2,A0_2] = cmpt_modeparameters(B_2,fm_2,fs);
delaylength_2 = floor(fs/fm_2);                                              % compute delay line length
delayline_2 = zeros(1,delaylength_2+1);                                       % initialize delayline
dloffset_2 = length(delayline_2)-1;                                            % compute delay line offset

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
signal = zeros(1,iterations);                                              % initialize output signal
y_1 = 0;                                                                   % initialize new sample
xm1_1 = 0;                                                                 % x[n-1]
xm2_1 = 0;                                                                 % x[n-2]
ym1_1 = 0;                                                                 % y[n-1]
ym2_1 = 0;                                                                 % y[n-2]

xm1_2 = 0;                                                                 % x[n-1]
xm2_2 = 0;                                                                 % x[n-2]
ym1_2 = 0;                                                                 % y[n-1]
ym2_2 = 0;                                                                 % y[n-2]
for n = 1:iterations
    
    % mode 1
    y_1 = A0_1 * (x(n) - xm2_1 + 2 * R_1 * cos(theta_1) * ym1_1 - R_1 * R_1 * ym2_1);               % differential equation
       
    % update variables for next iteration
    xm2_1 = xm1_1;                                                             % x[n-2]
    xm1_1 = x(n);                                                            % x[n-1]
    
    ym2_1 = ym1_1;                                                             % y[n-2]
    ym1_1 = y_1;                                                               % y[n-1]
    
    delayline_1 = [y_1, delayline_1(1:dloffset_1)];
    
    % mode 2
    y_2 = A0_2 * (x(n) - xm2_2 + 2 * R_2 * cos(theta_2) * ym1_2 - R_2 * R_2 * ym2_2);               % differential equation
       
    % update variables for next iteration
    xm2_2 = xm1_2;                                                             % x[n-2]
    xm1_2 = x(n);                                                            % x[n-1]
    
    ym2_2 = ym1_2;                                                             % y[n-2]
    ym1_2 = y_2;                                                               % y[n-1]
    
    delayline_2 = [y_2, delayline_2(1:dloffset_2)];
    
    
    signal(n) = y_1 + y_2;
end
% plot(signal(1:200))