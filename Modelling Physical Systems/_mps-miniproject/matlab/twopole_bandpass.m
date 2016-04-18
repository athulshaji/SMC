% Implementation of the band pass filter, equation 2 of the article:
% Essl, Georg, et al. "Theory of banded waveguides." Computer Music Journal 28.1 (2004): 37-50.
%
% Frequency response:
%                      1 - z^(-2)
% H(z) = ---------------------------------------
%          1 - (2*R*cos(theta))*z^-1 + R^2*z^-2
%
% freqz([1,0,-1],[1,-2*R*cos(theta),R*R])
%
% Differential equation:
% y[n] =  x[n] - x[n-2]  +  2 * R * cos(theta) * y[n-1]  -  R^2 * y[n-2]
%
% where R and theta are free parameters of the poles that relate to bandwidth B, 
% center frequency w, and gain A0 in the following way (Steiglitz 1996):
% R ? 1 - B/2;
%                2*R
% cos(theta) = --------- * cos(phi)
%               1 + R^2
% A0 = (1-R^2) * sin(theta)
%
% phi = 2*pi*(fm/fs)
%
% phi:  resonance frequency
% fm:   frequency of a mode
% fs:   sampling frequency
% A0:   gain of band pass
% B:    bandwidth, chosen to reject other modes
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function y = twopole_bandpass(x, fm,fs,B)

phi = 2 * pi * (fm/fs);                                                    % resonance frequency

% parameters for the poles
R = 0.99-B/2;                                                              % filter bandwidth
theta = acos((2 * R * cos(phi)) / (1 + R * R));                            %

A0 = (1 - R * R) * sin(theta);                                             % how to use this?

y = zeros(size(x));                                                        % initialize output signal
xm1 = 0;                                                                   % x[n-1]
xm2 = 0;                                                                   % x[n-2]
ym1 = 0;                                                                   % y[n-1]
ym2 = 0;                                                                   % y[n-2]

% Apply band pass filter
for n = 1:length(x)
    y(n) = x(n) - xm2 + 2 * R * cos(theta) * ym1 - R * R * ym2;            % differential equation
    
    % update variables for next iteration
    xm2 = xm1;                                                             % x[n-2]
    xm1 = x(n);                                                            % x[n-1]
    
    ym2 = ym1;                                                             % y[n-2]
    ym1 = y(n);                                                            % y[n-1]
    
end
% y should be the same as yy:
% yy = filter([1,0,-1],[1,-2*R*cos(theta),R*R],x);
end