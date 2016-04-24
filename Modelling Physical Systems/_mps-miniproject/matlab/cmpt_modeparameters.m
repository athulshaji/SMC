function [phi,R,theta,A0,d] = cmpt_modeparameters(B,fm,fs)
phi = 2 * pi * (fm/fs);                                                    % resonance frequency
% parameters for the poles
R = 0.99-B/2;                                                              % filter bandwidth
theta = acos((2 * R * cos(phi)) / (1 + R * R));                            %
A0 = (1 - R * R) * sin(theta);                                             % gain
d = floor(fs/fm);                                                          % compute delay line length
% delayline = zeros(1,d);                                                  % initialize delayline
end