clearvars -except x v; clc
%% Load files
if ~exist('x','var')
    [x,fs] = audioread('../../data/twoMaleTwoFemale20Seconds.wav');                 % desired signal
    [v,fs] = audioread('../../data/babble30Seconds.wav');                           % noise signal
end
%% Exercise 3:

desired_snr = 0;                                                           % define desired SNR
L = 40;                                                                    % define filter length
forgetFSignal = 0.985;                                                     % forgetting factor for speech
forgetFNoise = 0.995;                                                      % forgetting factor for noise

y = mixatsnr(x,v,desired_snr);                                             % mix signals at desired SNR

frameIndices = 1:L;
window = hann(L,'periodic');
CovMatSignal = eye(L);                                                     % initialize signal covariance matrix
CovMatNoise = eye(L);                                                      % initialize noise covariance matrix
filteredy = zeros(size(y));                                                % initialize output signal

while frameIndices(end) <= length(y)
    
    signalFrame = y(frameIndices);                                              % define frame
    noiseFrame = v(frameIndices);                                               % define frame
    
    % compute covariance matrices of signal and noise
    CovMatSignal = forgetFSignal * CovMatSignal + (1-forgetFSignal) * (signalFrame * signalFrame');
    CovMatNoise = forgetFNoise * CovMatNoise + (1-forgetFNoise) * (noiseFrame * noiseFrame');
       
	if mod(frameIndices(end),L/2) == 0,
        Hw = eye(L) - inv(CovMatSignal) * CovMatNoise;                          % compute filter for frame
        filteredy(frameIndices) = filteredy(frameIndices) + window .* (Hw' * signalFrame);    % filter signal
	end 
    
    frameIndices = frameIndices + 1;
end