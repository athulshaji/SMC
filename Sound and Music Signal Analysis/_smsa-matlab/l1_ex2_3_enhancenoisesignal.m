clearvars -except x v; clc
%% Load files
if ~exist('x','var')
    [x,fs] = audioread('../../data/twoMaleTwoFemale20Seconds.wav');        % desired signal
    [v,fs] = audioread('../../data/babble30Seconds.wav');                  % noise signal
end
%% Exercise 3:

desired_snr = 0;                                                           % define desired SNR
L = 40;                                                                    % define filter length
forgetFSignal = 0.2;%985;                                                  % forgetting factor for speech
forgetFNoise = 0.2;%995;                                                   % forgetting factor for noise

y = mixatsnr(x,v,desired_snr);                                             % mix signals at desired SNR

frameIndices = 1:L;                                                        % initialize iterator indices
window = hann(L,'periodic');                                               % create window to overlap
CovMatSignal = eye(L);                                                     % initialize signal covariance matrix
CovMatNoise = eye(L);                                                      % initialize noise covariance matrix
filteredy = zeros(size(y));                                                % initialize output signal
regulPar = 1e-10;                                                          % initialize regularization factor

while frameIndices(end) <= length(y)
    
    signalFrame = y(frameIndices);                                         % define signal frame
    noiseFrame = v(frameIndices);                                          % define noise frame
    
    % compute covariance matrices of signal and noise
    CovMatSignal = forgetFSignal * CovMatSignal + (1-forgetFSignal) * (signalFrame * signalFrame');
    CovMatNoise = forgetFNoise * CovMatNoise + (1-forgetFNoise) * (noiseFrame * noiseFrame');
    
    % regularize signal covariance matrix
    CovMatSignal = CovMatSignal*(1-regulPar)+...                           
        (regulPar)*trace(CovMatSignal)/(L)*...
        eye(L);

    % overlap 50% of frame length
	if mod(frameIndices(end),L/2) == 0,
        Hw = eye(L) - CovMatSignal\CovMatNoise; % compute filter for frame.
                                                % X = A\B is the solution to the equation AX = B
        filteredy(frameIndices) = filteredy(frameIndices) + window .* (Hw' * signalFrame);    % filter signal
	end 
    
    frameIndices = frameIndices + 1;                                       % increase iterator
end