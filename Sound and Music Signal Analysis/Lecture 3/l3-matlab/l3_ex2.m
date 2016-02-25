%% Lecture 2 - exercise 2
% Make a simple source filter model and use it for speech synthesis.
% Use the speech signal fa.wav. Process the signal in blocks of 25 ms. 
% Use lpc() with an order of 10 to obtain the prediction error (or excitation) for each segment (use filter()).
% Replace the the prediction error signal with one generated from white 
% Gaussian noise (use randn()) having the same variance as the prediction error. 
% Filter the synthetic signal with the all-pole filter with filter() to resynthesize the speech. 
% Remember to take care of the filter states (see 'help filter'). 

% Listen to the result and compare it to the original signal. How does it sound? 
% Look at the spectrograms of the two signals. What is different?
clearvars -except x fs;clc;close all
if ~exist('x','var')
    [x,fs] = audioread('../../data/fa.wav');                               % desired signal
end
%%
segmentlength = 0.025;                                                     % segment length in seconds (25 ms)
sgmtsamples = segmentlength * fs;                                          % segment length in samples
lpc_order = 10;                                                            % LPC order
idx = 1:sgmtsamples;                                                       % indices of sound segment
outputsignal = zeros(size(x));                                               % initialize output sound
predictionError = zeros(size(x));                                                % 
lpc_coeffs = [];

% Use lpc() with an order of 10 to obtain the prediction error (or excitation) for each segment (use filter())
while idx(end) < length(x) - sgmtsamples
    sgmt = x(idx); % get a segment of the sound
    
    [A,E] = lpc(sgmt, 10);
    
    % Replace the the prediction error signal with one generated from white 
    % Gaussian noise (use randn()) having the same variance as the prediction error. 
    noise = sqrt(E) .* randn(length(idx),1);
    
    % Filter the synthetic signal with the all-pole filter with filter() to resynthesize the speech. 
    % Remember to take care of the filter states (see 'help filter').
    outputsignal(idx) = filter(1,A,noise);
    
    idx = idx + sgmtsamples;
end