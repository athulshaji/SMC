%% Exercise 1 (Non-negative Matrix Factorization):
% 
% In this exercise, you will make your own implementation of simple NMF. 
% Compute FFTs (of the same size as the input segment) in steps of 25 ms and save the magnitude spectra in a matrix Y. 
% Implement the multiplicative updates in (7) and (8) in Virtanen et al. for 50 atoms. 
% Run the updates 100 times and plot the X and A matrices with imagesc in each iteration. 
% How does it look? After convergence, resynthesize the signal from the A and X matrices 
% by combining them with the phase of the original signal. 
% 
% Hint: If you use audio, I advice you downsample to 8kHz to speed everything up and make it simpler to plot things.

%% Import smsa library
clearvars -except x fs;clc;close all
addpath('../../_matlabscripts')
disp(smsa_library)

if ~exist('x','var')
%     [x,fs] = audioread('../../data/twoMaleTwoFemale20Seconds.wav');        % desired signal
    [x,fs] = audioread('xylophone.wav');        % desired signal
end
[x,fs]=myresample(x,fs,8000);

%% STFT Analysis: compute spectrogram
disp('STFT Analysis')
segmentlength = 0.025;                                                     % 30 ms
N = segmentlength * fs;                                                    % fft size
window = sqrt(hanning(N-1));                                               % initialize analysis window
soundlength = length(x);                                                   % compute sound lenght in samples

[mX, pX] = stft_Analysis(x,N,window);                                      % compute magnitude and phase spectrum as matrices

%% Non-negative Matrix Factorization
disp('Non-negative Matrix Factorization')
K = 50;                                                                    % number of atoms
iterations = 100;                                                          % iterations to find atoms and activation matrices
atomMatrix = rand(N/2+1,K);                                                % initialize atom matrix
activationMat = rand(K,size(mX,2));                                        % initialize activation matrix

it = 0;
while it <= iterations
    atomMatrix = atomMatrix .* (((mX ./ (atomMatrix * activationMat+eps)) * activationMat')./(ones(size(mX)) * activationMat'));
    activationMat = activationMat .* ((atomMatrix' * (mX./(atomMatrix*activationMat+eps))) ./ (atomMatrix' * ones(size(mX))));
    
    imagesc(atomMatrix);title('atom Matrix')
    drawnow
    
    it = it + 1;
end

mX = atomMatrix * activationMat;
pX = angle(mX);

%% STFT Synthesis
disp('STFT Synthesis')
y = stft_Synthesis(mX,pX,N,soundlength);
%% plot results
subplot(2,1,1)
plot(x)
title('original')
subplot(2,1,2)
plot(y)
title('reconstruction')