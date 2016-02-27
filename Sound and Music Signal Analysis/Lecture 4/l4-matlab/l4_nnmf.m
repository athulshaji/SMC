clearvars -except x fs;clc;close all
addpath('../../_matlabscripts')
disp(smsa_library)

if ~exist('x','var')
    [x,fs] = audioread('../../data/twoMaleTwoFemale20Seconds.wav');        % desired signal
end

segmentlength = 0.025;                                                     % 30 ms
N = segmentlength * fs;                                                    % fft size
window = sqrt(hanning(N-1));                                               %
soundlength = length(x);
K = 50;                                                                    % number of atoms
iterations = 100;

atomMatrix = rand(N/2+1,K);
activationMAt = rand(K,1);

[mX, pX] = stft_Analysis(x,N,window);
y = stft_Synthesis(mX,pX,N,soundlength);