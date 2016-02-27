% build a simple analysis/synthesis system in MATLAB that can be used 
% for modification and transformation of audio signals. 
% First, build a script in MATLAB that can read an audio file (of your choice), 
% extract segments of length N (corresponding to, say, 30 ms) with an overlap of 50%. 
% Then, for each segment, apply a square-root hanning window and compute the 
% Fourier transform (use fft()). The output of the Fourier transform can now 
% be modified, if pleased, but for now we leave it alone. 

% Next, your script should take the (un-)modified coefficients of the Fourier 
% transform and compute the inverse Fourier transform of them (use ifft()). 

% Finally, apply a square-root hanning window to the output and perform 
% overlap add and move on to the next segment. Once all segments have been 
% processed, compare the input and output of the system and verify that the 
% system has perfect reconstruction (up to the numerical noise floor). 

% Finally, implement some kind of modification in the Fourier domain, 
% like a low-pass filter that removes everything above fs/4.

clearvars -except x fs;clc;close all
if ~exist('x','var')
    [x,fs] = audioread('../../data/twoMaleTwoFemale20Seconds.wav');        % desired signal
end
%%
segmentlength = 0.03;                                                      % 30 ms
N = segmentlength * fs;
window = sqrt(hanning(N-1));

y = stft_AnalysisAndSynthesis(x,N,window);

%%
px1 = subplot(2,1,1);
p1 = plot(x);
px2 = subplot(2,1,2);
p2 = plot(y);
linkaxes([px1,px2],'xy')
