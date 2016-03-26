clear,clc
plotYN = 0;
%% Sound and Music Signal Analysis 9 (Dynamic Time Warping)
% Read in BeethovenPiano.wav and BeethovenOrchestra.wav? from the subfolder data in the TSM Toolbox root directory.
addpath(genpath('MATLAB_TSM-Toolbox_1.0'))
[x1 fs] = audioread('MATLAB_TSM-Toolbox_1.0/data/BeethovenPiano.wav');
[x2 fs] = audioread('MATLAB_TSM-Toolbox_1.0/data/BeethovenOrchestra.wav');

x1 = x1(1:min(length(x1),length(x2)));
x2 = x2(1:min(length(x1),length(x2)));
%% Generate a stereo sound file of the duration of the shorter of the two signals with each signal 
% (piano and orchestra) on a separate stereo channel.
% Listen to the stereo sound file. What is wrong with the example?
stereofile = [x1;x2];
% stereofile = zeros(2,min(length(x1),length(x2)));
% stereofile(1,:) = x1;
% stereofile(2,:) = x2;

%% 2. Calculate the optimal time warping path q for the alignment of the orchestra signal to the piano signal.
% (a) On each of the original signals (piano and orchestra) perform block processing with 
% window size w=1024 (fft size), hop size h=512, and a Hamming window, using FFT.
% i. How many frames do you have?
% ii. For each file plot the magnitude spectrum for each block and the first 100 FFT
% bins, using the Matlab imagesc.
N = 1024;
H = 512;
window = hamming(N);
framesnumber = round(length(x1)/H);
mag1 = zeros(framesnumber,N);
mag2 = zeros(framesnumber,N);
it = 1;

pin = 0;                                                                   % initialie fft buffer index
pend = length(x1) - N;                                                     % 
while pin <= pend
    fftbuffer1 = x1(pin+1:pin+N) .* window;                                % frame and window the input sound
	X1 = fft(fftbuffer1);                                                  % compute FFT
    mag1(it,:) = abs(X1);                                                  % get magnitude spectrum
    
    fftbuffer2 = x1(pin+1:pin+N) .* window;                                % frame and window the input sound
	X2 = fft(fftbuffer2);                                                  % compute FFT
    mag2(it,:) = abs(X2);                                                  % get magnitude spectrum
    
    pin = pin + H;                                                         % update fft buffer index
    it = it + 1;                                                           % update iteration
    
	% plot first 100 bins of magnitude spectrum
    if plotYN    
        subplot(2,1,1)
        imagesc(mag1(1:100))
        title('magnitude spectrum of signal 1')
        subplot(2,1,2)
        imagesc(mag2(1:100))
        title('magnitude spectrum of signal 2')
        drawnow
    end
end

% (b) For each combination of a frame of the piano signal and a frame of the orchestra signal, 
% calculate the cosine dissimilarity, Eq. 3.14 (MM, 2015), resulting in the cost matrix Eq. 3.13 (MM, 2015). 
% To do so, use the magnitude spectrum for the first 100 bins in each frame. 
% You can use the function simmx.m that calculates the similarity between frames, 
% the second term on the right side of Eq. 3.14 (MM, 2015). 
% simmx.m can be downloaded from Moodle. As first input argument to simmx use the piano signal, 
% as the second argument the orchestra signal. 
% Plot the similarity matrix, use colormap with parameter 1-gray.
costmat = zeros(framesnumber);
for n = 1:framesnumber
    for m = 1:framesnumber
        % use the magnitude spectrum for the first 100 bins in each frame
        m1 = mag1(n,1:100);
        m2 = mag2(m,1:100);
        % compute cosine dissimilarity
        cosdis = 1 - simmx(m1',m2');
        % add cosine dissimilarity to cost matrix
        costmat(n,m) = cosdis;
    end
end
if plotYN
    imagesc(costmat);
    colormap gray
    hold on
end
% (c) Calculate the optimal warping path according to Eq. 3.22-3.25, 3.27-3.29 (MM, 2015). 
% To do so, download dpfast.m from the Moodle and use it according to [n,m,D]=dpfast(C), 
% where the first input argument C is the dissimilarity or cost matrix, Eq. 3.13 (MM 2015), 
% n and m are the indices of the optimal warping path q calculated through Eq. 3.27-3.29 (MM, 2015), 
% and D is the cumulated cost matrix, calculated through Eq. 3.23-3.25 (MM, 2015). 
% On top of a plot of the similarity matrix C plot the optimal warping path q.

[n,m,D] = dpfast(costmat);

if plotYN
    
end

% (d) Plot the cumulative cost matrix D. What is the total accumulated cost of the optimal warping path?