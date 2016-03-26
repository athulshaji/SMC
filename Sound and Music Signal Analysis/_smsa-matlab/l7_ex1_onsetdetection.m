clear;clc
%% Onset Detection
% Implement onset detection following the paper by Bello et al.: A Tutorial on Onset Detection in Music Signals.
% 1. Annotate the onsets for the sound file 038 phrase disco simple slow sticks ds.wav, 
% using SonicVisualizer, and save them as a text file.
% 2. Load the sound file 038 phrase disco simple slow sticks ds.wav into Matlab as well as the annotations and 
% plot the waveform.% What is the sampling rate of the sound file? 
% In the same figure, also plot the onset annotation times as green crosses on the y=0 line.
%% Read audio and onset cvs
[x fs] = audioread('l7-038_phrase_disco_simple_slow_sticks_ds.wav');       % read audiofile
onsettime = csvread('l7-onsetgroundtruth_thr0,6300.csv');           % read ground truth onset (sonic visualizer)
%% plot wave and onset on top
t = 0:1/fs:(length(x)-1)/fs;                                               % time vector

% figure;
% plot(t,x);hold on                                                          %
% plot(onsettime,zeros(1,length(onsettime)),'yx')
% title('signal and ground truth onsets')

%% 3. Perform the STFT in Equation (3).
% (a) Choose a window function.
% (b) As a window length in samples choose a power of 2 that corresponds to around 15 ms. 
% Choose the hop size half the size of the window length.
N = round((15/1000)*fs);
window = hanning(N);
% (c) Note that in Equation (3) of the STFT, the sum goes from -N/2 to N/2-1.
% What is the time in s that corresponds to each analyzed frame?

Hfunction = @(x) (x+abs(x))/2;
counter = @(x) (x+1);
W = 1:N;                                                                   % with the high frequency content function.
H = N/2;                                                                   % hop size
pin = 0;
pend = length(x) - N;
Xprev = zeros(1,N);

it = 1;
SD = zeros(1,round(length(x)/H));
weightedenergy = zeros(1,round(length(x)/H));
while pin <= pend
    fftbuffer = x(pin+1:pin+N) .* window;                                  % frame and window the input sound
	X = fft(fftbuffer)';                                                   % compute FFT
    
    % (d) Implement the weighted energy (Equation (4)) with the high frequency content function.
    %     i. Compare the summation range in Equation (4) to the summation range of the fft implementation in Matlab
    %     and arrange the weights of the high frequency content function accordingly.
    weightedenergy(it) = sum(W .* power(abs(X),2) / N);                    % weighted energy
    
    pin = pin + H;
    it = it+1;
end

% figure;plot(weightedenergy)
% title('Weighted energy')

% (e) Substract the maximum of the absolute distance to the mean from the weighted energy,
%     yielding the normalized weighted energy.
normalized_energy = weightedenergy / max(weightedenergy);

% (f) Apply a moving average filter of order L to the normalized weighted energy. 
%     What are the times (in s) corresponding to the low pass filtered result?
L = 3;
averaged_energy = conv(normalized_energy, ones(1,L), 'valid')/3;


% figure;plot(normalized_energy,'b');hold on;plot(averaged_energy,'r')
% legend('normalized_energy','averaged_energy')
% title('Normalized and averaged energy')

% (g) Substract an adaptive threshold (median filter, Equation (21)) from the previous output.
threshold = 0;
scaling_factor = 0.5;
numberofsamples = 10;
thresholded_energy = threshold + scaling_factor * medfilt1(averaged_energy,numberofsamples);



peaks = max(averaged_energy - thresholded_energy,0);
% figure;plot(peaks)
% title('Substract adaptive threshold ot the averaged energy')

peak_threshold = 0.025;
ploc = find(((peaks(1:end-2) < peaks(2:end-1)) & (peaks(2:end-1) < peaks(3:end))) .* (peaks(1:end-2) > peak_threshold));


% from frames to samples
ploc_samples = ploc * H;

figure;
plot(x);hold on
plot(ploc_samples,zeros(size(ploc)),'*')
title('signal')
xlabel('Time (s)')
% (h) Pick local maxima that are above 0 as onsets.
%   Use diff.
%   What are the times (in s) corresponding to the local maxima?
% (i) To debug your code plot all relevant intermediate outputs in the same plot: waveform, weighted energy, 
% low pass filtered energy, median filter, result from thesholding, the
% true onsets and the onsets detected by your algorithm.Do annotated and detected onsets sync?
% (j) Transform the onsets into clicks and listen to a stereo file with the mono song on one channel 
% and the klicks (onsets) on the other channel.
% (k) Try different FIR filter orders L and different ?,? in Equation (21) and choose the parameter 
% configuration that works best.
% (l) Calculate the f-measure on the given .wav file for your best working configuration.










