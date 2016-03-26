clear; clc
plotYN = 1;
%% Exercise 1: Calculate the Log Magnitude Spectrum
% 1. Download o.wav from Moodle and extract a frame of window length 1024=FFT length, starting at sample 10001.
% 2. Apply a Hanning window on on the extracted frame and plot the time (in s) vs. 
% i) the orginal frame, 
% ii) the Hanning window, scaled by the absolute maximum of the orginal frame, and 
% iii) the frame, windowed by a Hanning window in the same plot, each plotted in a different color.
% 3. Calculate the log10 magnitude spectrum of the windowed frame. 
% Determine the frequenies (in Hz) associated with the first (window length)/2 FFT bins and 
% plot them versus the first (window length)/2 values of the log magnitude spectrum.

% 1. Download o.wav from Moodle and extract a frame of window length 1024=FFT length, starting at sample 10001.
[x fs] = audioread('o.wav');
startingsample = 10001;
N = 1024;
frame = x(startingsample:startingsample+N-1);

% 2. Apply a Hanning window on on the extracted frame and plot the time (in s) vs. 
window = hanning(N);                                                       % create hanning window
windowedframe = frame.*window;                                             % apply window

t = 0:1/fs:(length(frame)-1)/fs;                                           % create time vector

figure                                                                     % and plot the time (in s) vs.
plot(t,frame);hold on                                                      % i) the orginal frame, 
plot(t,window/max(abs(frame)));                                            % ii) the Hanning window, normalized
plot(t,windowedframe)                                                      % iii) the frame, windowed
legend('original frame','hanning','windowed frame')
xlabel('Time (s)')

% 3. Calculate the log10 magnitude spectrum of the windowed frame. 
% Determine the frequenies (in Hz) associated with the first (window length)/2 FFT bins and 
% plot them versus the first (window length)/2 values of the log magnitude spectrum.
magspec = log10(abs(fft(windowedframe)));
f_array = N;

f_vector = linspace(0,fs,N);
plot(f_vector,magspec)
title('magnitude spectrum')
xlabel('frequency (Hz)')

%% Calculate the Spectral Envelope with LPC
% 1. For the windowed frame from the previous exercise, calculate the 41 linear 
% prediction filter coeffients using matlab function lpc with filter order 40.

% 2. Generate the normalized frequency vector of N frequencies in steps of 2*pi/N 
% from 0 to 2*pi - (2*pi)/N and calculate the log10 absolute frequency response of 
% the prediction filter with the coefficients from the lpc and the normalized frequency vector.

% 3. Plot the log10 absolute frequency response versus the log10 magnitude spectrum. Describe the match.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% 1. Calculate 41 LPC filter coefficients.
lpc_order = 41;
[A,E] = lpc(windowedframe,lpc_order);                                      % A: lpc coefficients, E: variance of pred error

% 2. calculate the log10 absolute frequency response of the prediction filter
w = 0:2*pi/N:2*pi - (2*pi)/N;
[H,W] = freqz(1,A,w);
% H: complex frequency response
% W: frequency vector in radians/sample

% 3. Plot the log10 absolute frequency response versus the log10 magnitude spectrum. Describe the match.
diff = max(log10(abs(H))) - max(magspec);
plot(magspec);hold on
plot(log10(abs(H))-diff)
legend('magnitude spectrum','frequency response')

%% 3 Calculate the Spectral Envelope with the Cepstrum
% 3.1. - Take the real part of the inverse discrete Fourier transform of the log FFT magnitude and plot it.
% 3.2. - Get the first 41 values of the previously calculated vector, the cepstrum and plot it.
% 3.3. - Build a vector of length fftSize, with the cepstrum coefficeints in the beginning, 
%        cepstrum coefficients 2 to 41 in inverted order in the end and the rest zeros. 
%        Take the real part of the FFT of that vector and plot it in the same plot as the log Magnitude 
%        spectrum of the zero padded original frame.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 1. 
A = real(ifft(magspec)); % real part of the ifft of the log10 fft magnitude
plot(A);
title('real part of ifft of magnitude spectrum')

% 2. get cepstrum
p = 41;
C = A(1:p);
plot(C);
title('cepstrum');

% 3. vector of cepstrum and zeros
cepstrumpad = zeros(N,1);                                                  % initialize cepstrum padded vector
cepstrumpad(1:p) = C;                                                      % cepstrum coefficients in the beginning,
cepstrumpad(end-p+2:end) = C(p:-1:2);                                      % coeffs 2 to 41 in inverted order in the end

env = real(fft(cepstrumpad));                                              % find the envelope

% plot it
plot(magspec);hold on
plot(env)
legend('magnitude spectrum','envelope')


%% 4. Calculate the Spectral Envelope with the Iterative Cepstrum
% 4.1. - Implement the iterative spectrum method. Intialize A with the log FFT magnitude. 
% Initialize V with -inf. Plot the log magnitude spectrum of the frame, 
% plot the envelope estimate by LPC and by the Cepstrum. 
% 4.2. - In each iteration update
% a) A?max(A,V)
% b) V as the real of the FFT of the rearranged real part of the IFFT
% In each iteration plot both A and V in the same plot of the log magnitude spectrum
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

L = 10;                                                                   % number of iterations
p = 41;

A = magspec;                                                               % initialize A
V = -inf;                                                                  % initializee V

while L>0
    A = max(A,V);
    A2 = ifft(A);
    C = real(A2(1:p));
    
    cepstrumpad = zeros(N,1);                                              % initialize cepstrum padded vector
    cepstrumpad(1:p) = C;                                                  % cepstrum coefficients in the beginning,
    cepstrumpad(end-p+2:end) = C(p:-1:2);                                  % coeffs 2 to 41 in inverted order in the end
    
    V = real(fft(cepstrumpad));                                            % find the envelope
    
    
    % plot it
%     plot(magspec);hold on
%     plot(V)
%     legend('magnitude spectrum','envelope');
%     hold off
%     drawnow;
    
    L = L - 1;
end

%% 5. Implement the MFCC
% 1. - Use the same frame of o.wav as above.
% 2. - Run the pre-emphasis filter:
%           y[n] = x[n] - .97x[n-1] (1)
% over the frame. Plot it in the same figure as the original. What has changed?
% 3. - Apply a Hanning window on the pre-emphasized frame.
% 4. - Take the magnitude spectrum of the previous output.
% 5. - (Frequency spacing) Generate 42 frequencies according to the following rule: 
%   The lowest frequency should be 133.3333 Hz, the next 12 frequencies should 
%   be linearly spaced, i.e. have a frequency distance to the previous frequency 
%   of 66.66666666 Hz. From the 13th frequency on the spacing should be logarithmically, 
%   i.e. the next frequency should be derived from the previous frequency by multiplication 
%   with 1.0711703. Plot the Frequency points and comment on it.
% 6. - Now create a 40 x 512 matrix with a triangle filter in each row using 
%   the previously generated 42 frequencies:
%   (a) Each row should have 0s everywhere except between freq(i) until freq(i+2).
%   (b) At freq(i) the values should linearly increase from 0 up to 2/(freq(i+2)-freq(i)) 
%   at freq(i+1) then linearly decrease to 0 at position freq(i+2)
%   (c) Plot several rows into the same plot and comment on the graph.
% 7. - Multiply the triangle filter matrix with the magnitude vector and plot the result.
% 8. - Take the log10 of the previous and plot it.
% 9. - Perform the discrete cosine transform the previous output e (N=40).
% 10. - Plot the the first 13 values of the previous output, the Mel Frequency Cepstrum Coefficients.

fltrdframe = filter([1 -0.97],1,frame);                                    % apply pre-emphasis filter 
window = hanning(N);                                                       % create hanning window
windowedframe = fltrdframe.*window;                                        % apply window
magspec = (abs(fft(windowedframe)));                                       % magnitude spectrum

% Frequency spacing
% linear frequencies
f0 = 133.3333;                                                             % lowest frequency in Hz
dist = 66.66666666;                                                        % distance between frequencies
logratio = 1.0711703;                                                      % logarithmical distance ratio
nfreqs = 42;                                                               % number of frequencies in linear distance
nfreqslinear = 12;

freqs = [f0, zeros(1,nfreqs-1)];
for i = 2:nfreqs
    if i <= nfreqslinear
        freqs(i) = freqs(i-1)+dist;
    else
        freqs(i) = freqs(i-1)*logratio;
    end
end

% fftfreqs = 0:fs/N:N/2-fs/N;
ntr = 40;       % number of triangles
N2 = 512;        % fft size
fndx = 1;       % frequency index
trianglemat = zeros(ntr,N2);
bins = 1:512;
for i = 1:ntr
    trfreqs = freqs(fndx:fndx+2);                                          % frequencies for this triangle
    tramp = 2/(trfreqs(3)-trfreqs(1));                                     % amplitude of triangle
    
    trvalues = interp1(trfreqs,[0,tramp,0],bins*fs/N);
    trvalues(isnan(trvalues))=0;
    trianglemat(i,:) = trvalues;
    
    fndx = fndx + 1;
    
    if plotYN
        plot(trvalues);
        hold on
    end
end






