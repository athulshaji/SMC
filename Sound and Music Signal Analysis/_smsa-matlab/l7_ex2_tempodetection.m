%% Tempo estimation
% Build a BpM estimator.
% i. Apply it to 038 phrase disco simple slow sticks ds.wav.
% ii. Calculate the local energy (Bello et al. Equation (4) with Wk = 1?k).
% iii. Calculate the autocorrelation of the local energy.
% iv. Compare the results with your tapping results.

%% Read audio
[x fs] = audioread('l7-038_phrase_disco_simple_slow_sticks_ds.wav');       % read audiofile

%% 3. Perform the STFT in Equation (3).
N = round((15/1000)*fs);
window = hanning(N);

H = N/2;                                                                   % hop size
pin = 0;
pend = length(x) - N;

it = 1;
localenergy = zeros(1,round(length(x)/H));

while pin <= pend
    fftbuffer = x(pin+1:pin+N) .* window;                                  % frame and window the input sound
	X = fft(fftbuffer)';                                                   % compute FFT
    
    % (d) Implement the weighted energy (Equation (4)) with the high frequency content function.
    %     i. Compare the summation range in Equation (4) to the summation range of the fft implementation in Matlab
    %     and arrange the weights of the high frequency content function accordingly.
    localenergy(it) = sum(power(abs(X),2) / N);                       % weighted energy
    
    pin = pin + H;
    it = it+1;
end

[acf,lags] = xcorr(localenergy); % lags: time intervals, for each time interval, you have a autocorrelation value
plot(lags,acf)
% units of lags is block size divided by sampling rate
% multiply block size by hop size to get frequency, then get period and
% then get bpm