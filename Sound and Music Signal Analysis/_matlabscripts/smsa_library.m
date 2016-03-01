function message = smsa_library
  assignin('base','stft_AnalysisAndSynthesis',@stft_AnalysisAndSynthesis);
  assignin('base','stft_Analysis',@stft_Analysis);
  assignin('base','stft_Synthesis',@stft_Synthesis);
  assignin('base','mixatsnr',@mixatsnr);
  assignin('base','createVandermondeMat',@createVandermondeMat);
  assignin('base','vandermonde',@vandermonde);
  message='Done importing functions to workspace';
end

function y = stft_AnalysisAndSynthesis(x,N,window)
    % Analysis of a sound using the short-time fourier transform
    % x: input sound,
    % N: FFT size
    % y: output sound
    % window: window function
    M = length(window);                                                    % analysis window size
    H = floor(N/2);                                                        % hop size
    N2 = N/2+1;                                                            % size of positive spectrum, includes sample 0
    soundlength = length(x);                                               % compute input sound length
    pin = 1;                                                               % initialize sound pointer in
    pend = soundlength-N;                                                  % initialize sound pointer end
    y = zeros(soundlength,1);                                              % initialize output sound
    
    while pin <= pend
        %-----analysis-----%
        fftbuffer = x(pin:pin+M-1) .* window;                              % frame and window the input sound
        X = fft(fftbuffer);                                                % compute FFT
        mX = abs(X(1:N2));                                                 % magnitude spectrum (linear)
        pX = angle(X(1:N2));                                               % phase spectrum        
        % mX = 20*log10(abs(X(1:N2)));                                     % % magnitude spectrum (dB)

        %-----transformations: LOW PASS FILTER-----%
        mY = zeros(N2,1);                                                  % initialize output magnitude spectrum
        lpf_idx = 1:ceil(N2/2);                                            % low pass filter indeces (rest will be zeros)
        mY(lpf_idx) = mY(lpf_idx) + mX(lpf_idx);                           % use only samples in lpf indices
        % mY = mX;                                                         % no transformation
        pY = pX;                                                           % output phase spectrum
        
        %-----synthesis-----%        
        Y = zeros(N,1);                                                    % initialize output spectrum 
        Y(1:N2) = mY .* exp(1i.*pY);                                       % generate positive frequencies (linear)
        Y(N2+1:N) = mY(N2-1:-1:2) .* exp(-1i.*pY(N2-1:-1:2));              % generate negative frequencies (linear)
        % Y(1:N2) = 10.^(mY/20) .* exp(1i.*pY);                            % generate positive frequencies (dB)
        % Y(N2+1:N) = 10.^(mY(N2-1:-1:2)/20) .* exp(-1i.*pY(N2-1:-1:2));   % generate negative frequencies (dB)
        fftbuffer = real(ifft(Y));                                         % compute inverse FFT
        
        y(pin:pin+M) = y(pin:pin+M) + fftbuffer;                           % overlap-add
        
        pin = pin+H;                                                       % advance sound pointer
    end
end

% px1 = subplot(2,1,1);
% p1 = plot(abs(X));
% px2 = subplot(2,1,2);
% p2 = plot(abs(Y));
% linkaxes([px1,px2],'xy')


function [mX, pX] = stft_Analysis(x,N,window)
    % Analysis of a sound using the short-time fourier transform
    % x: input sound,
    % N: FFT size
    % y: output sound
    % window: window function
    M = length(window);                                                    % analysis window size
    H = floor(N/2);                                                        % hop size
    N2 = N/2+1;                                                            % size of positive spectrum, includes sample 0
    soundlength = length(x);                                               % compute input sound length
    pin = 1;                                                               % initialize sound pointer in
    pend = soundlength-N;                                                  % initialize sound pointer end
    mX = [];
    pX = [];
    
    while pin <= pend
        %-----analysis-----%
        fftbuffer = x(pin:pin+M-1) .* window;                              % frame and window the input sound
        X = fft(fftbuffer);                                                % compute FFT
        mX = [mX abs(X(1:N2))];                                            % magnitude spectrum (linear)
        pX = [pX angle(X(1:N2))];                                          % phase spectrum        
        % mX = 20*log10(abs(X(1:N2)));                                     % % magnitude spectrum (dB)
        
        pin = pin+H;                                                       % advance sound pointer
    end
end


function y = stft_Synthesis(mX,pX,N,soundlength)
    % Analysis of a sound using the short-time fourier transform
    % mX: input magnitude spectrum matrix
    % pX: input phase spectrum matrix
    % N: FFT size
    % soundlength: length of output sound in samples
    % y: output sound
    M = N-1;
    N2 = N/2+1;                                                            % size of positive spectrum, includes sample 0
    H = floor(N/2);                                                        % hop size
    pin = 1;                                                               % initialize sound pointer in
    pend = size(mX,2);                                                     % initialize sound pointer end
    y = zeros(soundlength,1);                                              % initialize output sound
    idx = 1;                                                               % initialize column index
    
    while idx <= pend
        mY = mX(:,idx);                                                    % get one frame of magnitude spectrum
        pY = pX(:,idx);                                                    % get one frame of phase spectrum
        %-----synthesis-----%        
        Y = zeros(N,1);                                                    % initialize output spectrum 
        Y(1:N2) = mY .* exp(1i.*pY);                                       % generate positive frequencies (linear)
        Y(N2+1:N) = mY(N2-1:-1:2) .* exp(-1i.*pY(N2-1:-1:2));              % generate negative frequencies (linear)
        % Y(1:N2) = 10.^(mY/20) .* exp(1i.*pY);                            % generate positive frequencies (dB)
        % Y(N2+1:N) = 10.^(mY(N2-1:-1:2)/20) .* exp(-1i.*pY(N2-1:-1:2));   % generate negative frequencies (dB)
        fftbuffer = real(ifft(Y));                                         % compute inverse FFT
        y(pin:pin+M) = y(pin:pin+M) + fftbuffer;                           % overlap-add
        pin = pin + H;                                                     % advance sound pointer
        idx = idx + 1;                                                     % go to next column in magnitude matrix
    end
end

function [output] = mixatsnr(clean,noise,desired_snr)
% Exercise 1: Make a MATLAB function that is able to mix a desired speech 
% signal with noise at a certain average input SNR level.

desired_snr = db2mag(desired_snr);                  % compute magnitude of desired SNR (db to mag)

desirednoise_var = var(clean) / desired_snr;        % compute variance of attenuated noise 

gain = sqrt(desirednoise_var / var(noise));         % compute attenuation gain

minlen = min(length(clean),length(noise));          % get the length of the shortest signal

output = clean(1:minlen) + gain * noise(1:minlen);  % mix signals

end

function [Z] = createVandermondeMat(w,L,M)
% input:
%   w: fundamental frequency
%   L:number of harmonics
%   M: number of samples
% output:
%   Z: MxL vandermonde matrix
    
    z = exp(1i*w * (1:L));          % compute harmonics
    Z = zeros(M-1,L);               % initialize Vandermonde matrix
    for n = 1:M                     % for each sample
        Z(n,:) = z.^(n-1);          % fill Vandermonde matrix
    end   
end