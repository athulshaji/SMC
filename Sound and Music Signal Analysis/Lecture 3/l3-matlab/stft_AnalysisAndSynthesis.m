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