function y = stft_AnalysisAndSynthesis(x,N,window)
    % Analysis of a sound using the short-time fourier transform
    % x: input sound,
    % N: FFT size
    % y: output sound
    % window: window function
    
    H = floor(N/2);                                                        % hop size
    N2 = N/2+1;                                                            % size of positive spectrum, includes sample 0
    soundlength = length(x);                                               % compute input sound length
    pin = 1;                                                               % initialize sound pointer in
    pend = soundlength-N;                                                  % initialize sound pointer end
    y = zeros(soundlength,1);                                              % initialize output sound
    
    while pin <= pend
        %-----analysis-----%
        fftbuffer = x(pin:pin+N-1) .* window;                         % frame and window the input sound
        X = fft(fftbuffer);                                                % compute FFT
        mX = 20*log10(abs(X(1:N2)));                                       % magnitude spectrum
        pX = angle(X(1:N2));                                               % phase spectrum        
        
        %-----transformations-----%
        mY = mX;
        pY = pX;
        
        %-----synthesis-----%
        Y = zeros(N,1);                                                    % initialize output spectrum
        Y(1:N2) = 10.^(mY/20) .* exp(1i.*pY);                              % generate positive frequencies
        Y(N2+1:N) = 10.^(mY(N2-1:-1:2)/20) .* exp(-1i.*pY(N2-1:-1:2));     % generate negative frequencies
        fftbuffer = real(ifft(Y));                                         % compute inverse FFT
        y(pin:pin+N-1) = y(pin:pin+N-1) + window .* fftbuffer;             % overlap-add to generate output sound
        
        pin = pin+H;                                                       % advance sound pointer
    end
end