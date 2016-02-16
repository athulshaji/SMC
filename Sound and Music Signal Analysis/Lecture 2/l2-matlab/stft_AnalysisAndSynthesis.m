% Exercise 1: 
% Write MATLAB functions that can compute the STFT and iSTFT of a signal 
% using a certain block length, 50% overlapping blocks and a Hanning window.

function y = stft(x,N)
    % Analysis of a sound using the short-time fourier transform
    % x: input sound,
    % N: FFT size
    % y: output sound
    
    M = floor(N/2)-1;                                                      % window length
    H = floor(M/2);                                                        % hop size
    N2 = N/2+1;                                                            % half of the FFT size
    soundlength = length(x);                                               % compute input sound length
    window = hanning(M);window = window/sum(window);                       % initialize window
    pin = 1;                                                               % initialize sound pointer in
    pend = soundlength-M;                                                  % initialize sound pointer end
    y = zeros(soundlength,1);                                              % initialize output sound
    
    while pin <= pend
        %-----analysis-----%
        fftbuffer = zeros(N,1);                                            % initialize buffer for fft
        fftbuffer(1:M) = x(pin+1:pin+M) .* window;                         % frame and window the input sound
        X = fft(fftbuffer);                                                % compute FFT
        mX = 20*log10(abs(X(1:N2)));                                       % magnitude spectrum
        pX = angle(X(1:N2));                                               % phase spectrum        
        
        %-----transformations-----%
        mY = mX;
        pY = pX;
        
        %-----synthesis-----%
        Y = zeros(N,1);                                                    % initialize buffer for ifft
        Y(1:N2) = 10.^(mY/20) .* exp(1i.*pY);                              % generate positive frequencies
        Y(N2+1:N) = 10.^(mY(N2-1:-1:2)/20) .* exp(-1i.*pY(N2-1:-1:2));     % generate negative frequencies
        fftbuffer = real(ifft(Y));                                         % compute inverse FFT
        y(pin:pin+M-1) = y(pin:pin+M-1) + H * fftbuffer(1:M);              % overlap-add to generate output sound
        
        pin = pin+H;                                                       % advance sound pointer
    end
end