% Exercise 1: 
% Write MATLAB functions that can compute the STFT and iSTFT of a signal 
% using a certain block length, 50% overlapping blocks and a Hanning window.

function y = wienerfilterfreqdomain(x,N)
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
    
    % noise estimation variables:
    noiseMag = zeros(N/2+1,1);                                             % initialize estimated noise magnitude
    prevmX = zeros(N/2+1,1);                                               % initialize previous input magnitude spectrum
    alpha_a = 0.1;                                                         % noise estimator attack coefficient
    alpha_d = 0.1;                                                         % noise estimator decay coefficient
    beta_a = 0.1;                                                          % smooth speech attack coefficient
    beta_d = 0.1;                                                          % smooth speech decay coefficient
    
    while pin <= pend
        %-----analysis-----%
        fftbuffer = zeros(N,1);                                            % initialize buffer for fft
        fftbuffer(1:M) = x(pin+1:pin+M) .* window;                         % frame and window the input sound
        X = fft(fftbuffer);                                                % compute FFT
        mX = 20*log10(abs(X(1:N2)));                                       % magnitude spectrum
        pX = angle(X(1:N2));                                               % phase spectrum        
        
        %-----Noise magnitude estimator-----%
        
        % to reduce its temporal fluctuation, the magnitude of the noisy 
        % speech spectrum is smoothed:
        mX = (mX >= prevmX) .* (beta_a * prevmX + (1-beta_a) * mX) + ...
            (mX < prevmX) .* (beta_d * prevmX + (1-beta_d) * mX);
        
        % estimate magnitude of the noise signal:
        noiseMag = (mX >= noiseMag) .* (alpha_a * noiseMag + (1-alpha_a) * mX) + ...
            (mX < noiseMag) .* (alpha_d * noiseMag + (1-alpha_d) * mX);
        
        prevmX = mX;
        %-----Compute filter-----%
        Hw = 1 - periodogram(noiseMag)/periodogram(mX);
%         Hw = 1 - power(noiseMag,2)/power(mX,2);
        
        %-----transformations-----%
        mY = Hw * mX;                                                      % filter input magnitude
%         mY = mX;
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