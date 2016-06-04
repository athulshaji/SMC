function y = stft_NNMF(x,N,ploton,K)
% Analysis and Synthesis of a sound using the short-time fourier transform
% x: input sound,
% N: FFT size
% y: output sound
% window: window function
window  = sqrt(hanning(N));
M       = length(window);                                       % analysis window size
H       = floor(N/2);                                           % hop size
N2      = N/2+1;                                                % size of positive spectrum, includes sample 0
sndlen  = length(x);                                            % compute input sound length
pin     = 1;                                                    % initialize sound pointer in
pend    = sndlen-N;                                             % initialize sound pointer end
mX      = [];
pX      = [];

%% %-----analysis-----%
while pin <= pend
    fftbuffer = x(pin:pin+M-1) .* window;                       % frame and window the input sound
    X = fft(fftbuffer);                                         % compute FFT
    mX = [mX abs(X(1:N2))];                                     % magnitude spectrum (linear)
    pX = [pX angle(X(1:N2))];                                   % phase spectrum        
    pin = pin+H;                                                % advance sound pointer
end

%% transformations: Non-negative Matrix Factorization
atomMatrix = abs(rand(N2,K));                                              % initialize atom matrix
activationMat = rand(K,size(mX,2));                                        % initialize activation matrix
BB = ones(size(mX));
for it = 0:iterations
    atomMatrix = atomMatrix .* (mX ./ (atomMatrix * activationMat+eps) * activationMat')./(BB  * activationMat');
    activationMat = activationMat .* (atomMatrix' * (mX./(atomMatrix * activationMat+eps)))./(atomMatrix' * BB + zeros(size(atomMatrix'*BB)));
    
    if ploton
        imagesc(atomMatrix);title('atom Matrix')
        drawnow
    end
end

mX = atomMatrix * activationMat;
pX = pX;% angle(mX);

%% %-----synthesis-----%  
pin = 1;                                                        % initialize sound pointer in
y = zeros(sndlen,1);                                            % initialize output sound
for idx = 1:size(mX,2)
    mY = mX(:,idx);                                             % get one frame of magnitude spectrum
    pY = pX(:,idx);                                             % get one frame of phase spectrum
    
    Y           = zeros(N,1);                                   % initialize output spectrum
    Y(1:N2)     = mY .* exp(1i.*pY);                            % generate positive frequencies (linear)
    Y(N2+1:N)   = mY(N2-1:-1:2) .* exp(-1i.*pY(N2-1:-1:2));     % generate negative frequencies (linear)
        
    fftbuffer       = real(ifft(Y)) .* window;                  % compute inverse FFT
    y(pin:pin+M-1)  = y(pin:pin+M-1) + fftbuffer;               % overlap-add
    pin             = pin + H;                                  % advance sound pointer
end
end