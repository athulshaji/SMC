clear;clc

sndtoload        = 'mix_09_10.flac';
inputdir         = '../testingset/';
dictionarytoload = 'atomDictionary.mat';

% load dictionary and create activation matrix
load(dictionarytoload)
% load sound
[x,fs] = audioread([inputdir,sndtoload]);x = mean(x(1:sgmtlen*fs,:),2);
[x,fs] = myresample(x,fs,8000);

nsndsindict = size(atomDictionary,2)/K;

%% fourier analysis of sndtoload
window  = sqrt(hanning(N));
M       = length(window);                                       % analysis window size
H       = floor(N/2);                                           % hop size
N2      = N/2+1;                                                % size of positive spectrum, includes sample 0
sndlen  = length(x);                                            % compute input sound length
pin     = 1;                                                    % initialize sound pointer in
pend    = sndlen-N;                                             % initialize sound pointer end
mX      = [];
pX      = [];

% %-----analysis-----%
while pin <= pend
    fftbuffer = x(pin:pin+M-1) .* window;                       % frame and window the input sound
    X = fft(fftbuffer);                                         % compute FFT
    mX = [mX abs(X(1:N2))];                                     % magnitude spectrum (linear)
    pX = [pX angle(X(1:N2))];                                   % phase spectrum
    pin = pin+H;                                                % advance sound pointer
end

%% transformations: Non-negative Matrix Factorization
atomMatrix = atomDictionary;                                    % initialize atom matrix
activationMat = rand(size(atomMatrix,2),size(mX,2));            % initialize activation matrix
BB = ones(size(mX));
for it = 0:iterations
    activationMat = activationMat .* (atomMatrix' * (mX./(atomMatrix * activationMat+eps)))./(atomMatrix' * BB);
end

%%
for i = 1:size(atomMatrix,2)
    mX = atomMatrix(:,i) * activationMat(i,:);
    % pX = angle(mX);
    
    %-----synthesis-----%
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
    xhat(:,i) = y;
end


%% create estimated signals
se = zeros(size(filesused_snd{1},1),size(xhat,2)/K);
for i = 1:size(xhat,2)/K
    se(:,i) = sum(xhat(:,(50*(i-1)+1:50*(i-1)+50)),2);
    se(:,i) = se(:,i) / max(se(:,i));    
end

s = cell2mat(filesused_snd);


[SDR,SIR,SAR,perm]=bss_eval_sources(se',s')

%%
% figure;
% for src=1:nsndsindict
% clear sound
% plot(se(:,src))
% soundsc(se(:,src),fs)
% pause
% end
% 
% 
% for src=1:nsndsindict
% audiowrite(sprintf('mix_46_22_2-_se%i.wav',src),se(:,src)/max(se(:,src)),fs)
% end
%%
% %% plot stuff
% % PLOT BASIS Vectors
% freq = linspace(0, fs/2, N/2+1);
% time = linspace(0, length(x)/fs, size(mX,2)); 
% 
% figure;
% for i=50:100%sum(K)
%     plot(time, (i-1)*max(max(activationMat))+(activationMat(i,:)),'LineWidth', 3)
%     hold on
% end
% ylabel('Activations')
% xlabel('Time (seconds)')
% 
% set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);