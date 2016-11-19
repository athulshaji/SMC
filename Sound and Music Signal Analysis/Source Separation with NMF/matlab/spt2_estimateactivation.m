clear;clc

sndtoload        = 'mix_08_37.flac';
fidx             = 1;
inputdir         = '../testingset/';
dictionarytoload = 'atomDictionary.mat';

% load dictionary and create activation matrix
load(dictionarytoload)
% load sound
[x,fs] = audioread([inputdir,sndtoload]);x = mean(x(1:sgmtlen*fs,:),2);
[input,fs] = myresample(x,fs,8000);

[v,fs] = audioread('../babble30Seconds.wav');v = mean(v(1:10*fs,:),2);
[v,fs] = myresample(v,fs,8000);

desired_snr_array = flip([100,30,25,20,15,10,5,0,-10,-50,-75,-100]);                    % define desired SNR

nsndsindict = size(atomDictionary,2)/K;

SDR     = zeros(nsndsindict,length(desired_snr_array));
SIR     = zeros(nsndsindict,length(desired_snr_array));
SAR     = zeros(nsndsindict,length(desired_snr_array));
perm    = zeros(nsndsindict,length(desired_snr_array));



for snr_i = 1:length(desired_snr_array)
    desired_snr = desired_snr_array(snr_i)

    x = l1_ex1_mixatsnr(input,v,desired_snr);                                      % mix signals at desired SNR
    x = x/max(x);
%     audiowrite(sprintf('../outsnd/snr_%04i_mix.wav',desired_snr),x,fs);
    
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
%     freq = linspace(0, fs/2, N/2+1);
%     time = linspace(0, length(mX)/fs, size(mX,2)); 
%     for i=50:100%sum(K)
%         subplot(1,2,1);title('atoms');
%         plot((i-1)*max(max(atomMatrix))+(1-atomMatrix(:,i)),freq,'LineWidth', 3)
%         ylabel('frequency')
%         hold on
%         subplot(1,2,2);title('activations');
%         plot(time, (i-1)*max(max(activationMat))+(activationMat(i,:)),'LineWidth', 3)
%         xlabel('time')
%         hold on
%     end
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
    % plot(x)
    % figure;plot(y)
    % soundsc(x,fs);
    snd_idx = 1;
    % soundsc(sum(xhat(:,(snd_idx:snd_idx+50)),2),fs)

    %% create estimated signals
    se = zeros(size(filesused_snd{1},1),size(xhat,2)/K);
    for i = 1:size(xhat,2)/K
        se(:,i) = sum(xhat(:,(50*(i-1)+1:50*(i-1)+50)),2);
        se(:,i) = se(:,i) / max(se(:,i));
        max(se(:,i))
%         audiowrite(sprintf('../outsnd/snr_%04i_s_%i.wav',desired_snr,i),se(:,i),fs);
    end

    s = cell2mat(filesused_snd);

    [SDR_tmp,SIR_tmp,SAR_tmp,perm_tmp] = bss_eval_sources(se',s');
    SDR(:,snr_i) = SDR_tmp;
    SIR(:,snr_i) = SIR_tmp;
    SAR(:,snr_i) = SAR_tmp;
%     perm_tmp

end

%%
for src = 1:nsndsindict
    figure;
    subplot(3,1,1)
    plot(SDR(src,:));title(sprintf('SDR signal %i',src));
    ax = gca;ax.XTickLabel = num2str(desired_snr_array');

    subplot(3,1,2)
    plot(SIR(src,:));title(sprintf('SIR signal %i',src));
    ax = gca;ax.XTickLabel = num2str(desired_snr_array');

    subplot(3,1,3)
    plot(SAR(src,:));title(sprintf('SAR signal %i',src));xlabel('SNR')
    ax = gca;ax.XTickLabel = num2str(desired_snr_array');
end

%%
figure;
for src=1:nsndsindict
clear sound
plot(se(:,src))
soundsc(se(:,src),fs)
pause
end
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