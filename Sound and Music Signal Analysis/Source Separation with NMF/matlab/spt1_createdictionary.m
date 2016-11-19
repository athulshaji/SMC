clear;clc

trnpath     = '../trainingset/';

ploton      = 0;                                                           % plot or not the NNMF
K           = 50;                                                          % number of atoms
iterations  = 100;                                                         % iterations to find atoms and activation matrices
atomlen     = 0.025;                                                       % 30 ms
sgmtlen     = 10;                                                          % get only the first sgmtlen seconds of the input sound

% read files in trnpath
trnroot         = dir(trnpath);
trnfolders      = trnroot(4:end);
atomDictionary  = [];
filesused_name  = {}; 
filesused_snd   = {};
for f = 1:length(trnfolders)
    fullpath     = [trnpath,trnfolders(f).name,'/'];
    filesinfolder= dir([fullpath,'*.flac']);
    disp(fullpath)
    for ff = 1:length(filesinfolder)
        sndfilename = [fullpath,filesinfolder(ff).name];
        disp(sndfilename)
        % load and resample sound
        [x,fs] = audioread(sndfilename);x = mean(x(1:sgmtlen*fs,:),2);       % desired signal
        [x,fs] = myresample(x,fs,8000);
        N = round(atomlen * fs/2) * 2;
        [y,atomMatrix,activationMat,mX,pX] = stft_NNMF(x,N,ploton,K,iterations);
        
        atomDictionary = [atomDictionary,atomMatrix];
        filesused_name{end+1} = sndfilename;
        filesused_snd{end+1}  = x;
    end 
end

save('atomDictionary.mat','atomDictionary','filesused_name','filesused_snd','K','iterations','sgmtlen','N','atomlen');