function [output] = mixatsnr(clean,noise,desired_snr)
% Exercise 1: Make a MATLAB function that is able to mix a desired speech 
% signal with noise at a certain average input SNR level.

desired_snr = db2mag(desired_snr);                  % compute magnitude of desired SNR (db to mag)

desirednoise_var = var(clean) / desired_snr;        % compute variance of attenuated noise 

gain = sqrt(desirednoise_var / var(noise));         % compute attenuation gain

minlen = min(length(clean),length(noise));          % get the length of the shortest signal

output = clean(1:minlen) + gain * noise(1:minlen);  % mix signals

end