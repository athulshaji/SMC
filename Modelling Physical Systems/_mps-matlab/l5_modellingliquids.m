%% Modelling Physical Systems: liquids
%% Single Bubble
% 1. Change values of the radius and verify that it affects the frequency of the bubble
% 2. Change the values of epsilon (slope) and verify that it corresponds to
% the speed of change of radius.
% 3. Try with the same parameteres from the paper
%           r =  [10,7,4,1,0.5,0.3] in mm 
%           and epsilon = [0,0.05, 0.1]
% 4. what is missing to have a perceptually more realistic bubble?
% The attack. The impact is really difficult to model.

clear; clc
fs = 44100;                                                                % sampling frequency
dt = 1/fs;                                                                 % time step
duration = 2;                                                              % duration in time
t = 0:dt:duration;                                                         % time vector 

a = 1;                                                                     % amplitude
r = 4/1000;                                                                % radius of the bubble 0.00015 < r < 0.15
f0 = 3/r;                                                                  % fundamental frequency 
epsilon = 0.1;                                                             % rise of the frequency

d = 0.043 * f0 + 0.0014 * f0^(-3/2);                                       % calculate damping factor

freqrise = epsilon * d;                                                    % compute frequency rise factor

f = f0 * (1+freqrise*t);                                                   % compute frequency change
i = a * sin(2*pi*f.*t).*exp(-d*t);                                         % compute impulse response of bubble

soundsc(i,fs)                                                              % play sound
%% Simulate the sound of a dense stream of bubbles
clear,clc

% General parameters: 
duration = 2;                                                              % duration in time
fs = 44100;                                                                % sampling frequency
dt = 1/fs;                                                                 % time step
t = 0:dt:duration;                                                         % time vector 
a = 1;                                                                     % amplitude of the bubble

% parameters of the bubbles:
N = 50;                                                                    % number of bubbles
rmin = 0.2;                                                                % minimum bubble radius in mm
rmax = 50;                                                                 % maximum bubble radius in mm

% computations
r = ((rmax-rmin).*rand(N,1) + rmin)/1000;                                  % N radius of random length between rmin and rmax
f0 = (3./r);                                                               % compute the frequencies of each bubble

epsilon = rand(N,1)*0.5;                                                   % rise of the frequency
d = 0.043 .* f0 + 0.0014 .* f0.^(-3/2);                                    % calculate damping factor
freqrise = epsilon .* d;                                                   % compute frequency rise factor

f = repmat(f0,1,length(t)) .* (1+freqrise*t);                              % compute frequency change
i = a * sin(2*pi*f.*repmat(t,size(f,1),1)).*exp(-d*t);                     % compute impulse response of bubbles

% move the bubbles in a random position.
% 1 - Define a threshold of amplitude of the sound.
% 2 - Find number of samples in which there is sound
% 3 - Use Circshift to shift the bubble to a random position on the buffer
ampthr = 0.001;                                                            % lowest amplitude for the bubble
samplesbubble = sum(abs(i) > ampthr,2);                                    % number of samples in where the bubble sounds

i_shifted = zeros(size(i));                                                % initialize matrix for time-shifted bubbles
for ndx = 1:size(i,1)                                                      % for each bubble,
    shift_samples = round(rand(1)*(length(i)-samplesbubble(ndx)));         % compute samples to shift
    i_shifted(ndx,:) = circshift(i(ndx,:),[0 shift_samples]);              % perform shift
end

% run this to see how the shift looks like:
% ndddx = 2;
% plot(i_shifted(ndddx,:));hold on
% plot(i(ndddx,:));
% plot([samplesbubble(ndddx),samplesbubble(ndddx)],[-1,1])
% legend('shifted','original','bubble length ')

% all bubbles!
outputsound = sum(i_shifted,1);                                            % put all bubbles in a vector
soundsc(outputsound,fs)                                                    % play it!

%% Simulate the sound of a sink leaking

