function [MSE, estimatedperiod] = l5_combfilter(x,maxperiod)
N = length(x);
a = 0.9;                                                                   % periodicity coefficient
MSE = zeros(1,maxperiod);                                                  % initialize MSE vector
for period = 1:maxperiod                                                   % from 1 to the maximum period value
    idx1 = 1:N-period;                                                     % indexes for not delayed signal
    idx2 = idx1+period;                                                    % indexes for delayed signal
    e = power(x(idx1) - a * x(idx2),2);                                    % compute the error between the n sample
    MSE(period) = sum(e)/(N - period);                                     % compute mean squared error for this tau (period)
%     MSE(period) = sum(e)/(N - period);                                   % compute mean squared error for this tau (period)
end
[~,estimatedperiod] = min(MSE);                                            % find the minimum mean squared error
end