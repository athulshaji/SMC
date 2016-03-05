function [Z] = l5_createVandermondeMat(w,L,M)
% input:
%   w: fundamental frequency
%   L:number of harmonics
%   M: number of samples
% output:
%   Z: MxL vandermonde matrix

    z = exp(1i*w * (1:L));          % compute harmonics
    Z = zeros(M-1,L);               % initialize Vandermonde matrix
    for n = 1:M                     % for each sample
        Z(n,:) = z.^(n-1);          % fill Vandermonde matrix
    end
end

