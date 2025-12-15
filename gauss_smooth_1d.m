function Y = gauss_smooth_1d(Y, sigma)
% GAUSS_SMOOTH_1D  1D Gaussian smoothing with replicate boundary
%
% Y     : (T,N) or (T,1) time series
% sigma : Gaussian sigma (in samples)
%
% Y     : smoothed output (same size as input)

    if sigma <= 0
        return
    end

    % --- build Gaussian kernel ---
    ksz = ceil(6*sigma);
    if mod(ksz,2)==0, ksz = ksz+1; end
    x = (-floor(ksz/2):floor(ksz/2))';
    g = exp(-x.^2/(2*sigma^2));
    g = g / sum(g);

    % --- replicate padding ---
    k = floor(numel(g)/2);
    Ypad = [ ...
        repmat(Y(1,:),  k, 1); ...
        Y; ...
        repmat(Y(end,:),k, 1) ...
    ];

    % --- 1D convolution along time ---
    Y = zeros(size(Y));
    for n = 1:size(Y,2)
        Y(:,n) = conv(Ypad(:,n), g, 'valid');
    end
end
