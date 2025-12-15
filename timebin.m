function [tb, Ys] = timebin(t, Y, bin, sigma)

    tb = (t(1):bin:t(end))';
    Yz = interp1(t, Y, tb, 'previous', 'extrap');

    if nargin < 4 || sigma <= 0
        Ys = Yz;
        return
    end

    sigma_bin = sigma / bin;
    ksz = ceil(6*sigma_bin);
    if mod(ksz,2)==0, ksz=ksz+1; end

    x = (-floor(ksz/2):floor(ksz/2))';
    g = exp(-x.^2/(2*sigma_bin^2));
    g = g / sum(g);

    % --- smoothing with replicate padding ---
    k = floor(numel(g)/2);
    Ypad = [repmat(Yz(1,:),k,1); Yz; repmat(Yz(end,:),k,1)];
    
    Ys = zeros(size(Yz));
    for n = 1:size(Yz,2)
        Ys(:,n) = conv(Ypad(:,n), g, 'valid');
    end

end
