function [tb, Yb] = timebin_sum(t, Y, bin)

    t = t(:);
    assert(size(Y,1) == numel(t), 'Y must be (T,N)');

    % ---- define bins ----
    edges = t(1):bin:t(end)+bin;
    tb = edges(1:end-1)' + bin/2;
    B  = numel(tb);
    N  = size(Y,2);

    % ---- assign samples to bins ----
    idx = discretize(t, edges);

    % ---- bin-wise sum ----
    Yb = zeros(B, N);
    for n = 1:N
        Yb(:,n) = accumarray(idx, Y(:,n), [B 1], @sum, 0);
    end
end
