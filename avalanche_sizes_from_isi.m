function sizes = avalanche_sizes_from_isi(isi, thr)
% AVALANCHE_SIZES_FROM_ISI
% 
% isi : 1D array of inter-spike intervals
% thr : threshold (same units as isi)
%
% sizes : vector of avalanche sizes (number of spikes per avalanche)

    isi = isi(:);              % ensure column
    assert(all(isi >= 0), 'ISI must be nonnegative');

    % boundaries: ISI > thr ends an avalanche
    boundaries = [true; isi > thr];   % true = start of new avalanche

    % count spikes per avalanche
    idx = cumsum(boundaries);         % avalanche index
    sizes = accumarray(idx, 1);       % count spikes

end