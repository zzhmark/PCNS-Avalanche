function isi = compute_isi(t, E)
    t = t(:);
    te = t(E ~= 0);   % event times
    isi = diff(te);
end
