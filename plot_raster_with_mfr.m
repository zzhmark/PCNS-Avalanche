function f = plot_raster_with_mfr(Tb, nspk, mfr, ttot)
% PLOT_RASTER_WITH_MFR
%
% Tb   : time vector (same length as rows of nspk), in ms
% nspk : (T × N) spike matrix (nonzero = spike)
% mfr  : mean firing rate over time (same length as Tb)
% ttot : total time for x-axis limit (ms)

    f = figure;

    % ---------- raster ----------
    yyaxis left
    hold on

    [t_ind, n_ind] = find(nspk);
    plot(Tb(t_ind), n_ind, '.', 'Color', [0.5 0.5 0.5]);

    xlim([0 ttot])
    xlabel('Time (ms)')

    % ---------- mean firing rate ----------
    yyaxis right
    plot(Tb, mfr, 'b-', 'LineWidth', 2);
    ylabel('Mean firing rate (Hz)')
    ylim([0 100])

    % ---------- axis styling ----------
    ax = gca;
    ax.Box = 'off';
    ax.XAxisLocation = 'bottom';

    % keep axes visible but consistent
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

    axis square
    hold off
end
