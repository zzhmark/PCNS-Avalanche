function f = plot_avalanche_isi_fit(avalanche_sizes, avalanche_sizes2, isi_sizes, p_hat, beta, lambda_hat, dt)

% PLOT_AVALANCHE_WITH_INSET
%
% avalanche_sizes : vector of avalanche sizes
% isi_sizes       : vector of ISIs
% p_hat           : geometric parameter
% beta            : power-law exponent
% x_min           : minimum size for power-law
% lambda_hat      : exponential rate for ISI
% dt              : time bin (seconds)
% outfile         : output filename (e.g. 'fig3/F.png')

    f = figure;
    hold on;

    % --- avalanche size distribution (scatter) ---
    [counts, edges] = histcounts(avalanche_sizes, 1:100, ...
        'Normalization','probability');
    centers = sqrt(edges(1:end-1).*edges(2:end));

    scatter(centers, counts, 10, 'k', 'filled');
    set(gca,'XScale','log','YScale','log')
    axis square

    ylabel('Probability');
    xlabel('Avalanche Size');
    xlim([1 100])
    ylim([1e-6 1]);

    % --- geometric curve ---
    s_fit = 1:max(avalanche_sizes);
    plot(s_fit, p_hat .* (1 - p_hat).^s_fit, ...
        'r', 'LineWidth', 2);

    % --- power-law curve ---
    s_fit = 1:max(avalanche_sizes);
    C = 1 / sum((avalanche_sizes2).^(-beta));
    plot(s_fit, s_fit.^(-beta) * C, 'b', 'LineWidth', 2);

    % --- annotation (anchored to figure) ---
    annotation('textbox', [0.7 0.8 0.2 0.1], ...
        'String', sprintf('$\\Delta t=%.3f$\n$\\beta=%.2f$', dt, beta), ...
        'Interpreter','latex', ...
        'FontSize',14, ...
        'FitBoxToText','on', ...
        'EdgeColor','none');

    hold off;

    % ---------- inset ----------
    ax1 = gca;
    p = ax1.Position;

    ax2 = axes('Position', [ ...
        p(1) + 0.18*p(3), ...
        p(2) + 0.08*p(4), ...
        0.35*p(3), ...
        0.35*p(4)]);

    box on
    hold(ax2, 'on')

    % --- ISI distribution ---
    [counts, edges] = histcounts(isi_sizes, 100, 'Normalization','count');
    centers = sqrt(edges(1:end-1).*edges(2:end));

    scatter(centers, counts, 5, 'b', 'filled');
    set(gca, 'YScale','log')
    axis(ax2, 'square')
    ylim([0.1 1e4]);
    x_fit = 0:max(isi_sizes);
    pdf_fit = lambda_hat * exp(-lambda_hat * x_fit);
    plot(ax2, x_fit, pdf_fit * numel(isi_sizes) * mean(diff(edges)), 'g', 'LineWidth', 2);

    xlabel(ax2, 'ISI');
    set(ax2, 'YScale','log', 'FontSize',8);

    hold(ax2, 'off')
end
