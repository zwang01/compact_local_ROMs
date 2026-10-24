% FIG1_SVD  Figure 1 of the paper from the results of run_fig1_svd: normalized
% singular values of the global snapshot set (solid black) and of three
% windows of w_max = 40 consecutive snapshots (dashed), for the cannula and the
% isolated body at Re = 40 and 200. Writes results/figures/svd_flow_<case>.pdf
% and prints the number of modes needed for 1e-3, 1e-6 and 1e-10 of sigma_1.
%   matlab -batch "setup_paths; fig1_svd"
set(0, 'DefaultAxesFontSize', 12);
cases = {'cannula', 'nse40', 'nse200'};
titles = {'Stokes cannula', 'Isolated body, $Re=40$', 'Isolated body, $Re=200$'};
cols = {[0 0.45 0.74], [0.85 0.33 0.1], [0.47 0.67 0.19]};
need = @(s, tol) find(s < tol, 1, 'first');
for c = 1:numel(cases)
    f = results_dir(sprintf('fig1_svd_%s.mat', cases{c}));
    if ~isfile(f), fprintf('%-8s: no data (run run_fig1_svd)\n', cases{c}); continue, end
    D = load(f);
    % the three panels share one row of the page: small figures with large fonts
    fh = figure('Position', [100 100 320 285], 'Visible', 'off');
    semilogy(1:numel(D.sg), D.sg, 'k-', 'LineWidth', 2.0); hold on
    leg = {'global'};
    for q = 1:numel(D.sw)
        semilogy(1:numel(D.sw{q}), D.sw{q}, '--', 'Color', cols{q}, 'LineWidth', 1.8);
        leg{end+1} = sprintf('$k=%d$', D.kw(q)); %#ok<SAGROW>
    end
    yline(1e-10, 'k:', 'LineWidth', 1.2);
    grid on; xlim([1 min(120, numel(D.sg))]); ylim([1e-16 1]);
    set(gca, 'FontSize', 17);
    xlabel('index $j$', 'Interpreter', 'latex', 'FontSize', 19);
    title(titles{c}, 'Interpreter', 'latex', 'FontSize', 14);
    legend(leg, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 15);
    % only the leftmost panel carries the y-label; the others take its plot box
    if c == 1
        ylabel('$\sigma_j/\sigma_1$', 'Interpreter', 'latex', 'FontSize', 19);
        drawnow; pos = get(gca, 'Position');
    elseif exist('pos', 'var')
        set(gca, 'Position', pos);
    end
    exportgraphics(fh, results_dir('figures', sprintf('svd_flow_%s.pdf', cases{c})), 'ContentType', 'vector');
    close(fh);
    fprintf('%-8s global %3s/%3s/%3s   windows', cases{c}, ...
        num2str(need(D.sg,1e-3)), num2str(need(D.sg,1e-6)), num2str(need(D.sg,1e-10)));
    for q = 1:numel(D.sw)
        fprintf('  k=%4d: %s/%s/%s', D.kw(q), num2str(need(D.sw{q},1e-3)), ...
            num2str(need(D.sw{q},1e-6)), num2str(need(D.sw{q},1e-10)));
    end
    fprintf('\n');
end
fprintf('(modes to reach 1e-3 / 1e-6 / 1e-10 of the leading singular value)\n');
