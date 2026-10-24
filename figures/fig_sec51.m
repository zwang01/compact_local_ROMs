% FIG_SEC51  Cannula geometry and error-study panels from run_sec51_cannula.
% Also writes an older two-method history panel; run fig_pc_main afterward
% to replace it with the published preconditioned-FOM/local/compact panel.
%   Figure 2: cannula_{phi,speed,stream}_{init,mid}.pdf, cannulaL_{phi,speed,stream}_final.pdf
%             (initial duct and the two designs kept by the full-order run, k = 3000
%             and k = 8314; the flow is recomputed at full order at each design)
%   Figure 3: body_cannula_Jhist.pdf, body_cannula_eta.pdf
%   Figure 4: body_errlaw.pdf, body_effectivity.pdf (error study)
%   matlab -batch "setup_paths; fig_sec51"
set(0, 'DefaultAxesFontSize', 12);
F = load(results_dir('sec51_fom.mat')); F = F.res;
R = load(results_dir('sec51_rom.mat')); R = R.res;

% ---------------- Table 1 and Figure 3 ----------------
table_rows('Fresh-LU comparison: cannula', F, R);
fprintf('  full order: viscous dissipation %.2f%% and J %.2f%% below the initial duct after %d iterations\n', ...
    100*(1 - F.Jvisc(end)/F.Jvisc(1)), 100*(1 - F.J(end)/F.J(1)), numel(F.J));
hist_panels(R, F, 'Stokes cannula', '$J(\phi_k)$', 'body_cannula');

% ---------------- Figure 2: designs and flows of the full-order run ----------------
par = F.par; S = brinkman_base(par); xc = F.xc; yc = F.yc;
designs = {F.phi0, F.phi_saved(1).phi, F.phi_saved(end).phi};    % k = 0, 3000, 8314 in the paper
files  = {'cannula_phi_init', 'cannula_phi_mid', 'cannulaL_phi_final'; ...
          'cannula_speed_init', 'cannula_speed_mid', 'cannulaL_speed_final'; ...
          'cannula_stream_init', 'cannula_stream_mid', 'cannulaL_stream_final'};
titles = {'Initial $\phi$: L-shaped duct', sprintf('After %d iterations', F.phi_saved(1).k), 'Stationary design'};
for p = 1:3
    phi = designs{p};
    [~, ~, o] = brinkman_stokes_solve_J(par, S, phi);
    if p == 1, Jv0 = o.Jvisc; end
    spd = sqrt((0.5*(o.u(1:end-1,:) + o.u(2:end,:))).^2 + (0.5*(o.v(:,1:end-1) + o.v(:,2:end))).^2);
    [psi, xn, yn] = mac_streamfunction(par, o.u, o.v);
    % phase field with the phi = 1/2 contour
    fh = figure('Position', [100 100 330 300], 'Visible', 'off');
    imagesc(xc, yc, phi'); axis xy equal tight; colorbar; hold on
    contour(xc, yc, phi', [0.5 0.5], 'r', 'LineWidth', 1.3);
    set(gca, 'FontSize', 14); title(titles{p}, 'Interpreter', 'latex', 'FontSize', 15);
    exportgraphics(fh, results_dir('figures', [files{1,p} '.pdf']), 'ContentType', 'vector'); close(fh);
    % speed with the solid (phi < 1/2) masked, so the colour scale spans the fluid
    sm = spd; sm(phi < 0.5) = NaN;
    fh = figure('Position', [100 100 330 300], 'Visible', 'off');
    imagesc(xc, yc, sm', 'AlphaData', ~isnan(sm')); axis xy equal tight
    set(gca, 'Color', [0.80 0.80 0.84]); colormap(gca, turbo); colorbar; hold on
    contour(xc, yc, phi', [0.5 0.5], 'r', 'LineWidth', 1.3);
    set(gca, 'FontSize', 14); title('Speed $|u|$', 'Interpreter', 'latex', 'FontSize', 15);
    exportgraphics(fh, results_dir('figures', [files{2,p} '.pdf']), 'ContentType', 'vector'); close(fh);
    % streamlines (equally spaced contours of the stream function), solid masked
    % by an alpha overlay (the phi = 1/2 contour is open here: the solid touches the box)
    fh = figure('Position', [100 100 330 300], 'Visible', 'off');
    contour(xn, yn, psi', linspace(min(psi(:)), max(psi(:)), 11), 'k', 'LineWidth', 0.7); hold on
    Im = imagesc(xc, yc, zeros(par.Ny, par.Nx)); set(Im, 'AlphaData', double(phi' < 0.5));
    colormap(gca, [0.80 0.80 0.84]);
    contour(xc, yc, phi', [0.5 0.5], 'r', 'LineWidth', 1.3);
    axis xy equal; xlim([xn(1) xn(end)]); ylim([yn(1) yn(end)]); set(gca, 'Layer', 'top'); box on
    cb = colorbar; cb.Visible = 'off';           % same width as the other panels
    set(gca, 'FontSize', 14); title('Streamlines', 'Interpreter', 'latex', 'FontSize', 15);
    exportgraphics(fh, results_dir('figures', [files{3,p} '.pdf']), 'ContentType', 'vector'); close(fh);
    fprintf('Figure 2, %-28s viscous dissipation %.6e (%.1f%% below the initial duct)\n', titles{p}, o.Jvisc, 100*(1 - o.Jvisc/Jv0));
end
fprintf('Figure 2, right column: the design at k = %d of %d iterations\n', F.phi_saved(end).k, numel(F.J));

% ---------------- Figure 4: error study ----------------
E = load(results_dir('sec51_errstudy.mat')); E = E.res; es = E.es;
ok = isfinite(es.err) & es.err > 0 & es.dJ > 0 & es.etaP > 0;
etaP = es.etaP(ok); dJ = es.dJ(ok); err = es.err(ok); erra = es.erra(ok);
etaJ = es.etaJ(ok); dJs = es.dJsigned(ok); JF = abs(dJs)./dJ;
corr = abs(dJs - etaJ)./JF;                       % objective error after the DWR correction
effJ = etaJ./dJs;                                 % effectivity of the dual-weighted estimate
P = polyfit(log(etaP), log(dJ), 1); C = exp(P(2)); p = P(1);
fh = figure('Position', [100 100 560 450], 'Visible', 'off');
loglog(etaP, dJ, '.', 'Color', [0 0.45 0.74], 'MarkerSize', 11); hold on
loglog(etaP, corr, '.', 'Color', [0.85 0.33 0.1], 'MarkerSize', 11);
ev = logspace(log10(min(etaP)), log10(max(etaP)), 50);
sel = etaP >= 1e-6 & corr > 0; C2 = median(corr(sel)./etaP(sel).^2);   % second-order reference above round-off
loglog(ev, C*ev.^p, 'k-', 'LineWidth', 1.8); loglog(ev, C2*ev.^2, 'k-.', 'LineWidth', 1.3);
xline(E.cfg.tol_eta, '-.', 'Color', [0.30 0.10 0.55], 'LineWidth', 2.0);
xlim([0.85*min(etaP), 1.4*max(max(etaP), E.cfg.tol_eta)]);
grid on; xlabel('reference-scaled indicator $\eta_P$', 'Interpreter', 'latex'); ylabel('relative objective error', 'Interpreter', 'latex');
legend({'$|J_h(U)-J_h(U_r)|/|J_h(U)|$', 'after DWR correction', sprintf('fit $%.2f\\,\\eta_P^{%.2f}$', C, p), ...
        'second-order reference', '$\tau_\eta$'}, 'Interpreter', 'latex', 'Location', 'west', 'FontSize', 12);
set(gca, 'FontSize', 15);
exportgraphics(fh, results_dir('figures', 'body_errlaw.pdf'), 'ContentType', 'vector'); close(fh);

kk = es.k(ok); rE = es.etaE(ok)./err; rC = es.etaC(ok)./err; rP = etaP./err;
fh = figure('Position', [100 100 560 450], 'Visible', 'off');
semilogy(kk, rE, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7); hold on
semilogy(kk, rC, '.', 'Color', [0 0.45 0.74], 'MarkerSize', 7);
semilogy(kk, effJ, '.', 'Color', [0 0 0], 'MarkerSize', 5);
semilogy(kk, rP, '.', 'Color', [0.85 0.33 0.1], 'MarkerSize', 6);
yline(1, 'k--', 'LineWidth', 1.0);
grid on; xlim([0 max(kk)]); ylim([0.2 3e6]); xlabel('iteration $k$', 'Interpreter', 'latex'); ylabel('effectivity', 'Interpreter', 'latex');
cols = {[0.6 0.6 0.6], [0 0.45 0.74], [0 0 0], [0.85 0.33 0.1]}; hL = gobjects(4,1);
for q = 1:4, hL(q) = semilogy(NaN, NaN, '.', 'Color', cols{q}, 'MarkerSize', 20); end   % legend-only markers
lg = legend(hL, {'Euclidean $\eta^{E}$', 'calibrated $\tilde\eta$', 'dual-weighted $\eta_J$', 'reference-scaled $\eta_P$'}, ...
            'Interpreter', 'latex', 'Location', 'north', 'FontSize', 13);
lg.NumColumns = 2; set(gca, 'FontSize', 15);
exportgraphics(fh, results_dir('figures', 'body_effectivity.pdf'), 'ContentType', 'vector'); close(fh);

fprintf(['\nError study (Figure 4): %d of %d iterations accepted\n' ...
         '  true state error median %.1e max %.1e; adjoint error median %.1e max %.1e\n' ...
         '  objective error = %.2f eta_P^%.2f; median %.1e, after the DWR correction %.1e\n' ...
         '  effectivity medians: Euclidean %.3g, calibrated %.3g (fraction below 1: %.3f), exact %.3f, DWR %.4f [%.3f, %.3f]\n'], ...
    numel(etaP), numel(E.J), median(err), max(err), median(erra), max(erra), C, p, median(dJ), median(corr), ...
    median(rE), median(rC), mean(rC < 1), median(rP), median(effJ), min(effJ), max(effJ));
