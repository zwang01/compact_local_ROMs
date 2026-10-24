function fig_pc_main()
% FIG_PC_MAIN  History figures of the paper from the runs of
% run_pc_baseline_chain (results_dir()/pc_*.mat), written to results_dir()/figures:
%   <stem>_Jhist.pdf       objective histories: preconditioned full order,
%                          local ROM, compact local ROM (indistinguishable)
%   <stem>_eta.pdf         calibrated indicator of the local ROM, first kmax
%                          iterations: accepted (blue), fallbacks (red), tau_eta
%   <stem>_eta_compact.pdf the same for the compact local ROM
% stems: body_cannula, body_nse40, body_nse200 (kmax = 1000); body3d_Jhist40 (histories only).
set(0, 'DefaultAxesFontSize', 12); LEG = 13;
cases = {'cannula', 'body_cannula', 'Stokes cannula', '$J(\phi_k)$'; 're40', 'body_nse40', '$Re = 40$', 'drag $\mathcal{D}(\phi_k)$'; ...
         're200', 'body_nse200', '$Re = 200$', 'drag $\mathcal{D}(\phi_k)$'; 'body3d', 'body3d', 'Pironneau''s body', 'drag $\mathcal{D}(\phi_k)$'};
for i = 1:size(cases, 1)
    c = cases{i,1}; f = results_dir(['pc_' c '_fom.mat']); if ~isfile(f), continue, end
    F = load(f); F = F.fom; R = load(results_dir(['pc_' c '_rom.mat'])); R = R.rom; C = load(results_dir(['pc_' c '_compact.mat'])); C = C.cpt;
    stem = cases{i,2}; ttl = cases{i,3}; ylab = cases{i,4};
    fh = figure('Position', [100 100 520 340], 'Visible', 'off');
    if strcmp(c, 'body3d'), pl = @plot; else, pl = @semilogy; end
    pl(0:numel(F.J)-1, F.J, 'k-', 'LineWidth', 1.0); hold on
    pl(0:numel(R.J)-1, R.J, 'r--', 'LineWidth', 1.6);
    pl(0:numel(C.J)-1, C.J, 'b:', 'LineWidth', 1.6);
    Jall = [F.J(:); R.J(:); C.J(:)]; ylim([0.95*min(Jall), 1.05*max(Jall)]);
    grid on; xlabel('iteration $k$', 'Interpreter', 'latex'); ylabel(ylab, 'Interpreter', 'latex');
    legend({'preconditioned full order', 'local ROM', 'compact local ROM'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', LEG);
    if ~strcmp(c, 'body3d'), title(ttl, 'Interpreter', 'latex'); end
    if strcmp(c, 'body3d'), fn = 'body3d_Jhist40.pdf'; else, fn = [stem '_Jhist.pdf']; end
    exportgraphics(fh, results_dir('figures', fn), 'ContentType', 'vector'); close(fh);
    if strcmp(c, 'body3d'), continue, end
    eta_panel(R, ttl, results_dir('figures', [stem '_eta.pdf']), LEG, 1000, 'local ROM');
    eta_panel(C, ttl, results_dir('figures', [stem '_eta_compact.pdf']), LEG, 1000, 'compact local ROM');
end
fprintf('figures written to %s\n', results_dir('figures'));
end

function eta_panel(R, ttl, file, LEG, kmax, which)
eta = R.eta; md = R.mode; nR = numel(eta); kk = (0:nR-1)'; kmax = min(kmax, nR-1); kp = kk <= kmax;
acc = kp & md == 0 & isfinite(eta); rej = kp & md == 1 & isfinite(eta);
fh = figure('Position', [100 100 520 340], 'Visible', 'off');
semilogy(kk(acc), eta(acc), 'b.', 'MarkerSize', 5); hold on
semilogy(kk(rej), eta(rej), 'rv', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
semilogy(kk(kp), R.tol(kp), 'k--', 'LineWidth', 1.4);
grid on; xlim([0 kmax]); xlabel('iteration $k$', 'Interpreter', 'latex'); ylabel('$\tilde\eta_k$', 'Interpreter', 'latex');
legend({'ROM accepted', 'fallback', 'tolerance $\tau_\eta$'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', LEG);
title([ttl ', ' which], 'Interpreter', 'latex');
exportgraphics(fh, file, 'ContentType', 'vector'); close(fh);
end
