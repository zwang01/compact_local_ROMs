function hist_panels(R, F, ttl, ylab, stem, kmax)
% HIST_PANELS  The two panels of Figures 3 and 7 of the paper for one case:
%   <stem>_Jhist.pdf : objective histories, Algorithm 2 (R) against Algorithm 1 (F)
%   <stem>_eta.pdf   : the screening indicator max(eta_s, eta_a) of Algorithm 2
%                      over the first kmax iterations: accepted reduced iterations
%                      (blue), iterations solved at full order (red), tolerance
%                      tau_eta (dashed)
% R, F: result structs of the shapeopt_* drivers. Files go to results/figures/.
if nargin < 6, kmax = 1000; end
nR = numel(R.J); nF = numel(F.J);
fh = figure('Position', [100 100 520 340], 'Visible', 'off');
semilogy(0:nR-1, R.J, 'r--', 'LineWidth', 1.6); hold on
semilogy(0:nF-1, F.J, 'k-', 'LineWidth', 1.0);
Jall = [R.J(:); F.J(:)]; ylim([0.95*min(Jall), 1.05*max(Jall)]);
grid on; xlabel('iteration $k$', 'Interpreter', 'latex'); ylabel(ylab, 'Interpreter', 'latex');
legend({'ROM', 'full order'}, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
title(ttl, 'Interpreter', 'latex');
exportgraphics(fh, results_dir('figures', [stem '_Jhist.pdf']), 'ContentType', 'vector'); close(fh);

eta = R.eta; md = R.mode; kk = (0:nR-1)'; kmax = min(kmax, nR-1); kp = kk <= kmax;
acc = kp & md == 0 & isfinite(eta); rej = kp & md == 1 & isfinite(eta);
fh = figure('Position', [100 100 520 340], 'Visible', 'off');
semilogy(kk(acc), eta(acc), 'b.', 'MarkerSize', 5); hold on
semilogy(kk(rej), eta(rej), 'rv', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
semilogy(kk(kp), R.tol(kp), 'k--', 'LineWidth', 1.4);
grid on; xlim([0 kmax]); xlabel('iteration $k$', 'Interpreter', 'latex'); ylabel('$\tilde\eta_k$', 'Interpreter', 'latex');
legend({'ROM accepted', 'fallback', 'tolerance $\tau_\eta$'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 13);
title(ttl, 'Interpreter', 'latex');
exportgraphics(fh, results_dir('figures', [stem '_eta.pdf']), 'ContentType', 'vector'); close(fh);
end
