% FIG_SEC52  Isolated-body geometry panels from run_sec52_body2d.
% Also writes older two-method history panels; run fig_pc_main afterward
% to replace them with the published preconditioned-FOM/local/compact panels.
%   Figures 5, 6: bodyNSE_init_bigdom.pdf (initial phi), dragconv_Re<Re>_stream_init.pdf,
%                 sub_Re<Re>_stream_final.pdf (streamlines at the initial disk and at
%                 the final design of the full-order run, recomputed at full order)
%   Figure 7    : body_nse<Re>_Jhist.pdf, body_nse<Re>_eta.pdf
% Also prints the drag reduction, the aspect ratio of the body (bounding box of
% the cells with phi < 1/2) and the reattachment point of the wake on the axis.
%   matlab -batch "setup_paths; fig_sec52"
set(0, 'DefaultAxesFontSize', 12);
if ~exist('Res', 'var'), Res = [40 200]; end
for Re = Res
    F = load(results_dir(sprintf('sec52_Re%d_fom.mat', Re))); F = F.res;
    R = load(results_dir(sprintf('sec52_Re%d_rom.mat', Re))); R = R.res;
    table_rows(sprintf('Fresh-LU comparison: Re = %d', Re), F, R);
    fprintf('  full order: drag %.4f -> %.6f (%.1f%% reduction) in %d iterations\n', ...
        F.J(1), F.J(end), 100*(1 - F.J(end)/F.J(1)), numel(F.J));
    hist_panels(R, F, sprintf('$Re = %d$', Re), 'drag $\mathcal{D}(\phi_k)$', sprintf('body_nse%d', Re));

    par = F.par; xc = F.xc; yc = F.yc; h = par.h; [Xc, ~] = ndgrid(xc, yc);
    if Re == Res(1)                               % initial phase field (the same for both Re)
        fh = figure('Position', [100 100 640 200], 'Visible', 'off');
        imagesc(xc, yc, F.phi0'); axis xy equal tight; colorbar; hold on
        contour(xc, yc, F.phi0', [0.5 0.5], 'r', 'LineWidth', 1.2);
        title('Initial $\phi$', 'Interpreter', 'latex');
        exportgraphics(fh, results_dir('figures', 'bodyNSE_init_bigdom.pdf'), 'ContentType', 'vector'); close(fh);
    end
    for q = 0:1
        if q == 0, phi = F.phi0; tag = 'initial'; fn = sprintf('dragconv_Re%d_stream_init.pdf', Re);
                   ttl = sprintf('Streamlines at the initial design, $Re=%d$', Re);
        else,      phi = F.phi;  tag = 'final';   fn = sprintf('sub_Re%d_stream_final.pdf', Re);
                   ttl = sprintf('Streamlines at the final design, $Re=%d$ (%d iterations)', Re, numel(F.J));
        end
        [D, o] = solve_at(par, phi, Re);
        [psi, xn, yn] = mac_streamfunction(par, o.u, o.v);
        stream_panel(results_dir('figures', fn), ttl, xc, yc, phi, psi, xn, yn);
        % wake: last zero crossing of the axial velocity on the centreline behind the body
        uc = 0.5*(o.u(1:end-1,:) + o.u(2:end,:));
        sol = phi < 0.5; [ii, jj] = find(sol); xback = xc(max(ii));
        [~, jm] = min(abs(yc)); ucl = 0.5*(uc(:,jm) + uc(:,max(jm-1,1)));
        neg = find(xc(:) > xback & ucl < 0);
        if isempty(neg), wake = 'no recirculation'; else, wake = sprintf('recirculation reattaching at x = %+.2f', xc(max(neg)) + h/2); end
        fprintf('  %-7s design: drag %.6f, aspect ratio %.2f (length %.3f, height %.3f), %s\n', ...
            tag, D, (max(ii)-min(ii)+1)/(max(jj)-min(jj)+1), (max(ii)-min(ii)+1)*h, (max(jj)-min(jj)+1)*h, wake);
    end
end

function [D, o] = solve_at(par, phi, Re)
% full-order state at the design phi (continuation in nu for Re > 100)
S = brinkman_base(par); Uw = [];
if Re > 100
    pc = par;
    for nuc = par.nu*[8 4 2]
        pc.nu = nuc; [~, ~, oc] = brinkman_nse_solve_J(pc, brinkman_base(pc), phi, Uw); Uw = oc.U;
    end
end
[D, ~, o] = brinkman_nse_solve_J(par, S, phi, Uw);
end

function stream_panel(fname, ttl, xc, yc, phi, psi, xn, yn)
% equally spaced contours of the stream function, the solid filled in grey
fh = figure('Position', [100 100 640 200], 'Visible', 'off');
contour(xn, yn, psi', linspace(min(psi(:)), max(psi(:)), 31), 'k', 'LineWidth', 0.6); hold on
Cc = contourc(xc, yc, phi', [0.5 0.5]); q = 1;
while q < size(Cc, 2)
    nseg = Cc(2,q); seg = Cc(:, q+1:q+nseg);
    fill(seg(1,:), seg(2,:), [0.80 0.80 0.84], 'EdgeColor', 'none');
    q = q + nseg + 1;
end
contour(xc, yc, phi', [0.5 0.5], 'r', 'LineWidth', 1.3);
axis xy equal; xlim([xn(1) xn(end)]); ylim([yn(1) yn(end)]); set(gca, 'Layer', 'top'); box on
cb = colorbar; cb.Visible = 'off';
title(ttl, 'Interpreter', 'latex');
exportgraphics(fh, fname, 'ContentType', 'vector'); close(fh);
end
