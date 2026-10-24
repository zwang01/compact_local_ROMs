% RUN_APPC_VALIDATION  Checks of Algorithm 2 in Appendix C of the paper, all on
% the isolated-body benchmark at Re = 40 (Section 5.2 setting) with a matched
% budget of 300 design iterations (k = 0..300) and the stopping rule off.
% Results in results/appC_<label>.mat; tables_appC prints Tables C.2-C.4.
%
%   part 'refine'   (Table C.2) h = 1/40, 1/80, 1/160 with eps = 2h and
%                   alpha_max = 2.5e4*nu*(80h)^{-1}: Algorithm 2; the full-order
%                   optimization with a fresh factorization at every solve
%                   (Algorithm 1 as in Section 5); and the full-order
%                   optimization through the retained factors and GMRES, the
%                   full-order solver of Algorithm 2         (about 70 min)
%   part 'ablate'   (Table C.3) h = 1/80, one component changed at a time:
%                   Euclidean gate at 3e-3, verification by the backstop only,
%                   fresh factorization at every fallback, window 2*w_max,
%                   cumulative history (the reference row is refine_rom_N80)
%                                                                  (about 6 min)
%   part 'validate' (Table C.4) every 5th accepted reduced iteration recomputed
%                   at full order                                  (about 3 min)
%
%   matlab -batch "setup_paths; run_appC_validation"
%   matlab -batch "setup_paths; parts = {'ablate'}; run_appC_validation"
%
% These runs use the current defaults (bound-constrained volume projection,
% QR least squares). Algorithm 2 and the retained-factor baseline reuse the
% last Newton factorization for the adjoint (adj_reuse = true); the
% fresh-factorization baseline factors the converged Jacobian (adj_reuse =
% false), as Algorithm 1 does in Section 5.
if ~exist('parts', 'var'), parts = {'refine', 'ablate', 'validate'}; end
base = struct('Re', 40, 'max_itr', 300, 'stop_tol', 0, 'adj_reuse', true);
runs = {};                                   % {label, cfg}
if any(strcmp(parts, 'refine'))
    for N = [40 80 160]
        mesh = struct('N', N, 'eps_fac', 2, 'amax_fac', N/80);
        runs(end+1, :) = {sprintf('refine_rom_N%d', N),      merge_struct(base, mesh)}; %#ok<SAGROW>
        runs(end+1, :) = {sprintf('refine_fom_N%d', N),      merge_struct(base, mesh, struct('tol_eta', 0, 'fallback', 'refactor', 'adj_reuse', false))}; %#ok<SAGROW>
        runs(end+1, :) = {sprintf('refine_fomkry_N%d', N),   merge_struct(base, mesh, struct('tol_eta', 0, 'fallback', 'precond'))}; %#ok<SAGROW>
    end
end
if any(strcmp(parts, 'ablate'))
    runs(end+1, :) = {'ablate_eucl',     merge_struct(base, struct('gate', 'eucl', 'tol_eta', 3e-3, 'tol_eta_max', 9e-3, 'tol_eta_min', 3e-7))};
    runs(end+1, :) = {'ablate_backstop', merge_struct(base, struct('verify', 'backstop'))};
    runs(end+1, :) = {'ablate_refactor', merge_struct(base, struct('fallback', 'refactor'))};
    runs(end+1, :) = {'ablate_window80', merge_struct(base, struct('w_max', 80))};
    runs(end+1, :) = {'ablate_cumul',    merge_struct(base, struct('w_max', 1e6))};
end
if any(strcmp(parts, 'validate'))
    runs(end+1, :) = {'validate_es5',    merge_struct(base, struct('errstudy', true, 'es_every', 5))};
end
for i = 1:size(runs, 1)
    fprintf('\n[appC] ===== %d/%d: %s =====\n', i, size(runs, 1), runs{i,1});
    res = shapeopt_body2d(runs{i,2}); %#ok<NASGU>
    save(results_dir(['appC_' runs{i,1} '.mat']), 'res');
end
