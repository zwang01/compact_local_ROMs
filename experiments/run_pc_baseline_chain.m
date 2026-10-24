function run_pc_baseline_chain(phases, smoke)
% RUN_PC_BASELINE_CHAIN  Every timed experiment of the paper with the
% preconditioned full-order optimization as the baseline: for each case the
% full-order run with the retained-factor GMRES of the fallbacks (warm-started
% state and adjoint, refactored past m_max), to the stopping rule; then the
% original local ROM (Algorithm 2) and the compact local ROM
% (experiments/rom_efficiency) for the same number of iterations. Strictly
% sequential, `cool` seconds idle before every timed leg; run on an otherwise
% idle machine with BRK_RESULTS set. Results: results_dir()/pc_<case>_{fom,rom,compact}.mat.
%
%   phases: subset of {'cannula','re40','re200','body3d','appC','ownrule','audit'}
%   smoke : true -> tiny grids and iteration counts (checks the chain only)
%
%   set BRK_RESULTS=...; matlab -batch "setup_paths; addpath('experiments/rom_efficiency'); run_pc_baseline_chain"
if nargin < 1 || isempty(phases), phases = {'cannula', 're40', 're200', 'body3d', 'appC', 'ownrule', 'audit'}; end
if nargin < 2, smoke = false; end
if ischar(phases), phases = {phases}; end
cool = 300; if smoke, cool = 0; end
pub = struct('vol_project', false, 'lsq', 'normal', 'verbose', false);      % settings of the published runs
cmp = struct('w_max', 8, 'r_max', 8, 'pod_tol', 1e-7, 'eff_rescue', 2, 'eff_goal_rescue', 2);   % compact variant
nse_cmp = struct('fallback', 'precond', 'adj_reuse', true, 'pc_tol', 1e-10, 'pc_maxit', 30, 'pc_refresh_its', 20, ...
                 'eff_variant', 'mapped', 'eff_max_newton', 2, 'eff_newton_tol', 1e-8, 'eff_rescue_mode', 'direct');
nse_pc = struct('tol_eta', 0, 'fallback', 'precond', 'adj_reuse', true, 'eff_warm_adjoint', true);
t_chain = tic; L = @(varargin) fprintf('[chain %s %.0f min] %s\n', datestr(now, 'HH:MM'), toc(t_chain)/60, sprintf(varargin{:}));
L('results in %s (smoke = %d)', results_dir(), smoke);

cases = struct();
cases.cannula = struct('mk_fom', @(c) shapeopt_cannula(m(pub, struct('tol_eta', 0, 'fallback', 'precond'), c)), ...
                       'mk_rom', @(c) shapeopt_cannula(m(pub, c)), ...
                       'mk_cmp', @(c) eff_compact_cannula(m(pub, cmp, struct('pod_tol', 1e-10), c)), 'sz', struct('Nx', 94));
cases.re40    = struct('mk_fom', @(c) eff_shapeopt_goal2d(m(pub, nse_pc, struct('Re', 40), c)), ...
                       'mk_rom', @(c) shapeopt_body2d(m(pub, struct('Re', 40), c)), ...
                       'mk_cmp', @(c) eff_shapeopt_goal2d(m(pub, cmp, nse_cmp, struct('Re', 40), c)), 'sz', struct('N', 80));
cases.re200   = struct('mk_fom', @(c) eff_shapeopt_goal2d(m(pub, nse_pc, struct('Re', 200), c)), ...
                       'mk_rom', @(c) shapeopt_body2d(m(pub, struct('Re', 200), c)), ...
                       'mk_cmp', @(c) eff_shapeopt_goal2d(m(pub, cmp, nse_cmp, struct('Re', 200), c)), 'sz', struct('N', 80));
cases.body3d  = struct('mk_fom', @(c) shapeopt_body3d(m(pub, struct('tol_eta', 0, 'fallback', 'precond'), c)), ...
                       'mk_rom', @(c) shapeopt_body3d(m(pub, c)), ...
                       'mk_cmp', @(c) eff_compact_body3d(m(pub, cmp, struct('pod_tol', 1e-10), c)), 'sz', struct('N', 40));
if smoke
    cases.cannula.sz = struct('Nx', 47, 'max_itr', 60); cases.re40.sz = struct('N', 20, 'max_itr', 40);
    cases.re200.sz = struct('N', 20, 'max_itr', 40);    cases.body3d.sz = struct('N', 16, 'max_itr', 12);
end

% ---------------- main cases: PC FOM to the rule, both ROMs matched ----------------
for name = {'cannula', 're40', 're200', 'body3d'}
    if ~any(strcmp(phases, name{1})), continue, end
    cs = cases.(name{1}); sz = cs.sz;
    idle(cool, L, [name{1} ' PC FOM']);
    fom = cs.mk_fom(sz); save(results_dir(['pc_' name{1} '_fom.mat']), 'fom');
    L('>>> %s PC FOM: %d its, %.1f s, %d fact, %d GMRES solves (%d failed), J %.7e', name{1}, fom.iters, fom.time, fom.nfact, fom.n_pc, fom.n_pc_fail, fom.J(end));
    len = struct('max_itr', fom.iters - 1, 'stop_tol', 0, 'stop_dJ_tol', 0);
    idle(cool, L, [name{1} ' original ROM']);
    rom = cs.mk_rom(m(sz, len)); save(results_dir(['pc_' name{1} '_rom.mat']), 'rom');
    report(L, [name{1} ' original ROM'], rom, fom);
    idle(cool, L, [name{1} ' compact ROM']);
    cpt = cs.mk_cmp(m(sz, len)); save(results_dir(['pc_' name{1} '_compact.mat']), 'cpt');
    report(L, [name{1} ' compact ROM'], cpt, fom);
end

% ---------------- Appendix C: coupled refinement, 301 iterations ----------------
if any(strcmp(phases, 'appC'))
    grids = [40 80 160]; if smoke, grids = 20; end
    appc = struct('Re', 40, 'eps_fac', 2, 'max_itr', 300, 'stop_tol', 0, 'stop_dJ_tol', 0, 'vol_project', true, 'lsq', 'qr', 'verbose', false);
    if smoke, appc.max_itr = 30; end
    for N = grids
        g = m(appc, struct('N', N, 'amax_fac', N/80));
        idle(min(cool, 60), L, sprintf('appC N=%d PC FOM', N));
        fom = eff_shapeopt_goal2d(m(g, nse_pc)); save(results_dir(sprintf('pc_appC_N%d_fom.mat', N)), 'fom');
        L('>>> appC N=%d PC FOM: %.1f s, %d fact', N, fom.time, fom.nfact);
        idle(min(cool, 60), L, sprintf('appC N=%d original ROM', N));
        rom = shapeopt_body2d(g); save(results_dir(sprintf('pc_appC_N%d_rom.mat', N)), 'rom'); report(L, sprintf('appC N=%d original ROM', N), rom, fom);
        idle(min(cool, 60), L, sprintf('appC N=%d compact ROM', N));
        cpt = eff_shapeopt_goal2d(m(g, cmp, nse_cmp, struct('pc_maxit', 20))); save(results_dir(sprintf('pc_appC_N%d_compact.mat', N)), 'cpt');
        report(L, sprintf('appC N=%d compact ROM', N), cpt, fom);
    end
end

% ---------------- compact ROM under its own stopping rule (untimed) ----------------
if any(strcmp(phases, 'ownrule'))
    for name = {'cannula', 're40', 're200', 'body3d'}
        cs = cases.(name{1}); f = results_dir(['pc_' name{1} '_fom.mat']); if ~isfile(f), continue, end
        F = load(f); fom = F.fom;
        own = cs.mk_cmp(cs.sz); save(results_dir(['pc_' name{1} '_compact_ownrule.mat']), 'own');
        n = min(numel(own.J), numel(fom.J));
        L('>>> %s compact ROM under its own rule: %d its (PC FOM %d), gap at that iteration %.2e', name{1}, own.iters, fom.iters, abs(own.J(n) - fom.J(n))/abs(fom.J(n)));
    end
end

% ---------------- same-design audits of the compact runs (measurement only, untimed) ----------------
if any(strcmp(phases, 'audit'))
    aud = {'cannula', struct('max_itr', 2999, 'errstudy', true, 'es_every', 1); 're40', struct('errstudy', true, 'es_every', 20); ...
           're200', struct('errstudy', true, 'es_every', 20); 'body3d', struct('errstudy', true, 'es_every', 50)};
    for i = 1:size(aud, 1)
        name = aud{i,1}; cs = cases.(name); f = results_dir(['pc_' name '_compact.mat']); if ~isfile(f), continue, end
        C = load(f); cpt = C.cpt; len = struct('max_itr', cpt.iters - 1, 'stop_tol', 0, 'stop_dJ_tol', 0);
        if smoke, aud{i,2}.max_itr = cpt.iters - 1; aud{i,2}.es_every = 1; end
        es = cs.mk_cmp(m(cs.sz, len, aud{i,2})); save(results_dir(['pc_' name '_compact_audit.mat']), 'es');
        n = min(numel(es.J), numel(cpt.J));
        L('>>> %s compact audit: %d samples, max state %.2e adjoint %.2e objective %.2e; replay path difference %.1e', name, numel(es.es.k), ...
          max(es.es.err), max(es.es.erra), max(es.es.dJ), max(abs(es.J(1:n) - cpt.J(1:n))./abs(cpt.J(1:n))));
    end
end
L('all phases done');
end

function s = m(varargin)
s = struct();
for i = 1:nargin
    fn = fieldnames(varargin{i});
    for j = 1:numel(fn), s.(fn{j}) = varargin{i}.(fn{j}); end
end
end

function idle(cool, L, what)
if cool > 0, L('--- idle %d s before %s', cool, what); pause(cool); end
end

function report(L, what, r, fom)
n = min(numel(r.J), numel(fom.J)); gap = abs(r.J(1:n) - fom.J(1:n))./abs(fom.J(1:n));
extra = '';
if isfield(r, 'n_rescue'), extra = sprintf(', corrections %d/%d, objective corrections %d/%d', r.n_rescued, r.n_rescue, r.n_goal_corrected, r.n_goal_attempts); end
L('>>> %s: %d its, %.1f s, %d fact, fallbacks %d + verify %d / accepted %d (%.1f%%), speedup %.2fx vs PC FOM | gap end %.2e max %.2e%s', ...
  what, r.iters, r.time, r.nfact, r.n_fom, r.n_verify, r.n_rom, 100*r.n_rom/r.iters, fom.time/r.time, gap(end), max(gap), extra);
end
