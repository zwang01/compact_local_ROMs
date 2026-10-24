function N = export_pc_numbers(outfile, fresh_dir)
% EXPORT_PC_NUMBERS  The numbers of the paper's tables and text from the results
% of run_pc_baseline_chain (results_dir()/pc_*.mat), written to one JSON file:
% per case the preconditioned full-order run (fom), the local ROM (rom) and the
% compact local ROM (compact) with times, counts, speedups and objective gaps;
% the Appendix C grids; the compact runs under their own stopping rule; the
% same-design audits. fresh_dir (optional): folder of the fresh-factorization
% full-order runs (sec51_fom.mat ...), whose times are quoted for comparison.
if nargin < 2, fresh_dir = ''; end
cases = {'cannula', 'sec51'; 're40', 'sec52_Re40'; 're200', 'sec52_Re200'; 'body3d', 'sec53'};
N = struct();
for i = 1:size(cases, 1)
    c = cases{i,1}; f = results_dir(['pc_' c '_fom.mat']); if ~isfile(f), continue, end
    F = load(f); F = F.fom;
    N.(c).fom = struct('iters', F.iters, 'time', F.time, 's_per_it', F.time/F.iters, 'nfact', F.nfact, 'n_pc', F.n_pc, ...
                       'n_pc_fail', F.n_pc_fail, 'J_end', F.J(end), 'pc_its_mean', mean(F.pc_its));
    if ~isempty(fresh_dir) && isfile(fullfile(fresh_dir, [cases{i,2} '_fom.mat']))
        G = load(fullfile(fresh_dir, [cases{i,2} '_fom.mat'])); G = G.res;
        N.(c).fresh = struct('iters', G.iters, 'time', G.time, 'nfact', G.nfact, 'J_end', G.J(end), ...
                             'traj_diff', max(abs(G.J(1:min(end,numel(F.J))) - F.J(1:min(end,numel(G.J))))./abs(F.J(1:min(end,numel(G.J))))));
    end
    for v = {'rom', 'compact'}
        g = results_dir(sprintf('pc_%s_%s.mat', c, v{1})); if ~isfile(g), continue, end
        R = load(g); R = R.(fn1(R));
        N.(c).(v{1}) = romstats(R, F);
    end
    g = results_dir(['pc_' c '_compact_ownrule.mat']);
    if isfile(g)
        O = load(g); O = O.own; n = min(numel(O.J), numel(F.J));
        N.(c).ownrule = struct('iters', O.iters, 'gap', abs(O.J(n) - F.J(n))/abs(F.J(n)));
    end
    g = results_dir(['pc_' c '_compact_audit.mat']);
    if isfile(g)
        A = load(g); A = A.es; es = A.es;
        N.(c).audit = struct('samples', numel(es.k), 'state_max', max(es.err), 'state_med', median(es.err), ...
            'adjoint_max', max(es.erra), 'adjoint_med', median(es.erra), 'obj_max', max(es.dJ), 'obj_med', median(es.dJ), ...
            'effC_med', median(es.etaC./es.err), 'effC_min', min(es.etaC./es.err), 'effC_below1', mean(es.etaC < es.err), ...
            'effJ_med', median(es.etaJ./es.dJsigned), 'iters_covered', A.iters);
    end
end
for Ng = [40 80 160]
    f = results_dir(sprintf('pc_appC_N%d_fom.mat', Ng)); if ~isfile(f), continue, end
    F = load(f); F = F.fom; key = sprintf('N%d', Ng);
    N.appC.(key).fom = struct('time', F.time, 'nfact', F.nfact, 'J_end', F.J(end));
    for v = {'rom', 'compact'}
        g = results_dir(sprintf('pc_appC_N%d_%s.mat', Ng, v{1})); if ~isfile(g), continue, end
        R = load(g); R = R.(fn1(R)); N.appC.(key).(v{1}) = romstats(R, F);
    end
end
fid = fopen(outfile, 'w'); fprintf(fid, '%s', jsonencode(N, 'PrettyPrint', true)); fclose(fid);
fprintf('wrote %s\n', outfile);
end

function s = romstats(R, F)
n = min(numel(R.J), numel(F.J)); gap = abs(R.J(1:n) - F.J(1:n))./abs(F.J(1:n)); [gmax, kg] = max(gap);
acc = R.mode(:) == 0; d = acc'; st = diff([0 d 0]); lens = find(st == -1) - find(st == 1);
s = struct('iters', R.iters, 'time', R.time, 'nfact', R.nfact, 'n_fom', R.n_fom, 'n_verify', R.n_verify, 'n_rom', R.n_rom, ...
           'accept', R.n_rom/R.iters, 'speedup', F.time/R.time, 'J_end', R.J(end), 'gap_end', gap(end), 'gap_max', gmax, 'k_gap_max', kg - 1, ...
           'n_pc', R.n_pc, 'n_pc_fail', R.n_pc_fail, 'stretch_max', max(lens), 'stretch_mean', mean(lens), 'n_relax', R.n_relax, ...
           'r_med', median(R.r(R.r > 0)), 'r_max', max(R.r));
if isfield(R, 'n_rescue')
    s.n_rescue = R.n_rescue; s.n_rescued = R.n_rescued; s.n_goal_attempts = R.n_goal_attempts; s.n_goal_corrected = R.n_goal_corrected;
end
end

function f = fn1(S)
f = fieldnames(S); f = f{1};
end
