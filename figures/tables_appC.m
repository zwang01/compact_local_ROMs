% TABLES_APPC  Tables C.2-C.4 of the paper (Appendix C) from the results of
% run_appC_validation (isolated body, Re = 40, 300-iteration budget).
%   Table C.2: coupled sensitivity study over h = 1/40, 1/80, 1/160
%   Table C.3: ablation at h = 1/80
%   Table C.4: indicators against true errors on the sampled accepted iterations
%   matlab -batch "setup_paths; tables_appC"
ld = @(lab) getfield(load(results_dir(['appC_' lab '.mat'])), 'res'); %#ok<GFLD>
gap = @(R, F) abs(R.J - F.J)./abs(F.J);           % along the trajectory (all runs have 301 iterations)

% ---------------- Table C.2 ----------------
fprintf('Table C.2 (300 iterations; times in s: screened ROM, FOM with a fresh factorization, FOM through the retained factors)\n');
fprintf('%7s %8s %-14s %5s %10s %9s %9s %7s %7s %7s\n', 'h', 'dB/h', 'FOM/ROM', 'fact', 'final D', '|dD|/D', 'max', 'ROM', 'FOM', 'FOMret');
for N = [40 80 160]
    f = results_dir(sprintf('appC_refine_rom_N%d.mat', N));
    if ~isfile(f), continue, end
    R = ld(sprintf('refine_rom_N%d', N)); F = ld(sprintf('refine_fom_N%d', N)); K = ld(sprintf('refine_fomkry_N%d', N));
    g = gap(R, F); dB = sqrt(R.par.nu/R.par.alpha_max)/R.par.h;    % Brinkman penetration depth over h
    fprintf('%7s %8.2f %-14s %5d %10.6f %9.1e %9.1e %7.0f %7.0f %7.0f\n', sprintf('1/%d', N), dB, ...
        sprintf('%d+%d / %d', R.n_fom, R.n_verify, R.n_rom), R.nfact, R.J(end), g(end), max(g), R.time, F.time, K.time);
end

% ---------------- Table C.3 ----------------
if isfile(results_dir('appC_refine_fom_N80.mat'))
    F = ld('refine_fom_N80');
    rows = {'refine_rom_N80', 'Algorithm 2 as stated'; 'ablate_eucl', 'Euclidean gate at 3e-3'; ...
            'ablate_backstop', 'backstop verification only'; 'ablate_refactor', 'fresh factorization at every fallback'; ...
            'ablate_window80', 'window 2 w_max'; 'ablate_cumul', 'cumulative history'};
    fprintf('\nTable C.3 (h = 1/80, 300 iterations; fallbacks include the warm-up)\n');
    fprintf('%-38s %9s %6s %5s %5s %11s %9s %9s %8s\n', 'variant', 'fallbacks', 'verif.', 'ROM', 'fact', 'final D', '|dD|/D', 'max', 'time (s)');
    for i = 1:size(rows, 1)
        if ~isfile(results_dir(['appC_' rows{i,1} '.mat'])), continue, end
        R = ld(rows{i,1}); g = gap(R, F);
        fprintf('%-38s %9d %6d %5d %5d %11.7f %9.1e %9.1e %8.0f\n', rows{i,2}, R.n_fom, R.n_verify, R.n_rom, ...
            R.nfact, R.J(end), g(end), max(g), R.time);
    end
end

% ---------------- Table C.4 ----------------
if isfile(results_dir('appC_validate_es5.mat'))
    es = ld('validate_es5'); es = es.es;
    rat = {'Euclidean eta^E_s',        'state error',     es.etaE ./es.err; ...
           'preconditioned eta_s',     'state error',     es.etaP ./es.err; ...
           'calibrated eta~_s',        'state error',     es.etaC ./es.err; ...
           'preconditioned eta_a',     'adjoint error',   es.etaPa./es.erra; ...
           'calibrated eta~_a',        'adjoint error',   es.etaCa./es.erra; ...
           'dual-weighted |eta_J|',    'objective error', abs(es.etaJ./es.dJsigned)};
    fprintf('\nTable C.4 (%d sampled accepted iterations; ratio indicator / true error)\n', numel(es.k));
    fprintf('%-26s %-16s %8s %8s %8s\n', 'indicator', 'against', 'min', 'median', 'max');
    for i = 1:size(rat, 1)
        r = rat{i,3};
        fprintf('%-26s %-16s %8.3g %8.3g %8.3g\n', rat{i,1}, rat{i,2}, min(r), median(r), max(r));
    end
    fprintf(['true relative errors: state median %.1e, max %.1e; adjoint max %.1e; objective max %.1e;\n' ...
             'calibrated state indicator below the true error at %d of %d samples\n'], ...
        median(es.err), max(es.err), max(es.erra), max(es.dJ), nnz(es.etaC < es.err), numel(es.k));
end
