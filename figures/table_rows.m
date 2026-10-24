function table_rows(label, F, R)
% TABLE_ROWS  The two rows of Tables 1-3 of the paper for one case, and the
% derived numbers quoted in the text, from the results of Algorithm 1 (F) and
% Algorithm 2 (R) run for the same number of iterations.
%   FOM/ROM solves: (full-order iterations: warm-up + fallbacks) + verifications
%                   / accepted reduced iterations
%   fact.         : factorizations of the saddle-point system
nF = numel(F.J); nR = numel(R.J); n = min(nF, nR);
gap = abs(R.J(1:n) - F.J(1:n))./abs(F.J(1:n));
fprintf('\n%s\n', label);
fprintf('  %-13s %6s  %-20s %6s  %-14s %9s  %s\n', 'method', 'iters', 'FOM/ROM solves', 'fact.', 'final J', 'time (s)', '|dJ|/J');
fprintf('  %-13s %6d  %-20s %6d  %-14.7e %9.1f  --\n', 'full order', nF, sprintf('%d / --', F.n_fom), F.nfact, F.J(end), F.time);
fprintf('  %-13s %6d  %-20s %6d  %-14.7e %9.1f  %.1e\n', 'screened ROM', nR, ...
    sprintf('%d+%d / %d', R.n_fom, R.n_verify, R.n_rom), R.nfact, R.J(end), R.time, abs(R.J(end) - F.J(end))/abs(F.J(end)));
Nev = R.n_fom + R.n_verify;
fprintf(['  speedup %.1fx; largest gap along the trajectory %.1e (k = %d); reduced iterations %.1f%%;\n' ...
         '  full-order events %d (%d estimate-triggered verifications), mean GMRES iterations %.1f, %d refreshes;\n' ...
         '  iterations per full-order event N_it/N_ev = %.1f; tau_eta raised %d times%s, lowered %d times\n'], ...
    F.time/R.time, max(gap), find(gap == max(gap), 1) - 1, 100*R.n_rom/nR, Nev, R.n_verify_dwr, ...
    mean(R.pc_its), R.n_refresh, nR/Nev, nnz(diff(R.tol) > 0), raised_at(R.tol), nnz(diff(R.tol) < 0));
end

function s = raised_at(tol)
k = find(diff(tol) > 0);                 % tol(k+1) > tol(k): raised at iteration k
if isempty(k), s = ''; else, s = sprintf(' (at k = %s)', strjoin(arrayfun(@num2str, k(:)', 'UniformOutput', false), ', ')); end
end
