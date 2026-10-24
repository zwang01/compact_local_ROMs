% RUN_SEC51_CANNULA  Stokes inflow cannula: fresh-factorization reference,
% local ROM, and same-design error study. Three runs, results in results/:
%   step 1  full-order optimization (Algorithm 1) to the stopping rule,
%           keeping the designs at k = 3000 and 8314 for Figure 2  -> sec51_fom.mat
%           (paper's fresh-factorization comparison: 11181 iterations, 1638 s)
%   step 2  Algorithm 2 for the same number of iterations       -> sec51_rom.mat
%           (historical local-ROM timing; the main table uses the
%            preconditioned-baseline chain instead)
%   step 3  error study: the first 3000 iterations of Algorithm 2, every
%           accepted reduced pair recomputed at full order       -> sec51_errstudy.mat
%           (about 10 min)
% Then fig_sec51 draws the geometry and error-study panels. Run fig_pc_main
% afterward for the final three-method history panels; export_pc_numbers
% supplies Table 1 numbers.
%
%   matlab -batch "setup_paths; run_sec51_cannula"
%   matlab -batch "setup_paths; steps = 2; run_sec51_cannula"     (selected steps)
%
% The runs reported in the paper used the clip-after-update design step
% (vol_project = false, see ac_update_vol.m) and solved the reduced least-squares
% problems by the normal equations (lsq = 'normal', see rom_lsq.m); both are set
% here. Timings are meaningful only on an otherwise idle machine.
if ~exist('steps', 'var'), steps = 1:3; end
common = struct('vol_project', false, 'lsq', 'normal');   % settings of the published runs

if any(steps == 1)
    res = shapeopt_cannula(merge_struct(common, struct('tol_eta', 0, 'save_phi_at', [3000 8314])));
    save(results_dir('sec51_fom.mat'), 'res');
end
if any(steps == 2)
    F = load(results_dir('sec51_fom.mat'));          % matched length: as many iterations as the full-order run
    res = shapeopt_cannula(merge_struct(common, struct('max_itr', F.res.iters - 1, 'stop_tol', 0)));
    save(results_dir('sec51_rom.mat'), 'res');
end
if any(steps == 3)
    res = shapeopt_cannula(merge_struct(common, struct('max_itr', 2999, 'stop_tol', 0, 'errstudy', true)));
    save(results_dir('sec51_errstudy.mat'), 'res');
end
