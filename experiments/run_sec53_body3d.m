% RUN_SEC53_BODY3D  Pironneau's body on the quarter duct, h = 1/40:
% fresh-factorization reference, local ROM, and own-rule run. Results in results/:
%   step 1  full-order optimization (Algorithm 1) to the stopping rule -> sec53_fom.mat
%           (fresh-factorization comparison: 2796 iterations, 35335 s)
%   step 2  Algorithm 2 for the same number of iterations              -> sec53_rom.mat
%           (historical timing; the main table uses the preconditioned chain)
%   step 3  Algorithm 2 under its own stopping rule                    -> sec53_rom_rule.mat
%           (paper: 2581 iterations)
% Then fig_sec53 draws the 3D geometry and flow panels. Run fig_pc_main
% afterward for the final three-method history panel; export_pc_numbers
% supplies Table 3 numbers.
%
%   matlab -batch "setup_paths; run_sec53_body3d"
%   matlab -batch "setup_paths; steps = [2 3]; run_sec53_body3d"
%
% The runs reported in the paper solved the reduced least squares by the normal
% equations (lsq = 'normal', set here). A coarser grid for a quick look:
% N3d = 32.
if ~exist('steps', 'var'), steps = 1:3; end
if ~exist('N3d', 'var'), N3d = 40; end
common = struct('N', N3d, 'lsq', 'normal', 'vol_project', false);   % settings of the published runs

if any(steps == 1)
    res = shapeopt_body3d(merge_struct(common, struct('tol_eta', 0)));
    save(results_dir('sec53_fom.mat'), 'res');
end
if any(steps == 2)
    F = load(results_dir('sec53_fom.mat'));
    res = shapeopt_body3d(merge_struct(common, struct('max_itr', F.res.iters - 1, 'stop_tol', 0)));
    save(results_dir('sec53_rom.mat'), 'res');
end
if any(steps == 3)
    res = shapeopt_body3d(common);
    save(results_dir('sec53_rom_rule.mat'), 'res');
end
