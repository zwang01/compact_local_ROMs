% RUN_SEC52_BODY2D  Isolated-body drag minimization at Re = 40 and 200:
% fresh-factorization references and local ROMs. Results in results/:
%   step 1  full-order optimization (Algorithm 1) to the stopping rule -> sec52_Re<Re>_fom.mat
%           (fresh-factorization comparison: Re = 40, 5738 iterations and
%            9473 s; Re = 200, 4697 iterations and 8210 s)
%   step 2  Algorithm 2 for the same number of iterations              -> sec52_Re<Re>_rom.mat
%           (historical timings; the main table uses the preconditioned chain)
% Then fig_sec52 draws the geometry panels. Run fig_pc_main afterward for
% the final three-method history panels; export_pc_numbers supplies Table 2.
%
%   matlab -batch "setup_paths; run_sec52_body2d"
%   matlab -batch "setup_paths; Res = 40; steps = 2; run_sec52_body2d"
%
% Settings: h = 1/80 on D = (-1,3) x (-1/2,1/2), disk of radius 0.2, solid
% area pi*0.2^2, inflow 2.5(1-4y^2), nu = 1/Re (the defaults of
% shapeopt_body2d). The runs reported in the paper used the clip-after-update
% design step (vol_project = false) and the normal equations for the reduced
% least squares (lsq = 'normal'); both are set here. Algorithm 1 factors the
% converged Jacobian for the adjoint (adj_reuse = false, the default for
% tol_eta = 0).
if ~exist('Res', 'var'), Res = [40 200]; end
if ~exist('steps', 'var'), steps = 1:2; end
common = struct('vol_project', false, 'lsq', 'normal');   % settings of the published runs

for Re = Res
    tag = sprintf('sec52_Re%d', Re);
    if any(steps == 1)
        res = shapeopt_body2d(merge_struct(common, struct('Re', Re, 'tol_eta', 0)));
        save(results_dir([tag '_fom.mat']), 'res');
    end
    if any(steps == 2)
        F = load(results_dir([tag '_fom.mat']));
        res = shapeopt_body2d(merge_struct(common, struct('Re', Re, 'max_itr', F.res.iters - 1, 'stop_tol', 0)));
        save(results_dir([tag '_rom.mat']), 'res');
    end
end
